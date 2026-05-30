#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
sonar2_calc_singlefile_v1_0.py

Single-file SONAR-Calc batch engine for Protein Acoustome Atlas / SONAR2.

Purpose
-------
This is a clean, distributable single-program version of the current H-AEND / SONAR-Calc
workflow. It does not import local project scripts and does not call nested pipeline files.
It directly performs, in one file:

1. mmCIF/PDB structure loading, including .cif.gz
2. residue-node construction for one selected chain
3. chemically annotated residue-contact network construction
4. reversible H-AEND acoustic response simulation
5. destructive H-AEND epsilon sweep
6. compact SONAR label table generation
7. multiprocessing, checkpoint-like part outputs, resume/skip support

Main input
----------
A chain table CSV, e.g.:
    E:\\03_Research_Data\\pdb\\sonar2_stage1A\\outputs\\sonar2_stage1A_selected_chains.csv

Required input columns:
    pdb_id, file_path, chain_id

Typical run
-----------
    python sonar2_calc_singlefile_v1_0.py --params E:\\03_Research_Data\\pdb\\sonar2_stage1B\\params_sonar2_calc_v1_0.csv

Dependencies
------------
    pip install numpy pandas biopython

Notes
-----
- This program is intentionally table-first. It does not generate figures.
- For very large screening, use run_destructive=0 first if speed is limiting.
- For final paper labels, keep one frozen params CSV and report it as the SONAR-Calc v1.0 configuration.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import shutil
import sys
import tempfile
import time
import traceback
from collections import Counter, deque
from dataclasses import dataclass
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    from Bio.PDB import MMCIFParser, PDBParser
except ImportError as exc:
    raise ImportError(
        "Biopython is required. Install it with:\n"
        "  pip install biopython\n"
    ) from exc

# =============================================================================
# Chemistry definitions
# =============================================================================

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
}
STANDARD_AA = set(AA3_TO_1.keys())

HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}
AROMATIC = {"PHE", "TYR", "TRP", "HIS"}
POLAR_OR_CHARGED = {"SER", "THR", "ASN", "GLN", "TYR", "CYS", "ASP", "GLU", "LYS", "ARG", "HIS"}
POSITIVE = {"LYS", "ARG", "HIS"}
NEGATIVE = {"ASP", "GLU"}
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}
POLAR_ATOM_PREFIX = ("N", "O", "S")

CHARGE_ATOMS = {
    "LYS": {"NZ"},
    "ARG": {"NH1", "NH2", "NE"},
    "HIS": {"ND1", "NE2"},
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
}

RING_ATOMS = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TRP": {"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "HIS": {"CG", "ND1", "CD2", "CE1", "NE2"},
}

CONTACT_ORDER = [
    "backbone",
    "disulfide",
    "salt_bridge",
    "aromatic",
    "hydrophobic",
    "polar_hbond_candidate",
    "generic_weak",
]

# =============================================================================
# General utilities
# =============================================================================

def read_params(path: str) -> Dict[str, str]:
    df = pd.read_csv(path)
    if {"key", "value"}.issubset(df.columns):
        return {str(k).strip(): str(v).strip() for k, v in zip(df["key"], df["value"]) if str(k).strip()}
    if {"param", "value"}.issubset(df.columns):
        return {str(k).strip(): str(v).strip() for k, v in zip(df["param"], df["value"]) if str(k).strip()}
    raise ValueError("Parameter CSV must contain key,value or param,value columns.")


def p_get(params: Dict[str, str], key: str, default=None, cast=str):
    val = params.get(key, default)
    if val is None:
        return None
    if cast is bool:
        return str(val).strip().lower() in {"1", "true", "yes", "y", "on"}
    try:
        return cast(val)
    except Exception:
        return default


def ensure_dir(path: str | Path) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def bool_int(x) -> int:
    return 1 if bool(x) else 0


def safe_float(x, default=np.nan) -> float:
    try:
        if pd.isna(x):
            return default
        return float(x)
    except Exception:
        return default


def safe_int(x, default=0) -> int:
    try:
        if pd.isna(x):
            return default
        return int(float(x))
    except Exception:
        return default


def dist(a, b) -> float:
    if a is None or b is None:
        return np.nan
    return float(np.linalg.norm(np.asarray(a, dtype=float) - np.asarray(b, dtype=float)))


def safe_min_distance(coords_a, coords_b) -> float:
    if len(coords_a) == 0 or len(coords_b) == 0:
        return np.nan
    arr_a = np.asarray(coords_a, dtype=float)
    arr_b = np.asarray(coords_b, dtype=float)
    dmin = np.inf
    for x in arr_a:
        dd = np.linalg.norm(arr_b - x, axis=1)
        m = float(np.min(dd))
        if m < dmin:
            dmin = m
    return dmin if np.isfinite(dmin) else np.nan


def centroid(coords):
    if len(coords) == 0:
        return (np.nan, np.nan, np.nan)
    arr = np.asarray(coords, dtype=float)
    c = np.mean(arr, axis=0)
    return tuple(float(x) for x in c)


def remove_com(r, v=None):
    r2 = r - np.mean(r, axis=0, keepdims=True)
    if v is None:
        return r2
    v2 = v - np.mean(v, axis=0, keepdims=True)
    return r2, v2


def minmax_normalize_series(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    if x.notna().sum() == 0:
        return pd.Series(np.nan, index=s.index)
    mn = x.min()
    mx = x.max()
    if not np.isfinite(mn) or not np.isfinite(mx):
        return pd.Series(np.nan, index=s.index)
    if abs(mx - mn) < 1e-15:
        out = pd.Series(0.5, index=s.index)
        out[x.isna()] = np.nan
        return out
    return (x - mn) / (mx - mn)

# =============================================================================
# Structure parsing and node construction
# =============================================================================

def atom_coord(residue, atom_name):
    if atom_name not in residue:
        return None
    atom = residue[atom_name]
    return tuple(float(x) for x in atom.get_coord())


def get_atom_records(residue):
    records = []
    for atom in residue.get_atoms():
        name = atom.get_name().strip()
        element = (atom.element or "").strip().upper()
        if element == "H" or name.startswith("H"):
            continue
        records.append({
            "name": name,
            "element": element,
            "coord": tuple(float(x) for x in atom.get_coord()),
        })
    return records


def is_standard_residue(residue):
    hetfield, resseq, icode = residue.get_id()
    return hetfield == " " and residue.get_resname().strip().upper() in STANDARD_AA


def residue_sort_key(residue):
    hetfield, resseq, icode = residue.get_id()
    return int(resseq), str(icode).strip()


def load_structure_any(input_file: str, protein_id: str):
    """Load .cif, .mmcif, .pdb, or .cif.gz/.pdb.gz using Biopython."""
    input_file = str(input_file)
    lower = input_file.lower()
    tmp_path = None
    parse_path = input_file

    if lower.endswith(".gz"):
        suffix = ".cif" if lower.endswith(".cif.gz") or lower.endswith(".mmcif.gz") else ".pdb"
        fd, tmp_path = tempfile.mkstemp(prefix="sonar2_", suffix=suffix)
        os.close(fd)
        with gzip.open(input_file, "rb") as fin, open(tmp_path, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        parse_path = tmp_path
        lower = lower[:-3]

    try:
        if lower.endswith(".cif") or lower.endswith(".mmcif"):
            parser = MMCIFParser(QUIET=True)
        elif lower.endswith(".pdb") or lower.endswith(".ent"):
            parser = PDBParser(QUIET=True)
        else:
            raise ValueError(f"Unsupported structure format: {input_file}")
        return parser.get_structure(protein_id, parse_path)
    finally:
        if tmp_path and os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except Exception:
                pass


def find_chain(model, chain_id: str):
    chains = list(model.get_chains())
    for c in chains:
        if str(c.id) == str(chain_id):
            return c
    # fallback for occasional whitespace differences
    for c in chains:
        if str(c.id).strip() == str(chain_id).strip():
            return c
    available = [str(c.id) for c in chains]
    raise ValueError(f"Chain {chain_id} not found. Available chains: {available[:30]}")


def build_nodes_from_structure(structure, protein_id: str, chain_id: str, params: Dict[str, str]):
    model_index = p_get(params, "model_index", 0, int)
    models = list(structure.get_models())
    if model_index >= len(models):
        raise ValueError(f"model_index={model_index} out of range; structure has {len(models)} models")
    model = models[model_index]
    chain = find_chain(model, chain_id)

    residues = [r for r in chain.get_residues() if is_standard_residue(r)]
    residues = sorted(residues, key=residue_sort_key)

    node_rows = []
    payloads = []

    for idx, res in enumerate(residues):
        resname = res.get_resname().strip().upper()
        one = AA3_TO_1.get(resname, "X")
        hetfield, resseq, icode = res.get_id()
        icode_clean = str(icode).strip()

        ca = atom_coord(res, "CA")
        cb = atom_coord(res, "CB")
        is_gly = resname == "GLY"
        if is_gly:
            node_atom = "CA"
            node = ca
        else:
            node_atom = "CB" if cb is not None else "CA"
            node = cb if cb is not None else ca
        if node is None:
            continue

        atom_records = get_atom_records(res)
        heavy_coords = [a["coord"] for a in atom_records]
        sidechain_records = [a for a in atom_records if a["name"] not in BACKBONE_ATOMS]
        sidechain_coords = [a["coord"] for a in sidechain_records]

        sc_cent = centroid(sidechain_coords)
        ring_names = RING_ATOMS.get(resname, set())
        ring_coords = [a["coord"] for a in atom_records if a["name"] in ring_names]
        ring_cent = centroid(ring_coords)

        polar_coords = [
            a["coord"] for a in atom_records
            if a["element"].startswith(POLAR_ATOM_PREFIX) or a["name"].startswith(POLAR_ATOM_PREFIX)
        ]
        charge_names = CHARGE_ATOMS.get(resname, set())
        charge_coords = [a["coord"] for a in atom_records if a["name"] in charge_names]
        sg_coord = atom_coord(res, "SG") if resname == "CYS" else None

        node_x, node_y, node_z = node
        ca_x, ca_y, ca_z = ca if ca is not None else (np.nan, np.nan, np.nan)
        cb_x, cb_y, cb_z = cb if cb is not None else (np.nan, np.nan, np.nan)

        node_rows.append({
            "protein_id": protein_id,
            "chain_id": chain_id,
            "node_index": idx,
            "res_seq": int(resseq),
            "icode": icode_clean,
            "res_name": resname,
            "one_letter": one,
            "node_atom": node_atom,
            "node_x": node_x,
            "node_y": node_y,
            "node_z": node_z,
            "ca_x": ca_x,
            "ca_y": ca_y,
            "ca_z": ca_z,
            "cb_x": cb_x,
            "cb_y": cb_y,
            "cb_z": cb_z,
            "sidechain_centroid_x": sc_cent[0],
            "sidechain_centroid_y": sc_cent[1],
            "sidechain_centroid_z": sc_cent[2],
            "has_sidechain": bool_int(len(sidechain_coords) > 0),
            "has_missing_key_atom": 0,
            "is_gly": bool_int(is_gly),
            "is_cys": bool_int(resname == "CYS"),
            "is_aromatic": bool_int(resname in AROMATIC),
            "is_hydrophobic": bool_int(resname in HYDROPHOBIC),
            "is_polar": bool_int(resname in POLAR_OR_CHARGED),
            "is_charged_pos": bool_int(resname in POSITIVE),
            "is_charged_neg": bool_int(resname in NEGATIVE),
        })

        payloads.append({
            "node_index": idx,
            "residue": res,
            "chain_id": chain_id,
            "res_seq": int(resseq),
            "icode": icode_clean,
            "res_name": resname,
            "one_letter": one,
            "node_coord": node,
            "ca_coord": ca,
            "cb_coord": cb,
            "heavy_coords": heavy_coords,
            "sidechain_coords": sidechain_coords,
            "polar_coords": polar_coords,
            "charge_coords": charge_coords,
            "ring_centroid": ring_cent,
            "ring_has": not any(np.isnan(ring_cent)),
            "ring_coords": ring_coords,
            "sg_coord": sg_coord,
        })

    return pd.DataFrame(node_rows), payloads

# =============================================================================
# Network construction
# =============================================================================

def type_params(params: Dict[str, str], prefix: str):
    k0 = p_get(params, f"{prefix}_k0", 1.0, float)
    c_ac = p_get(params, f"{prefix}_c_ac", 0.5, float)
    phi_pi = p_get(params, f"{prefix}_phi_pi", 0.0, float)
    return {"k0": k0, "c_ac": c_ac, "phi_pi": phi_pi, "phi_rad": phi_pi * math.pi}


def make_param_map(params: Dict[str, str]):
    return {t: type_params(params, t) for t in CONTACT_ORDER}


def choose_primary_type(flags: Dict[str, bool]) -> str:
    for t in CONTACT_ORDER:
        if flags.get(f"is_{t}", False):
            return t
    return ""


def detect_pair_features(ri, rj):
    seq_sep = abs(int(ri["res_seq"]) - int(rj["res_seq"]))
    node_distance = dist(ri["node_coord"], rj["node_coord"])
    min_heavy = safe_min_distance(ri["heavy_coords"], rj["heavy_coords"])
    min_sc = safe_min_distance(ri["sidechain_coords"], rj["sidechain_coords"])
    min_polar = safe_min_distance(ri["polar_coords"], rj["polar_coords"])
    min_charge = safe_min_distance(ri["charge_coords"], rj["charge_coords"])
    ring_d = dist(ri["ring_centroid"], rj["ring_centroid"]) if ri["ring_has"] and rj["ring_has"] else np.nan
    min_ring_atom_d = safe_min_distance(ri.get("ring_coords", []), rj.get("ring_coords", []))
    sg_d = dist(ri["sg_coord"], rj["sg_coord"]) if ri["res_name"] == "CYS" and rj["res_name"] == "CYS" else np.nan
    return {
        "seq_sep": seq_sep,
        "node_distance": node_distance,
        "min_heavy_distance": min_heavy,
        "min_sidechain_distance": min_sc,
        "min_polar_distance": min_polar,
        "min_charge_distance": min_charge,
        "ring_center_distance": ring_d,
        "min_ring_atom_distance": min_ring_atom_d,
        "sg_sg_distance": sg_d,
    }


def build_edges_v1_1(nodes_df: pd.DataFrame, payloads: List[dict], protein_id: str, chain_id: str, params: Dict[str, str]) -> pd.DataFrame:
    include_backbone = p_get(params, "include_backbone_edges", 1, bool)
    include_disulfide = p_get(params, "include_disulfide_edges", 1, bool)
    include_nonlocal = p_get(params, "include_nonlocal_edges", 1, bool)
    seq_excl = p_get(params, "sequence_exclusion_nonlocal", 2, int)

    disulfide_cutoff = p_get(params, "disulfide_cutoff_A", 2.5, float)
    salt_cutoff = p_get(params, "saltbridge_cutoff_A", 4.0, float)
    aromatic_cutoff = p_get(params, "aromatic_cutoff_A", 5.5, float)
    aromatic_atom_cutoff = p_get(params, "aromatic_atom_cutoff_A", 4.8, float)
    hydrophobic_cutoff = p_get(params, "hydrophobic_cutoff_A", 5.0, float)
    polar_cutoff = p_get(params, "polar_cutoff_A", 3.5, float)
    generic_cutoff = p_get(params, "generic_heavy_cutoff_A", 4.5, float)
    node_prefilter = p_get(params, "node_distance_prefilter_A", 14.0, float)

    param_map = make_param_map(params)
    edge_rows = []
    n_payloads = len(payloads)

    for a in range(n_payloads - 1):
        ri = payloads[a]
        for b in range(a + 1, n_payloads):
            rj = payloads[b]
            i = int(ri["node_index"])
            j = int(rj["node_index"])
            seq_sep = abs(int(ri["res_seq"]) - int(rj["res_seq"]))

            # Fast prefilter for non-adjacent/non-disulfide pairs.
            node_d_quick = dist(ri["node_coord"], rj["node_coord"])
            if seq_sep != 1 and node_d_quick > node_prefilter:
                if not (ri["res_name"] == "CYS" and rj["res_name"] == "CYS" and node_d_quick < 20.0):
                    continue

            aa_i = ri["res_name"]
            aa_j = rj["res_name"]
            feat = detect_pair_features(ri, rj)
            seq_sep = feat["seq_sep"]
            flags = {f"is_{t}": False for t in CONTACT_ORDER}

            if include_backbone and seq_sep == 1:
                flags["is_backbone"] = True

            if include_disulfide and aa_i == "CYS" and aa_j == "CYS":
                sg_d = feat["sg_sg_distance"]
                if not np.isnan(sg_d) and sg_d < disulfide_cutoff:
                    flags["is_disulfide"] = True

            if include_nonlocal and seq_sep > seq_excl:
                if ((aa_i in POSITIVE and aa_j in NEGATIVE) or (aa_i in NEGATIVE and aa_j in POSITIVE)):
                    d = feat["min_charge_distance"]
                    if not np.isnan(d) and d < salt_cutoff:
                        flags["is_salt_bridge"] = True

                if aa_i in AROMATIC and aa_j in AROMATIC:
                    rc = feat["ring_center_distance"]
                    ra = feat["min_ring_atom_distance"]
                    if (not np.isnan(rc) and rc < aromatic_cutoff) or (not np.isnan(ra) and ra < aromatic_atom_cutoff):
                        flags["is_aromatic"] = True

                if aa_i in HYDROPHOBIC and aa_j in HYDROPHOBIC:
                    d = feat["min_sidechain_distance"]
                    if not np.isnan(d) and d < hydrophobic_cutoff:
                        flags["is_hydrophobic"] = True

                d = feat["min_polar_distance"]
                if not np.isnan(d) and d < polar_cutoff:
                    flags["is_polar_hbond_candidate"] = True

                d = feat["min_heavy_distance"]
                if not np.isnan(d) and d < generic_cutoff:
                    flags["is_generic_weak"] = True

            if not any(flags.values()):
                continue

            contact_type = choose_primary_type(flags)
            p = param_map[contact_type]
            evidence_score = int(sum(bool(v) for v in flags.values()))
            d0 = feat["node_distance"]
            if not np.isfinite(d0) or d0 <= 0:
                continue

            edge_rows.append({
                "protein_id": protein_id,
                "chain_i": chain_id,
                "i_node_index": i,
                "res_i": ri["res_seq"],
                "icode_i": ri["icode"],
                "aa_i": aa_i,
                "chain_j": chain_id,
                "j_node_index": j,
                "res_j": rj["res_seq"],
                "icode_j": rj["icode"],
                "aa_j": aa_j,
                "seq_sep": seq_sep,
                "node_distance": feat["node_distance"],
                "min_heavy_distance": feat["min_heavy_distance"],
                "min_sidechain_distance": feat["min_sidechain_distance"],
                "min_polar_distance": feat["min_polar_distance"],
                "min_charge_distance": feat["min_charge_distance"],
                "ring_center_distance": feat["ring_center_distance"],
                "min_ring_atom_distance": feat["min_ring_atom_distance"],
                "sg_sg_distance": feat["sg_sg_distance"],
                "contact_type": contact_type,
                "is_backbone": bool_int(flags["is_backbone"]),
                "is_disulfide": bool_int(flags["is_disulfide"]),
                "is_salt_bridge": bool_int(flags["is_salt_bridge"]),
                "is_aromatic": bool_int(flags["is_aromatic"]),
                "is_hydrophobic": bool_int(flags["is_hydrophobic"]),
                "is_polar_hbond_candidate": bool_int(flags["is_polar_hbond_candidate"]),
                "is_generic_weak": bool_int(flags["is_generic_weak"]),
                "evidence_score": evidence_score,
                "k0": p["k0"],
                "c_ac": p["c_ac"],
                "phi_rad": p["phi_rad"],
                "phi_pi": p["phi_pi"],
                "d0": d0,
            })

    return pd.DataFrame(edge_rows)


def build_network_summary(nodes_df: pd.DataFrame, edges_df: pd.DataFrame) -> Dict[str, float]:
    n_nodes = len(nodes_df)
    n_edges = len(edges_df)
    degree = Counter()
    if n_edges:
        for _, row in edges_df.iterrows():
            degree[int(row["i_node_index"])] += 1
            degree[int(row["j_node_index"])] += 1
    degree_values = np.array([degree[i] for i in range(n_nodes)], dtype=float) if n_nodes else np.array([])
    out = {
        "n_nodes": n_nodes,
        "n_edges": n_edges,
        "mean_degree": float(np.mean(degree_values)) if n_nodes else np.nan,
        "median_degree": float(np.median(degree_values)) if n_nodes else np.nan,
        "max_degree": int(np.max(degree_values)) if n_nodes else 0,
        "n_isolated_nodes": int(np.sum(degree_values == 0)) if n_nodes else 0,
        "contact_density_per_residue": float(n_edges / n_nodes) if n_nodes else np.nan,
        "graph_edge_density": float(2.0 * n_edges / (n_nodes * (n_nodes - 1))) if n_nodes > 1 else np.nan,
    }
    for t in CONTACT_ORDER:
        n = int((edges_df["contact_type"] == t).sum()) if n_edges else 0
        out[f"n_{t}_edges"] = n
        out[f"{t}_fraction"] = float(n / n_edges) if n_edges else 0.0
    return out

# =============================================================================
# Reversible H-AEND
# =============================================================================

def prepare_edge_data_reversible(nodes_df: pd.DataFrame, edges_df: pd.DataFrame):
    r0 = nodes_df[["node_x", "node_y", "node_z"]].values.astype(float)
    if np.isnan(r0).any():
        raise ValueError("NaN found in node coordinates.")
    i_idx = edges_df["i_node_index"].values.astype(int)
    j_idx = edges_df["j_node_index"].values.astype(int)
    k0 = edges_df["k0"].values.astype(float)
    c_ac = edges_df["c_ac"].values.astype(float)
    phi = edges_df["phi_rad"].values.astype(float)
    d0 = edges_df["d0"].values.astype(float)
    if np.any(d0 <= 0) or np.isnan(d0).any():
        d0 = np.linalg.norm(r0[i_idx] - r0[j_idx], axis=1)
    return r0, {"i_idx": i_idx, "j_idx": j_idx, "k0": k0, "c_ac": c_ac, "phi": phi, "d0": d0, "contact_type": edges_df["contact_type"].astype(str).values}


def compute_forces_reversible(r, phase, edge_data, epsilon0):
    F = np.zeros_like(r)
    i_idx = edge_data["i_idx"]
    j_idx = edge_data["j_idx"]
    k0 = edge_data["k0"]
    c_ac = edge_data["c_ac"]
    phi = edge_data["phi"]
    d0 = edge_data["d0"]

    ri = r[i_idx]
    rj = r[j_idx]
    rij = ri - rj
    d = np.linalg.norm(rij, axis=1)
    d_safe = np.where(d < 1e-12, 1e-12, d)
    d0_t = d0 * (1.0 - c_ac * epsilon0 * np.sin(phase + phi))
    extension = d - d0_t
    scalar = -k0 * extension / d_safe
    fij = rij * scalar[:, None]
    np.add.at(F, i_idx, fij)
    np.add.at(F, j_idx, -fij)
    rel_static = np.abs(d - d0) / d0
    rel_driven = np.abs(d - d0_t) / d0
    energy_edge = 0.5 * k0 * extension * extension
    return F, float(np.sum(energy_edge)), energy_edge, rel_static, rel_driven


def simulate_reversible(nodes_df: pd.DataFrame, edges_df: pd.DataFrame, params: Dict[str, str]) -> Dict[str, object]:
    r0, edge_data = prepare_edge_data_reversible(nodes_df, edges_df)
    epsilon0 = p_get(params, "epsilon0", 0.015, float)
    n_cycles = p_get(params, "reversible_n_cycles", p_get(params, "n_cycles", 8, int), int)
    steps_per_cycle = p_get(params, "reversible_steps_per_cycle", p_get(params, "steps_per_cycle", 100, int), int)
    dt = p_get(params, "dt_dimensionless", 0.0025, float)
    mass = p_get(params, "mass", 1.0, float)
    gamma = p_get(params, "gamma", 2.0, float)
    contact_break_threshold = p_get(params, "contact_break_threshold", 0.08, float)
    remove_center = p_get(params, "remove_center_of_mass", 1, bool)

    total_steps = int(n_cycles * steps_per_cycle)
    n_edges = len(edge_data["i_idx"])
    n_nodes = r0.shape[0]
    r = r0.copy().astype(float)
    r = r - np.mean(r, axis=0, keepdims=True)
    r_native = r.copy()
    v = np.zeros_like(r)
    start_avg = total_steps // 2

    edge_rel_sum = np.zeros(n_edges, dtype=float)
    edge_rel_max = np.zeros(n_edges, dtype=float)
    edge_count = 0
    residue_resp_sum = np.zeros(n_nodes, dtype=float)
    residue_count = np.zeros(n_nodes, dtype=float)

    mean_rel_records = []
    rmsd_records = []
    broken_records = []
    energy_records = []

    for step in range(total_steps):
        phase = 2.0 * math.pi * (step % steps_per_cycle) / steps_per_cycle
        F, energy, energy_edge, rel_static, rel_driven = compute_forces_reversible(r, phase, edge_data, epsilon0)
        a = (F - gamma * v) / mass
        v = v + dt * a
        r = r + dt * v
        if remove_center:
            r, v = remove_com(r, v)
        F2, energy2, energy_edge2, rel_static2, rel_driven2 = compute_forces_reversible(r, phase, edge_data, epsilon0)
        rmsd = float(np.sqrt(np.mean(np.sum((r - r_native) ** 2, axis=1))))
        mean_rel = float(np.mean(rel_static2)) if n_edges else np.nan
        broken = float(np.mean(rel_static2 > contact_break_threshold)) if n_edges else np.nan
        mean_rel_records.append(mean_rel)
        rmsd_records.append(rmsd)
        broken_records.append(broken)
        energy_records.append(energy2)
        if step >= start_avg:
            edge_rel_sum += rel_static2
            edge_rel_max = np.maximum(edge_rel_max, rel_static2)
            edge_count += 1
            np.add.at(residue_resp_sum, edge_data["i_idx"], rel_static2)
            np.add.at(residue_resp_sum, edge_data["j_idx"], rel_static2)
            np.add.at(residue_count, edge_data["i_idx"], 1.0)
            np.add.at(residue_count, edge_data["j_idx"], 1.0)

    denom = max(1, edge_count)
    edge_avg = edge_rel_sum / denom
    residue_avg = residue_resp_sum / np.maximum(residue_count, 1.0)
    half_mean = float(np.nanmean(mean_rel_records[start_avg:])) if len(mean_rel_records) else np.nan

    contact_type = edge_data["contact_type"]
    type_response = {}
    for t in CONTACT_ORDER:
        mask = contact_type == t
        type_response[f"{t}_reversible_response"] = float(np.nanmean(edge_avg[mask])) if np.any(mask) else np.nan

    top_idx = int(np.nanargmax(residue_avg)) if len(residue_avg) else -1
    if top_idx >= 0:
        top_row = nodes_df.iloc[top_idx]
        top_hotspot_residue = f"{str(top_row['res_name']).upper()}{int(top_row['res_seq'])}"
        top_hotspot_score = float(residue_avg[top_idx])
    else:
        top_hotspot_residue = ""
        top_hotspot_score = np.nan

    out = {
        "epsilon_reversible": epsilon0,
        "mean_reversible_edge_deformation": half_mean,
        "final_mean_reversible_edge_deformation": float(mean_rel_records[-1]) if mean_rel_records else np.nan,
        "max_rmsd_reversible": float(np.nanmax(rmsd_records)) if rmsd_records else np.nan,
        "max_broken_fraction_proxy_reversible": float(np.nanmax(broken_records)) if broken_records else np.nan,
        "mean_elastic_energy_reversible": float(np.nanmean(energy_records[start_avg:])) if energy_records else np.nan,
        "max_elastic_energy_reversible": float(np.nanmax(energy_records)) if energy_records else np.nan,
        "top_hotspot_residue": top_hotspot_residue,
        "top_hotspot_score": top_hotspot_score,
    }
    out.update(type_response)
    return out

# =============================================================================
# Destructive H-AEND
# =============================================================================

def get_type_damage_params(params: Dict[str, str]):
    default_soft = {
        "backbone": 0.20,
        "disulfide": 0.16,
        "salt_bridge": 0.035,
        "aromatic": 0.040,
        "hydrophobic": 0.045,
        "polar_hbond_candidate": 0.030,
        "generic_weak": 0.025,
    }
    default_break = {
        "backbone": 0.50,
        "disulfide": 0.35,
        "salt_bridge": 0.090,
        "aromatic": 0.100,
        "hydrophobic": 0.110,
        "polar_hbond_candidate": 0.070,
        "generic_weak": 0.055,
    }
    soft = {t: p_get(params, f"{t}_softening_delta", default_soft[t], float) for t in CONTACT_ORDER}
    brk = {t: p_get(params, f"{t}_break_threshold", default_break[t], float) for t in CONTACT_ORDER}
    return soft, brk


def largest_connected_component_fraction(n_nodes, i_idx, j_idx, active_mask):
    adj = [[] for _ in range(n_nodes)]
    for e, active in enumerate(active_mask):
        if not active:
            continue
        i = int(i_idx[e]); j = int(j_idx[e])
        adj[i].append(j); adj[j].append(i)
    visited = np.zeros(n_nodes, dtype=bool)
    largest = 0
    for start in range(n_nodes):
        if visited[start]:
            continue
        q = deque([start])
        visited[start] = True
        size = 0
        while q:
            u = q.popleft(); size += 1
            for v in adj[u]:
                if not visited[v]:
                    visited[v] = True; q.append(v)
        largest = max(largest, size)
    return largest / max(1, n_nodes)


def prepare_edge_data_damage(nodes_df: pd.DataFrame, edges_df: pd.DataFrame):
    r0, rev = prepare_edge_data_reversible(nodes_df, edges_df)
    return r0, {
        "i_idx": rev["i_idx"],
        "j_idx": rev["j_idx"],
        "k0_base": rev["k0"],
        "c_ac": rev["c_ac"],
        "phi": rev["phi"],
        "d0": rev["d0"],
        "contact_type": rev["contact_type"],
    }


def compute_forces_damage(r, phase, edge_data, epsilon0, active_mask, softening_delta, use_softening=True):
    F = np.zeros_like(r)
    i_idx = edge_data["i_idx"]
    j_idx = edge_data["j_idx"]
    k0_base = edge_data["k0_base"]
    c_ac = edge_data["c_ac"]
    phi = edge_data["phi"]
    d0 = edge_data["d0"]
    ri = r[i_idx]
    rj = r[j_idx]
    rij = ri - rj
    d = np.linalg.norm(rij, axis=1)
    d_safe = np.where(d < 1e-12, 1e-12, d)
    d0_t = d0 * (1.0 - c_ac * epsilon0 * np.sin(phase + phi))
    rel_static = np.abs(d - d0) / d0
    if use_softening:
        delta = np.maximum(softening_delta, 1e-8)
        soft_factor = np.exp(-((rel_static / delta) ** 2))
    else:
        soft_factor = np.ones_like(rel_static)
    k_eff = k0_base * soft_factor * active_mask.astype(float)
    extension = d - d0_t
    scalar = -k_eff * extension / d_safe
    fij = rij * scalar[:, None]
    np.add.at(F, i_idx, fij)
    np.add.at(F, j_idx, -fij)
    energy_edge = 0.5 * k_eff * extension * extension
    return F, float(np.sum(energy_edge)), rel_static


def simulate_destructive_single_epsilon(nodes_df, edges_df, r0, edge_data, params, epsilon0: float) -> Dict[str, object]:
    n_cycles = p_get(params, "destructive_n_cycles", 8, int)
    steps_per_cycle = p_get(params, "destructive_steps_per_cycle", 80, int)
    dt = p_get(params, "dt_dimensionless", 0.0025, float)
    mass = p_get(params, "mass", 1.0, float)
    gamma = p_get(params, "gamma", 2.0, float)
    remove_center = p_get(params, "remove_center_of_mass", 1, bool)
    use_softening = p_get(params, "use_softening", 1, bool)
    use_breaking = p_get(params, "use_irreversible_breaking", 1, bool)
    break_hold_steps = p_get(params, "break_hold_steps", 5, int)
    damage_record_start_cycle = p_get(params, "damage_record_start_cycle", 2, int)
    max_broken_fraction_stop = p_get(params, "max_broken_fraction_stop", 0.70, float)

    total_steps = int(n_cycles * steps_per_cycle)
    damage_start_step = int(damage_record_start_cycle * steps_per_cycle)
    n_nodes = r0.shape[0]
    n_edges = len(edge_data["i_idx"])
    soft_map, break_map = get_type_damage_params(params)
    ctype = edge_data["contact_type"]
    softening_delta = np.array([soft_map.get(t, 0.04) for t in ctype], dtype=float)
    break_threshold = np.array([break_map.get(t, 0.10) for t in ctype], dtype=float)

    r = r0.astype(float).copy()
    r = r - np.mean(r, axis=0, keepdims=True)
    r_native = r.copy()
    v = np.zeros_like(r)
    active_mask = np.ones(n_edges, dtype=bool)
    break_counter = np.zeros(n_edges, dtype=int)
    broken_step = np.full(n_edges, -1, dtype=int)
    first_break = None
    max_rmsd = 0.0

    for step in range(total_steps):
        phase = 2.0 * math.pi * (step % steps_per_cycle) / steps_per_cycle
        F, energy, rel_static = compute_forces_damage(r, phase, edge_data, epsilon0, active_mask, softening_delta, use_softening)
        a = (F - gamma * v) / mass
        v = v + dt * a
        r = r + dt * v
        if remove_center:
            r, v = remove_com(r, v)
        F2, energy2, rel_static2 = compute_forces_damage(r, phase, edge_data, epsilon0, active_mask, softening_delta, use_softening)
        if use_breaking and step >= damage_start_step:
            over = (rel_static2 > break_threshold) & active_mask
            break_counter[over] += 1
            break_counter[~over] = 0
            newly = (break_counter >= break_hold_steps) & active_mask
            if np.any(newly):
                new_indices = np.where(newly)[0]
                for e in new_indices:
                    active_mask[e] = False
                    broken_step[e] = step
                    if first_break is None:
                        erow = edges_df.iloc[int(e)]
                        first_break = {
                            "edge_index": int(e),
                            "contact_type": str(ctype[e]),
                            "pair": f"{erow['aa_i']}{int(erow['res_i'])}-{erow['aa_j']}{int(erow['res_j'])}",
                        }
        rmsd = float(np.sqrt(np.mean(np.sum((r - r_native) ** 2, axis=1))))
        max_rmsd = max(max_rmsd, rmsd)
        if float(np.mean(~active_mask)) >= max_broken_fraction_stop:
            break

    broken_fraction = float(np.mean(~active_mask)) if n_edges else np.nan
    lcc = largest_connected_component_fraction(n_nodes, edge_data["i_idx"], edge_data["j_idx"], active_mask) if n_edges else np.nan
    type_broken = {}
    for t in CONTACT_ORDER:
        mask = ctype == t
        type_broken[f"{t}_broken_fraction_at_max_epsilon"] = float(np.mean(~active_mask[mask])) if np.any(mask) else np.nan

    return {
        "epsilon0": epsilon0,
        "n_broken_final": int(np.sum(~active_mask)),
        "broken_fraction_final": broken_fraction,
        "lcc_fraction_final": lcc,
        "max_rmsd": max_rmsd,
        "first_break": first_break,
        "type_broken": type_broken,
    }


def simulate_destructive(nodes_df: pd.DataFrame, edges_df: pd.DataFrame, params: Dict[str, str]) -> Dict[str, object]:
    r0, edge_data = prepare_edge_data_damage(nodes_df, edges_df)
    eps_min = p_get(params, "epsilon_min", 0.01, float)
    eps_max = p_get(params, "epsilon_max", 0.30, float)
    eps_points = p_get(params, "epsilon_points", 8, int)
    eps_values = np.linspace(eps_min, eps_max, eps_points)
    results = []
    for eps in eps_values:
        results.append(simulate_destructive_single_epsilon(nodes_df, edges_df, r0, edge_data, params, float(eps)))

    damaging = [r for r in results if r["n_broken_final"] > 0]
    onset = damaging[0]["epsilon0"] if damaging else np.nan
    first = damaging[0]["first_break"] if damaging and damaging[0].get("first_break") else None
    final = results[-1] if results else {}

    out = {
        "epsilon_min": eps_min,
        "epsilon_max": eps_max,
        "epsilon_points": eps_points,
        "first_damaging_epsilon_bin": onset,
        "damage_onset_bin_index": int(np.where(np.isclose(eps_values, onset))[0][0]) if np.isfinite(onset) else -1,
        "first_break_contact_type": first["contact_type"] if first else "",
        "first_break_pair": first["pair"] if first else "",
        "first_break_edge_index": first["edge_index"] if first else -1,
        "final_epsilon0": float(final.get("epsilon0", np.nan)),
        "n_broken_final_at_max_epsilon": int(final.get("n_broken_final", 0)),
        "final_broken_fraction_at_max_epsilon": float(final.get("broken_fraction_final", np.nan)),
        "lcc_fraction_final_at_max_epsilon": float(final.get("lcc_fraction_final", np.nan)),
        "max_rmsd_at_max_epsilon": float(final.get("max_rmsd", np.nan)),
    }
    out.update(final.get("type_broken", {f"{t}_broken_fraction_at_max_epsilon": np.nan for t in CONTACT_ORDER}))
    return out

# =============================================================================
# SONAR label row and worker
# =============================================================================

def build_compact_sonar_row(input_row: dict, params: Dict[str, str]) -> Dict[str, object]:
    pdb_id = str(input_row.get("pdb_id") or input_row.get("protein_id") or "").strip()
    chain_id = str(input_row.get("chain_id") or "").strip()
    file_path = str(input_row.get("file_path") or input_row.get("input_file") or "").strip()
    rank = input_row.get("stage1A_rank", input_row.get("rank", ""))

    row = {
        "pdb_id": pdb_id,
        "protein_id": pdb_id,
        "chain_id": chain_id,
        "file_path": file_path,
        "stage_rank": rank,
        "status": "failed",
        "reason": "",
    }

    if not file_path or not os.path.exists(file_path):
        row["reason"] = f"input file not found: {file_path}"
        return row

    try:
        structure = load_structure_any(file_path, pdb_id)
        nodes_df, payloads = build_nodes_from_structure(structure, pdb_id, chain_id, params)
        n_nodes = len(nodes_df)
        min_res = p_get(params, "min_chain_residues", 40, int)
        max_res = p_get(params, "max_chain_residues", 2000, int)
        if n_nodes < min_res:
            raise ValueError(f"too few parsed residues: {n_nodes} < {min_res}")
        if n_nodes > max_res:
            raise ValueError(f"too many parsed residues: {n_nodes} > {max_res}")

        edges_df = build_edges_v1_1(nodes_df, payloads, pdb_id, chain_id, params)
        if len(edges_df) == 0:
            raise ValueError("no edges constructed")

        net = build_network_summary(nodes_df, edges_df)
        row.update(net)

        run_reversible = p_get(params, "run_reversible", 1, bool)
        run_destructive = p_get(params, "run_destructive", 1, bool)

        if run_reversible:
            rev = simulate_reversible(nodes_df, edges_df, params)
            row.update(rev)
        else:
            row.update({
                "epsilon_reversible": np.nan,
                "mean_reversible_edge_deformation": np.nan,
                "final_mean_reversible_edge_deformation": np.nan,
                "max_rmsd_reversible": np.nan,
                "max_broken_fraction_proxy_reversible": np.nan,
                "mean_elastic_energy_reversible": np.nan,
                "max_elastic_energy_reversible": np.nan,
                "top_hotspot_residue": "",
                "top_hotspot_score": np.nan,
            })
            for t in CONTACT_ORDER:
                row[f"{t}_reversible_response"] = np.nan

        if run_destructive:
            des = simulate_destructive(nodes_df, edges_df, params)
            row.update(des)
        else:
            row.update({
                "epsilon_min": np.nan,
                "epsilon_max": np.nan,
                "epsilon_points": np.nan,
                "first_damaging_epsilon_bin": np.nan,
                "damage_onset_bin_index": -1,
                "first_break_contact_type": "",
                "first_break_pair": "",
                "first_break_edge_index": -1,
                "final_epsilon0": np.nan,
                "n_broken_final_at_max_epsilon": np.nan,
                "final_broken_fraction_at_max_epsilon": np.nan,
                "lcc_fraction_final_at_max_epsilon": np.nan,
                "max_rmsd_at_max_epsilon": np.nan,
            })
            for t in CONTACT_ORDER:
                row[f"{t}_broken_fraction_at_max_epsilon"] = np.nan

        # raw SONAR components; dataset-level normalized fields are added after merging.
        row["sonar_vulnerability_raw"] = (
            row.get("final_broken_fraction_at_max_epsilon", np.nan)
            * (row.get("epsilon_max", np.nan) / row.get("first_damaging_epsilon_bin", np.nan))
        ) if np.isfinite(safe_float(row.get("first_damaging_epsilon_bin"))) and safe_float(row.get("first_damaging_epsilon_bin")) > 0 else np.nan
        row["sonar_integrity_index_raw"] = safe_float(row.get("lcc_fraction_final_at_max_epsilon")) * (1.0 - safe_float(row.get("backbone_broken_fraction_at_max_epsilon"), 0.0))
        row["status"] = "success"
        row["reason"] = ""
        return row

    except Exception as exc:
        row["status"] = "failed"
        row["reason"] = str(exc)
        if p_get(params, "write_traceback_in_failure_reason", 0, bool):
            row["reason"] = traceback.format_exc()
        return row


def worker(args):
    input_row, params = args
    return build_compact_sonar_row(input_row, params)

# =============================================================================
# Batch orchestration
# =============================================================================

def add_dataset_level_indices(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    success = df["status"].astype(str) == "success" if "status" in df.columns else pd.Series(False, index=df.index)
    df["sonar_reversible_index"] = np.nan
    df["sonar_vulnerability_index"] = np.nan
    df.loc[success, "sonar_reversible_index"] = minmax_normalize_series(df.loc[success, "mean_reversible_edge_deformation"])
    df.loc[success, "sonar_vulnerability_index"] = minmax_normalize_series(df.loc[success, "sonar_vulnerability_raw"])
    # Composite placeholder: high reversible response + high vulnerability + low integrity. Keep transparent.
    df["sonar_composite_index"] = np.nan
    comp = (
        pd.to_numeric(df["sonar_reversible_index"], errors="coerce")
        + pd.to_numeric(df["sonar_vulnerability_index"], errors="coerce")
        + (1.0 - pd.to_numeric(df["sonar_integrity_index_raw"], errors="coerce"))
    ) / 3.0
    df.loc[success, "sonar_composite_index"] = comp.loc[success]

    rev_med = df.loc[success, "sonar_reversible_index"].median(skipna=True)
    vul_med = df.loc[success, "sonar_vulnerability_index"].median(skipna=True)
    classes = []
    for _, r in df.iterrows():
        if str(r.get("status")) != "success" or pd.isna(r.get("sonar_reversible_index")) or pd.isna(r.get("sonar_vulnerability_index")):
            classes.append("unclassified")
        else:
            a = "high_reversible" if float(r["sonar_reversible_index"]) >= rev_med else "low_reversible"
            b = "high_vulnerability" if float(r["sonar_vulnerability_index"]) >= vul_med else "low_vulnerability"
            classes.append(f"{a}__{b}")
    df["sonar_response_class"] = classes
    return df


def read_input_rows(params: Dict[str, str]) -> List[dict]:
    input_table = p_get(params, "input_table", None, str)
    if not input_table or not os.path.exists(input_table):
        raise FileNotFoundError(f"input_table not found: {input_table}")
    df = pd.read_csv(input_table)
    required = {"file_path", "chain_id"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"input_table missing required columns: {sorted(missing)}")
    if "pdb_id" not in df.columns and "protein_id" not in df.columns:
        raise ValueError("input_table must contain pdb_id or protein_id")

    start_rank = p_get(params, "start_rank", 1, int)
    end_rank = p_get(params, "end_rank", 0, int)
    max_jobs = p_get(params, "max_jobs", 0, int)
    if "stage1A_rank" in df.columns:
        df = df[pd.to_numeric(df["stage1A_rank"], errors="coerce") >= start_rank]
        if end_rank and end_rank > 0:
            df = df[pd.to_numeric(df["stage1A_rank"], errors="coerce") <= end_rank]
    if max_jobs and max_jobs > 0:
        df = df.head(max_jobs)
    return df.to_dict(orient="records")


def already_done_keys(output_dir: Path) -> set:
    keys = set()
    for p in sorted(output_dir.glob("sonar2_labels_part_*.csv")):
        try:
            df = pd.read_csv(p, usecols=lambda c: c in {"pdb_id", "protein_id", "chain_id", "status"})
        except Exception:
            continue
        id_col = "pdb_id" if "pdb_id" in df.columns else "protein_id"
        for _, r in df.iterrows():
            keys.add((str(r.get(id_col, "")), str(r.get("chain_id", ""))))
    return keys


def stable_field_order(rows: List[dict]) -> List[str]:
    preferred = [
        "pdb_id", "protein_id", "chain_id", "stage_rank", "file_path", "status", "reason",
        "n_nodes", "n_edges", "contact_density_per_residue", "graph_edge_density", "mean_degree", "median_degree", "max_degree", "n_isolated_nodes",
    ]
    for t in CONTACT_ORDER:
        preferred += [f"n_{t}_edges", f"{t}_fraction"]
    preferred += [
        "epsilon_reversible", "mean_reversible_edge_deformation", "final_mean_reversible_edge_deformation",
        "max_rmsd_reversible", "max_broken_fraction_proxy_reversible", "mean_elastic_energy_reversible", "max_elastic_energy_reversible",
        "top_hotspot_residue", "top_hotspot_score",
    ]
    for t in CONTACT_ORDER:
        preferred.append(f"{t}_reversible_response")
    preferred += [
        "epsilon_min", "epsilon_max", "epsilon_points", "first_damaging_epsilon_bin", "damage_onset_bin_index",
        "first_break_contact_type", "first_break_pair", "first_break_edge_index", "final_epsilon0",
        "n_broken_final_at_max_epsilon", "final_broken_fraction_at_max_epsilon", "lcc_fraction_final_at_max_epsilon", "max_rmsd_at_max_epsilon",
    ]
    for t in CONTACT_ORDER:
        preferred.append(f"{t}_broken_fraction_at_max_epsilon")
    preferred += [
        "sonar_vulnerability_raw", "sonar_integrity_index_raw", "sonar_reversible_index", "sonar_vulnerability_index", "sonar_composite_index", "sonar_response_class",
    ]
    all_keys = []
    seen = set()
    for k in preferred:
        if k not in seen:
            all_keys.append(k); seen.add(k)
    for r in rows:
        for k in r.keys():
            if k not in seen:
                all_keys.append(k); seen.add(k)
    return all_keys


def write_csv_rows(path: Path, rows: List[dict], fieldnames: List[str]):
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def merge_parts(output_dir: Path):
    parts = sorted(output_dir.glob("sonar2_labels_part_*.csv"))
    if not parts:
        return None
    dfs = []
    for p in parts:
        try:
            dfs.append(pd.read_csv(p))
        except Exception:
            pass
    if not dfs:
        return None
    df = pd.concat(dfs, axis=0, ignore_index=True)
    df = add_dataset_level_indices(df)
    merged_path = output_dir / "sonar2_labels_merged.csv"
    failed_path = output_dir / "sonar2_failed_chains.csv"
    df.to_csv(merged_path, index=False, encoding="utf-8-sig")
    df[df["status"].astype(str) != "success"].to_csv(failed_path, index=False, encoding="utf-8-sig")
    return df


def run_batch(params: Dict[str, str]) -> None:
    output_dir = Path(p_get(params, "output_dir", "sonar2_stage1B_outputs", str))
    log_dir = Path(p_get(params, "log_dir", str(output_dir / "logs"), str))
    ensure_dir(output_dir); ensure_dir(log_dir)

    n_cpu = max(1, min(p_get(params, "n_cpu", 1, int), cpu_count()))
    chunk_size = p_get(params, "chunk_size", 100, int)
    resume = p_get(params, "resume", 1, bool)
    rows = read_input_rows(params)

    if resume:
        done = already_done_keys(output_dir)
        rows = [r for r in rows if (str(r.get("pdb_id") or r.get("protein_id") or ""), str(r.get("chain_id") or "")) not in done]
    else:
        done = set()

    print("=== SONAR2 single-file SONAR-Calc v1.0 ===")
    print(f"Input jobs after filtering/resume: {len(rows)}")
    print(f"Already done keys: {len(done)}")
    print(f"Output dir: {output_dir}")
    print(f"n_cpu: {n_cpu}")
    print(f"chunk_size: {chunk_size}")
    print(f"run_reversible: {p_get(params, 'run_reversible', 1, bool)}")
    print(f"run_destructive: {p_get(params, 'run_destructive', 1, bool)}")

    t0 = time.time()
    buffer = []
    part_index = 1
    existing_parts = sorted(output_dir.glob("sonar2_labels_part_*.csv"))
    if existing_parts:
        nums = []
        for p in existing_parts:
            try:
                nums.append(int(p.stem.split("_")[-1]))
            except Exception:
                pass
        part_index = max(nums or [0]) + 1

    processed = 0
    fieldnames = None

    def flush(buf, idx):
        if not buf:
            return
        tmp_df = add_dataset_level_indices(pd.DataFrame(buf))
        out_rows = tmp_df.to_dict(orient="records")
        fields = stable_field_order(out_rows)
        out_path = output_dir / f"sonar2_labels_part_{idx:04d}.csv"
        write_csv_rows(out_path, out_rows, fields)
        print(f"[WRITE] {out_path} ({len(out_rows)} rows)")

    try:
        if n_cpu == 1:
            iterator = (worker((r, params)) for r in rows)
        else:
            pool = Pool(processes=n_cpu)
            iterator = pool.imap_unordered(worker, [(r, params) for r in rows], chunksize=1)

        for result in iterator:
            buffer.append(result)
            processed += 1
            if processed % 10 == 0:
                elapsed = time.time() - t0
                rate = processed / elapsed if elapsed > 0 else 0
                ok = sum(1 for x in buffer if x.get("status") == "success")
                print(f"Processed {processed}/{len(rows)} | {rate:.3f} jobs/s | current-buffer successes={ok}/{len(buffer)}")
            if len(buffer) >= chunk_size:
                flush(buffer, part_index)
                buffer = []
                part_index += 1
        if n_cpu != 1:
            pool.close(); pool.join()
    except KeyboardInterrupt:
        if n_cpu != 1:
            pool.terminate(); pool.join()
        print("[INTERRUPTED] Writing current buffer before exit.")
        flush(buffer, part_index)
        raise

    if buffer:
        flush(buffer, part_index)

    merged = merge_parts(output_dir)
    elapsed = time.time() - t0
    log_path = log_dir / "sonar2_singlefile_run_log.txt"
    with open(log_path, "w", encoding="utf-8") as f:
        f.write("SONAR2 single-file SONAR-Calc v1.0 run log\n")
        f.write(f"elapsed_seconds={elapsed:.2f}\n")
        f.write(f"new_jobs_processed={processed}\n")
        f.write(f"n_cpu={n_cpu}\n")
        if merged is not None:
            f.write(f"merged_rows={len(merged)}\n")
            f.write(str(merged["status"].value_counts()) + "\n")
    print("\n=== DONE ===")
    print(f"Elapsed seconds: {elapsed:.2f}")
    if merged is not None:
        print(f"Merged labels: {output_dir / 'sonar2_labels_merged.csv'}")
        print(merged["status"].value_counts())
    print(f"Run log: {log_path}")


def main():
    ap = argparse.ArgumentParser(description="Single-file SONAR2 SONAR-Calc batch engine")
    ap.add_argument("--params", required=True, help="Path to key,value parameter CSV")
    args = ap.parse_args()
    params = read_params(args.params)
    run_batch(params)


if __name__ == "__main__":
    main()
