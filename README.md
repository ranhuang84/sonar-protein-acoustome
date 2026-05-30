# SONAR-Protein-Acoustome

**SONAR** (*Structural Oscillation and Network Acoustic Response*) is a structure-based acoustic-elastic network calculator for standardized protein acoustome profiling.

Given a protein structure and a selected chain, SONAR-Calc builds a chemically annotated residue-contact network, simulates reversible and destructive Heterogeneous Acoustic-Elastic Network Dynamics (H-AEND), and outputs standardized SONAR descriptors for protein acoustome analysis.

This repository provides the public single-file implementation corresponding to:

**SONAR: a structure-based acoustic-elastic network calculator for standardized protein acoustome profiling**

---

## Main features

* Builds residue-level acoustic-elastic networks from PDB/mmCIF structures.
* Uses Cβ residue-node representation, with Gly represented by Cα.
* Annotates residue contacts into chemically interpretable layers:

  * backbone
  * disulfide-labeled
  * salt-bridge candidate
  * aromatic
  * hydrophobic
  * polar/H-bond candidate
  * generic weak
* Runs low-amplitude reversible H-AEND simulations.
* Runs high-amplitude destructive H-AEND sweeps.
* Generates standardized report outputs:

  * **SONAR Core Report v1.0**: 44 fields
  * **SONAR Extended Label Matrix v1.0**: 87 fields
* Supports batch processing and CPU control through a single parameter file.

---

## Repository structure

```text
sonar-protein-acoustome/
  README.md
  LICENSE
  sonar_calc_singlefile_v1_0.py
  params_sonar_calc_manuscript_v1_0.csv
  params_sonar_calc_quick_test.csv
  examples/
    hewl_193L/
```

During execution, `sonar_calc_singlefile_v1_0.py` writes embedded frozen SONAR/H-AEND modules into:

```text
_sonar_calc_v1_0_embedded/scripts/
```

This directory is generated automatically and does not need to be edited by users.

---

## Installation

SONAR-Calc requires Python 3 and common scientific Python packages.

Install dependencies with:

```bash
pip install numpy pandas scipy biopython matplotlib
```

If you use a clean Python environment, the following minimal setup is recommended:

```bash
python -m venv sonar_env
sonar_env\Scripts\activate
pip install numpy pandas scipy biopython matplotlib
```

On Linux/macOS:

```bash
python -m venv sonar_env
source sonar_env/bin/activate
pip install numpy pandas scipy biopython matplotlib
```

---

## Quick start

A quick test can be run with:

```bash
python sonar_calc_singlefile_v1_0.py --params params_sonar_calc_quick_test.csv
```

The manuscript-parameter version can be run with:

```bash
python sonar_calc_singlefile_v1_0.py --params params_sonar_calc_manuscript_v1_0.csv
```

The number of CPU workers is controlled in the parameter file:

```csv
n_cpu,8
```

For small tests, use:

```csv
n_cpu,1
```

For larger batch runs, use an appropriate value based on available CPU cores and memory.

---

## Input panel format

The input protein panel should be a CSV file with the following columns:

```csv
protein_id,protein_name,input_file,chain_id,category,notes
```

Example:

```csv
193L,Hen egg-white lysozyme,examples/hewl_193L/input/193L.cif,A,disulfide-rich enzyme,HEWL parity validation
```

The `input_file` field can point to either a PDB or mmCIF structure file.

---

## Output structure

By default, SONAR-Calc writes outputs to:

```text
output/
  networks/
  reversible/
  destructive/
  summary/
  SONAR_report_v1_0/
```

The final standardized report files are written under:

```text
output/SONAR_report_v1_0/
```

Key files include:

```text
SONAR_core_report_v1_0.csv
SONAR_extended_label_matrix_v1_0.csv
SONAR_core_report_schema_v1_0.csv
SONAR_extended_label_matrix_schema_v1_0.csv
```

---

## HEWL parity validation

The public single-file implementation was validated against the manuscript HEWL/193L exemplar.

Expected HEWL parity values include:

```text
n_nodes = 129
n_edges = 498
mean_degree = 7.72093

n_disulfide_edges = 4
n_salt_bridge_edges = 2
n_aromatic_edges = 6

mean_reversible_edge_deformation ≈ 2.78 × 10^-4
top_hotspot_residue = ARG5

first_damaging_epsilon_bin = 0.184
first_break_contact_type = generic_weak
first_break_pair = ARG61-SER72
final_broken_fraction_at_max_epsilon ≈ 0.1888
```

The output schema is fixed as:

```text
SONAR Core Report v1.0 = 44 fields
SONAR Extended Label Matrix v1.0 = 87 fields
```

---

## Interpretation notes

SONAR-Calc is a coarse-grained residue-level acoustic-elastic network model. Edge failure in destructive H-AEND denotes failure of a model-defined residue-level elastic contact. It should not be interpreted as atomistic covalent bond rupture or experimentally calibrated chemical damage.

For proteins with no irreversible contact removal within the scanned destructive range, destructive onset is treated as right-censored at the maximum tested amplitude for report-index computation. These cases are labeled as no damage observed within range.

---

## Citation

If you use SONAR-Calc, please cite the corresponding manuscript:

```text
Huang R. et al. SONAR: a structure-based acoustic-elastic network calculator for standardized protein acoustome profiling.
```

A formal citation and DOI will be added after publication or software release archiving.

---

## License

This project is licensed under the Apache License 2.0.
