# Certifying the false-positive rate of a frozen classifier at a new site from a few annotated fields

Code, sealed preregistration documents and result files for the manuscript *"Certifying the false-positive rate of a frozen classifier at a new site from a few annotated fields: limits and rules under cluster sampling"* (submitted to *Biometrical Journal*).

## Layout
| Folder | Content |
|---|---|
| `results/results_master.json.gz` | Master result file (gzip, deterministic). Every number in the manuscript and the Supporting Information is a key of this file (`values`), with its source file in `sources`. The public copy is redacted with the placeholders listed below: `MANIFEST.tsv` gives the sha256 of the original uncompressed file (column `original_sha256`), and `REDACTIONS.tsv` gives the original and the public sha256 of the uncompressed file. |
| `results/results_src/` | The per-analysis result files that the master is built from (aggregates per site, cell and rule). |
| `code/results_src/` | Analysis and simulation scripts, in the folder of the result they produce, plus their shared libraries (`e_batch`, `m3`, `c5b`, `e10`, `e12`, `simw3`, `sim50`, `sizefpr`). |
| `code/build/` | Number filling (`fill_numbers_v3.py`, which the build uses, and the earlier `fill_numbers_v2.py`), the four build checks, the key builders (`add_*_keys.py`), the derivation scripts (`d_*.py`) and the table and figure renderers. The manuscript assembly and its text-edit modules carry manuscript prose and are released after acceptance (`EXCLUDED.tsv`). |
| `code/figures/` | Scripts for Figures 1, 2, 5 and the funnel figure. |
| `preregistration/` | Sealed specification documents of every evaluation, each fixed before the numbers it governs, with their `.sha256` files where they exist. Most are written in Chinese. |

## Reproducing tables and figures
1. Every number in the manuscript and the Supporting Information is filled from a key of `results_master.json` by `code/build/fill_numbers_v3.py`, which substitutes `[[key :: format]]` placeholders. The manuscript sources themselves are released on acceptance.
2. The simulation study (Section 4) is reproduced by `code/results_src/e10/e10_synth.py`, `code/results_src/e12/e12_stress.py`, `code/results_src/simw3/simw3_baselines.py` and `code/results_src/simw3/sim50_extend.py` (50 synthetic sites per cell). These scripts use synthetic data only and need no external data set.
3. The development-site and confirmation analyses read per-cell classifier scores of public image data sets. Those score files are **not** included (see Data and licences). The scripts document the expected inputs and pin them by sha256.
4. Scripts pin the files they read by sha256. In this repository, local paths and host names were replaced by placeholders (`<P5_ROOT>`, `<VM_HOME>`, `<DATA>`, `<LOCAL_PATH>`, `<LOCAL_TMP>`, `<VM_HOST>`, `<USER>`). The 40 affected files, the master result file included, are listed in `REDACTIONS.tsv` with their original and public sha256. A pin in one script can therefore refer to the original, unredacted copy of another file.

## Data and licences
- No image and no per-cell or per-image derived value is redistributed. The image data sets are available from their original providers (BBBC041, NIH malaria cell images, MP-IDB, the Lacuna Fund Uganda thin-film archive, PlasmoCount, IML Malaria, M5), cited in the manuscript.
- IML Malaria has no licence file and was used for internal evaluation only. Only aggregate results are reported.
- Files left out are listed in `EXCLUDED.tsv` with the reason.

## Software
Python 3.10.19 with NumPy 2.2.5, SciPy 1.15.3 and PyTorch 2.5.1 (analysis server). Python 3.9.6 and a Python environment with NumPy 2.0.2, SciPy 1.16.3 and Matplotlib 3.10.8 (local checks and figures). pdfTeX from TeX Live 2026 (manuscript build).

## Licence
Code: MIT licence (`LICENSE`). Result files and documents: Creative Commons Attribution 4.0 International (`LICENSE-RESULTS.md`).

## Integrity
`MANIFEST.tsv` lists every file with its size and sha256. `REDACTIONS.tsv` lists the redacted files. `EXCLUDED.tsv` lists the files left out.
