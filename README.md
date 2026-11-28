# Certifying the false-positive rate of a frozen classifier at a new site from a few annotated fields

Code, sealed preregistration documents and result files for the manuscript *"Certifying the false-positive rate of a frozen classifier at a new site from a few annotated fields: limits and rules under cluster sampling"* (submitted to *Biometrical Journal*).

## Layout
| Folder | Content |
|---|---|
| `results/results_master.json.gz` | Master result file (gzip, deterministic). Every number in manuscript versions up to 6.31 and their Supporting Information is a key of this file (`values`), with its source file in `sources`. The public copy is redacted with the placeholders listed below: `MANIFEST.tsv` gives the sha256 of the original uncompressed file (column `original_sha256`), and `REDACTIONS.tsv` gives the original and the public sha256 of the uncompressed file. |
| `results/results_master_r3.json.gz` | Superset of `results/results_master.json.gz`: every key of that file is present with the same value (0 keys changed), plus the keys of the two preregistered robustness addenda (number of fields, small-sample corrections and tuning constants, and the site-calibrated generators with their calibration) and the keys derived from them. Every number in manuscript version 6.32 and its Supporting Information is a key of this file. It is redacted with the same placeholders, and `MANIFEST.tsv` and `REDACTIONS.tsv` give its original sha256. The addenda are `preregistration/PREREG_P5A_S2_ADDENDUM_20261128.md` and `preregistration/PREREG_P5A_S3_ADDENDUM_20261128.md`, and the key builders are `code/build/add_s2_keys_r2.py`, `add_s3_keys_r3.py`, `add_s3_derived_r3a.py` and `add_s3calib_keys_r3b.py`. The site-calibrated simulation of the second addendum needs the development-site calibration file, whose sha256 is in this master (`s3sim/calib|npz`), so it cannot be re-run from the repository alone. |
| `results/results_src/` | The per-analysis result files that the master is built from (aggregates per site, cell and rule). |
| `code/results_src/` | Analysis and simulation scripts, in the folder of the result they produce, plus their shared libraries (`e_batch`, `m3`, `c5b`, `e10`, `e12`, `simw3`, `sim50`, `sizefpr`). |
| `code/build/` | Number filling (`fill_numbers_v3.py`, which the build uses, and the earlier `fill_numbers_v2.py`), the four build checks, the key builders (`add_*_keys.py`), the derivation scripts (`d_*.py`) and the table and figure renderers. The manuscript assembly and its text-edit modules carry manuscript prose and are released after acceptance (`EXCLUDED.tsv`). |
| `code/figures/` | Scripts for the nine main-text figures, the style module `p5a_style_ieee.py` and the input loader `public_inputs.py` (see Figures). |
| `preregistration/` | Sealed specification documents of every evaluation, each fixed before the numbers it governs, with their `.sha256` files where they exist. Most are written in Chinese. The status headers of the two 2026-11-28 addenda predate review. Both were sealed and run on 2026-11-28, and their run records are appended. |

## Reproducing tables and figures
1. Every number in manuscript versions up to 6.31 and their Supporting Information is filled from a key of `results_master.json` (version 6.32: `results_master_r3.json.gz`) by `code/build/fill_numbers_v3.py`, which substitutes `[[key :: format]]` placeholders. The manuscript sources themselves are released on acceptance.
2. The simulation study (Section 4) is reproduced by `code/results_src/e10/e10_synth.py`, `code/results_src/e12/e12_stress.py`, `code/results_src/simw3/simw3_baselines.py` and `code/results_src/simw3/sim50_extend.py` (50 synthetic sites per cell). These scripts use synthetic data only and need no external data set.
3. The development-site and confirmation analyses read per-cell classifier scores of public image data sets. Those score files are **not** included (see Data and licences). The scripts document the expected inputs and pin them by sha256.
4. Scripts pin the files they read by sha256. In this repository, local paths and host names were replaced by placeholders (`<P5_ROOT>`, `<VM_HOME>`, `<DATA>`, `<LOCAL_PATH>`, `<LOCAL_TMP>`, `<VM_HOST>`, `<USER>`). The 44 affected files, the master result file included, are listed in `REDACTIONS.tsv` with their original and public sha256. A pin in one script can therefore refer to the original, unredacted copy of another file.

## Figures
Each script reads the master through `public_inputs.load_master`, which accepts `results/results_master.json.gz` (the public, redacted copy) and checks its sha256. Run from `code/figures/` with Python 3, NumPy, SciPy and Matplotlib, for example `python3 make_fig_simheat_v66.py --master ../../results/results_master.json.gz --eps 0.1 --out <new directory>`. Each script writes the PDF, a PNG and a `_numbers.json` file that lists every number shown with its master key.

| Figure | Script | Inputs | From this repository alone |
|---|---|---|---|
| 1 | `make_fig_simblocks_v67r3.py` | master | yes |
| 2 | `make_fig_simheat_v66.py` (`--eps 0.1`) | master | yes |
| 3 | `make_fig_simbox_v66.py` (`--eps 0.1`) | master | yes |
| 4 | `make_fig_drift_v42.py` (`--s6a ../../results/results_src/s6/s6a_survey.json`) | master and the per-checkpoint AUROC and FPR rows of panels a and b | yes |
| 5 | `make_fig_mechanism_v39.py` | master, `results_src/figmech_v39/figmech_v39_numbers.json` and `figmech_v39_arrays.npz` (per-cell scores of the drawn fields) | no: the per-cell scores are not redistributed (`EXCLUDED.tsv`) |
| 6 | `make_fig_funnel_v65.py` | master | yes |
| 7 | `make_fig_devforest_v65.py` | master | yes |
| 8 | `make_fig_dev_v39.py` | master | yes |
| 9 | `make_fig_confforest_v65.py` (`--de ../../results/results_src/deconfirm/d_e_confirm_intervals.json`) | master and the D-e interval file | yes; per-draw values were not stored for PlasmoCount or for the sensitivity difference on M5, so those points have no interval |

A script whose input is not in this repository stops with exit code 3 and a message naming the file. The rendered PDFs are not byte-identical to the submitted ones, because Matplotlib writes the creation date into each PDF; the `_numbers.json` files can be compared instead.

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
