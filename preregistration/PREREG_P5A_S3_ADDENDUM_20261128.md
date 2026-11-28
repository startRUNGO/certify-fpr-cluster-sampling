# P5a prereg addendum S3: site-calibrated simulation and Split-Select constants on the development sites (DRAFT v1, 2026-11-28, for line-by-line review by 实验管理)

**Source.** User approval relayed by 实验管理 on 2026-11-28.

**STATUS: DRAFT.** Nothing has run. Both scripts and the wrapper are written and checked with ast.parse and bash -n only. Every number below "Fixed now" is fixed before any run.

## 0. Weakness, change and check
- **Weakness 1.** The simulation generators (E10, E12, size-FPR) are stylised. A referee can say that coverage on them does not transfer to a real site.
  - **Change:** a generator calibrated to one development site, BBBC-ALL, plus a second cell calibrated to Lacuna-B, where the cluster baselines fail (Table 5) and DE-Grid's coverage is lowest.
  - **Check:** prediction 1 below.
- **Weakness 2.** The variance-ratio condition κ_W² ≥ (z_{1−2α}/z_{1−α})² is checked only on the population (vardom/kappaw: minimum 0.794 over 90 development cells). How often a single draw falls below the threshold is unknown.
  - **Change:** the per-draw κ̂² defined in §1.
  - **Check:** descriptive (§1).
- **Weakness 3.** The tuning-constant sensitivity of S2-C has no Split-Select row on real data, because dev C ran DE-Grid and B1u only.
  - **Change:** Split-Select on the development sites through its own path, with the stored Platt maps.
  - **Check:** prediction 2.
- **Design correction to the request (2026-11-28).**
  - The Split-Select development m-grid (S3(2)(a)) is not run.
  - It already exists for SPLIT6, S3 and the two cell-level baselines, at m ∈ {6, 10, 15, 22, 30, 40, 60} on the three development sites. The keys are ebatch_e2p/cells|<site>|<m>|<eps>|<rule>|…, and main Figure 8 shows them.
  - The round-3 statement that Split-Select had no m-grid was wrong. Rerunning it would only re-confirm a published result.

## 1. Fixed now

### Part 1: site-calibrated simulation (`p5a_s3_sim.py`)
- **Step calib** (a VM job, development-site scores only, read through m3_replay_lib.Inputs with manifest checks). Per site (BBBC-ALL, Lacuna-B) and A2 seed (41–45), it computes:
  - the fields with at least one negative and the empirical (n_neg, n_pos) pair of each field
  - the field negative means m_f and μ0 = mean(m_f)
  - the pooled within-field negative residuals e = s − m_f
  - the pooled positive offsets d = s_pos − m_f
  - τ_b² = max(0, var(m_f) − mean(s_w²/n_f)), a one-way moment estimate with the pooled within-field variance s_w²
  - **Output:** s3calib.npz and s3calib.json. The per-cell pools stay on the VM. Both shas go into this addendum (r1) before the sim step runs.
- **Generator** (per realisation r = 41–90, calibration seed cs = 41 + (r − 41) mod 5, stream [20261129, 3, site code, r, 0]):
  - F equals the site's field count.
  - Each field's (n_neg, n_pos) pair is resampled with replacement from the site's pairs.
  - μ_f ~ N(0, τ_b²).
  - Negatives are μ0 + μ_f + e and positives μ0 + μ_f + d, with e and d resampled from their pools.
  - So the site's joint field-size and prevalence distribution, its between-field shift variance and its within-field score shape are kept.
- **Rules.**
  - DE-Grid (S3) and Split-Select (SPLIT6) through e12_stress.run_site_seed_e12 with the uniform sampler (sim50 pass 1).
  - m = 22, 200 draws, ε ∈ {0.10, 0.05}, 50 realisations per cell.
  - Split-Select's control variate uses the stored Platt map of the calibration seed (sp_dev.json 364ebc08, key platt).
- **Per-draw κ̂², descriptive** (per draw, ε and candidate k):
  - τ* is the largest threshold with population FPR > ε: the ⌊εN⌋ + 1-th largest population negative, as in kappaw_margin.
  - V_true = (1 − m/F)/m · Σ_f (x_f − p* n_f)² / ((F − 1) n̄²) at τ*.
  - κ̂² = p̂(1 − p̂)(1 + (n̄_w − 1)ρ_U)/(m n̄_d) / V_true, where:
    - p̂ is the draw's pooled FPR
    - ρ_U is the ICC upper limit that b1u uses at level δ/3 (the script's copy is checked against m3_replay_lib.b1u)
    - n̄_d and n̄_w = Σn²/Σn are the draw's own
  - **Reported:**
    - the share of draws with κ̂² < 0.66997 (vardom/kappa2_threshold), per realisation (mean and MCSE) and pooled
    - the share with p̂ = 0
    - the coverage of S3 and SPLIT6 among draws with κ̂²(c1s) below and above the threshold
  - **Judgement:** descriptive, with no bound. The population result (minimum κ_W² 0.794) already states the condition. The per-draw share shows how often a single draw's design-effect variance falls short of the true variance, and whether coverage is lower in those draws. A conditional coverage below 1 − δ is reported as a measured limit.
- **Gates (exit 3, nothing written):**
  - module and file pins, including the calibration shas given to the sim step
  - feasibility: F ≥ 22 for every site and seed
  - the ICC-copy check (50 fixed vectors)
  - anchor c rerun on the calibrated generator: realisation 41, 5 draws, run_site_seed_e12 (uniform) equals SP.run_site_seed for S3, B1u_c0, B1u_c1s and SPLIT6
  - There is no stored-row reproduction gate, because the generator is new.
- **Prediction 1, fixed now:** at ε 0.10, coverage ≥ 1 − δ − 2·MCSE in both cells for both rules. MCSE is the SD of the 50 per-realisation coverages divided by √50. ε 0.05 is reported only. If a cell fails, it is reported as a measured limit and no rule is changed.
- **Calibration diagnostics, reported, no gate.** Per site and seed: τ_b, s_w, F, prevalence and the pool sizes (s3calib.json).

### Part 2(b): Split-Select constant variants on the development sites (`p5a_s3_devsplit.py`)
- **Path.** m3_split_dev_eval.run_site_seed (cbf1e7f4), the routine behind sp_dev.json.
  - The three rebuilt sites, A2 seeds 41–45, m = 22, 200 draws, both ε.
  - The Platt map is the one sp_dev.json stores.
  - L.DEV is extended with mpidb_old after the three sites, exactly as in m3_split_dev_eval.main, so the stream indices are those of sp_dev.
- **Variants (one at a time):** default, n_s ∈ {10, 40}, r_low ∈ {0.60, 0.80}, G ∈ {100, 400}. These are applied through p5a_s2_common.set_constants; at the default the frozen cstat object is kept.
- **Gates (exit 3):**
  - m3_replay_lib.self_check of the default rows
  - every default S3 and SPLIT6 cell (3 sites × 2 ε) equals sp_dev.json in n, coverage, fpr_mean, cannot_certify, sens, ds_eps and ds_own, within 1e-12
  - check_cstat_copy
- **Prediction 2:** the same form as S2-C, descriptive on five seeds. Split-Select coverage stays ≥ 0.90 − 2·(seed MCSE) at ε 0.10 for every variant and site. Sensitivity and cannot-certify may move. Five seeds give a coarse MCSE, so an exception is reported, not judged.

## 2. Run (VM only, after review and 实验管理's OK)
- **Wrapper.** `run_p5a_s3.sh <calib|sim|dev> <tag> [workers] [calib_dir calib_npz_sha calib_json_sha]`. It is a copy of S2 wrapper r1 with the same gates:
  - package sha256sum -c
  - house RUN.log format, with a tail-1 ' START ' busy gate
  - nvidia-smi failure counts as busy
  - MemAvailable ≥ 32 GB
  - nice 19, taskset 56–63, 8 workers, 1 BLAS thread
  - watchdog RSS 24 GB and wall 4 h
  - job directory <DATA>/p5s3_<kind>_<tag>
- **Package p5s3_pkg_<tag>:** results_src from S2 package e681ec81 (all pinned), s2/p5a_s2_common.py (r1 0b02fb98), s3/ with the two scripts, and inputs/sp_dev.json.
- **Order:** calib → fix the calibration shas in addendum r1 → sim → dev. One job at a time.
- **Time estimate:**
  - calib: minutes
  - sim: 100 tasks at ≈ 15–25 s on 8 workers, ≈ 5 min, plus the anchor
  - dev: 105 tasks with the K1/K1vB bootstraps of run_site_seed. A timed first task will show the per-task cost; the estimate is ≤ 1 h.
- **Write volume:** under 120 MB. The calib pools are about 60 MB compressed and stay on the VM.
- **Results:** an add-only P5a master r3 with s3sim/ and s3dev/ keys. They go to the main text only as one sentence each, with details in the SI.


## r1 (2026-11-28T09:40:40Z; 实验管理 line-by-line review; replaces the Part 1 generator, κ̂² and summary text above where they differ; sim script r1 10b8ba89)

### Generators: two, neither primary
- **A (as drafted).** Field pairs are resampled. μ_f ~ N(0, τ_b²) is independent of the pair, and residuals and offsets come from the pooled pools. Stream [20261129, 3, site code, r, 0].
  - A keeps the site's size and prevalence distribution, its between-field variance and its within-field shape.
  - A cannot make a field's exceedance rate depend on its size. That coupling is the mechanism behind κ < 1 (SI S8).
- **B (source-field coupled, new).** Stream [20261129, 3, site code, r, 1]. Each synthetic field copies one source field f*, drawn uniformly with replacement from the site's fields (per site and calibration seed):
  - its (n_neg, n_pos)
  - its negative mean m_f*
  - its **own** negative residuals and its **own** positive offsets, each resampled with replacement
  - no extra shift
  - So the size–mean–shape coupling of each real field is kept. Every field with positives has its own offsets, so no pooled positive draw is needed.
- The calibration npz additionally stores m_f and the per-field start indices into the residual and offset pools. The volume is about the same.
- **Cells:** 2 generators × 2 sites × 50 realisations. Prediction 1 applies to both generators and is fixed now. A failure in either is a measured limit, and no rule changes. Anchor c runs on both generators.

### κ̂², primary and secondary
- **Primary:** κ̂²_W = p*(1 − p*)[1 + (n̄_w,d − 1)ρ_U,d]/(m n̄_d) / V_true, where p* is the population FPR at τ*. This matches SI S8, which computes κ_W² with the site's boundary FPR.
- **Secondary:** the same with the draw's p̂(1 − p̂). It is 0 whenever p̂ = 0, so it is reported as a column with the p̂ = 0 share.
- **Conditional coverage** (S3 and SPLIT6, below versus above 0.66997) splits on the **minimum over the three candidates** of κ̂²_W, as in the min_k of the DE-Grid familywise bound. The split on c1s's κ̂²_W is a second report.
- **V_true = 0:** κ̂² is not defined. Such draws are counted as n_nan in every κ̂² summary and are excluded from the below and above shares.

### Calibration diagnostics
- s3calib.json additionally records n_fields_dropped and n_pos_dropped per site and seed. These are fields with no negative, and their positives.

### Order
- calib (VM, after the busy-gate and MemAvailable check), then s3calib.json goes to 实验管理.
- Then both calibration shas go into this addendum (r2), 实验管理 OKs the sim step, and sim runs, then dev.
- One job at a time.


## r2 (2026-11-28T09:44:26Z; calibration sealed before any sim run)
- **Calib run:** TAG 20261128T094156Z, package p5s3_pkg_20261128T094156Z (SHA256SUMS b962712954cda56f…), rc 0, 09:42:47Z to 09:43:07Z.
- **Calibration files (sealed):**
  - ad11c5712a1915b1aa11c4ae45bdbac4de7721fb8a0795bd89ef1e0623fcca9d  s3calib.json
  - 3c21b040b34e83a56328c12c7f3003279f1b8b1a425c29da05e93b42bf4325ab  s3calib.npz
- **The sim step runs with:** `run_p5a_s3.sh sim 20261128T094156Z 8 <DATA>/p5s3_calib_20261128T094156Z/s3calib_20261128T094156Z 3c21b040b34e83a56328c12c7f3003279f1b8b1a425c29da05e93b42bf4325ab ad11c5712a1915b1aa11c4ae45bdbac4de7721fb8a0795bd89ef1e0623fcca9d` (wrapper order: calib_dir, npz sha, json sha).
- **Diagnostics (descriptive):**
  - BBBC-ALL: F = 1,328, prevalence 0.0287, τ_b 0.66–1.18, s_w 0.71–1.02 across seeds.
  - Lacuna-B: F = 483, prevalence 0.0287, τ_b 0.45–0.56, s_w 1.05–1.19.
  - No field had zero negatives, so 0 fields and 0 positives were dropped.
  - The two prevalences agree to 3 decimals: 2,452/85,486 = 0.02868 and 1,517/52,859 = 0.02870.
