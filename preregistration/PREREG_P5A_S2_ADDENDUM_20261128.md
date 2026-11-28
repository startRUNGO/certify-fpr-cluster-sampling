# P5a prereg addendum S2: pre-submission robustness (DRAFT v1, 2026-11-28, for line-by-line review by 实验管理)

**Source.** User approval relayed by 实验管理 on 2026-11-28. All three benchmark reviewers flagged these points as likely referee requests.

**STATUS: DRAFT.** Nothing has run. The scripts and the wrapper are written, compiled and statically checked only, as the no-local-runs rule requires. Every number below the line "Fixed now" is fixed before any run.

## 0. Weakness, change and check

| Part | Weakness found | What S2 adds | Check fixed now |
|---|---|---|---|
| A | The simulation never varies m. All E10, E12, sim50, simw3 and size-FPR settings use m = 22. Development sites have an m grid for DE-Grid only (rev1 m ∈ {3, 5, 10, 22, 30, 50}; m3_stress m ∈ {10, 15, 22}). Split-Select has no m grid anywhere. | Coverage of DE-Grid and Split-Select as a function of m, on the main-study generator | Gate A plus prediction A |
| B | The cluster-aware baselines (CB, GEE, BB) are large-sample. A referee will ask for a small-sample cluster-robust bound. | CR2 (Bell–McCaffrey) Wald with t(m−1), and a wild cluster bootstrap-t | Gate B. Reported only: no prediction on the baselines |
| C | n_s = 20, r_low = 0.70 and G = 200 are fixed in the freeze file without a written rationale (PROTOCOL_FREEZE f2565f30; first stated in PREREG_ITER10 l.7, ITER9 l.9, ITER7 l.9). | One-at-a-time sensitivity of these three constants | Gate C plus prediction C |

**m_s is not part of C.** m_s = 6 has a written design rationale (problem card draft7.1 d2c5a182 §2.4). m_s = 4 and 8 were preregistered as its sensitivity and are already reported (protocol/split_ms_alt_lo/hi).

## 1. Fixed now

### Part A: coverage as a function of m (synthetic, `p5a_s2_sim.py --part A`)
- **Code path.**
  - Sites, samplers and Platt fits come from `sim50_extend.site_platt` (ea3ea087), with the size-FPR registration `sizefpr_check.register` (7d4727be).
  - DE-Grid (S3) and Split-Select (SPLIT6) come from `e12_stress.run_site_seed_e12` (aef35708), pass 1 of sim50. m enters only through `m3_replay_lib.M`.
  - The draw stream `default_rng([20261022, key8, ai, seed, m, sidx])` contains m. So m ≠ 22 gets fresh draws, and m = 22 reuses the main-study draws.
- **Cells.** Fixed from the main-study results, before any S2 number exists:
  - e10:0: the base.
  - e12:19, C2 two-stage, S = 12, ω = 0.8, s = 2. This is the lowest Split-Select cell of the lowest block (sim50 all50 coverage 0.880 at ε 0.10).
  - e12:10, B field sizes, σ = 2.0. The lowest cell of the second-lowest block (0.947).
  - sc:101, size-FPR, σ = 1.0, ρ_c = 0.8. The lowest DE-Grid size-FPR cell (0.989). This is "size-FPR P2".
- **m ∈ {10, 15, 22, 33, 45, 60}.**
  - Split-Select keeps m_s = 6, so m_c = m − 6. At m = 10, m_c = 4.
  - 33 is the per-candidate distribution-free minimum (prop_minm/n_min|delta_over_K).
- **Realisations.** 50 (seeds 41–90), 200 draws each, as in the main study. ε = 0.10 is judged. ε = 0.05 is reported only.
- **Gate A (exit 3, nothing written).** At m = 22, every S3 and SPLIT6 row equals the stored main-study row exactly.
  - Fields compared: cert, cand, fpr, sens, dse, dso.
  - Count: 4 cells × 50 × 200 × 2 ε × 2 rules = 80,000 rows.
  - Sources:
    - e12_rows.json (ece7c768), for e12 seeds 41–45
    - sim50 batches b1 88baad0d, b2 6a04d4b9, b3 c82330c2
    - size-FPR 101 (f89a5730)
  - Stored row counts were checked read-only: sc101 40,000; e12:10 seeds 46–90 36,000.
- **Prediction A (written before the run).**
  - (i) For m ≥ 22, both rules have coverage ≥ 1 − δ − 2·MCSE in every cell at ε 0.10. MCSE is the SD of the 50 per-realisation coverages divided by √50, the sim50 rule.
  - (ii) For m < 22, coverage may fall short, and the cannot-certify rate rises as m falls.
  - Judgement: (i) is supported only if it holds in all 4 cells × 4 values of m ≥ 22 × 2 rules. Any exception is reported as a measured limit, with no change to the rules.
- **Reported per cell, m, rule and ε.** Coverage, MCSE, binomial SE, cannot-certify, ΔSens_ε, ΔSens_own, mean FPR and, for Split-Select, the bound chosen.

### Part B: small-sample cluster-robust baselines (`p5a_s2_sim.py --part B`, `p5a_s2_dev.py --part B`)
- **New bounds** (p5a_s2_common). Single candidate c1s, level δ = 0.10, the same descending grid and fixed-sequence test as the W3 and simw3 baselines.
  - **CR2_c1s.** V = Σ_g (x_g − p̂ n_g)² / (1 − n_g/N) / N², and UB = p̂ + t_{m−1, 1−δ} √V, clipped to [0, 1].
    - This is the CR2 adjustment for the cell-level intercept model, with fields as clusters: H_gg = J/N.
    - X = 0 gives UB = 0, as in GEE. That behaviour is kept and reported.
  - **WCB_c1s.** Wild cluster bootstrap-t with unrestricted residuals, Webb six-point weights and B = 999.
    - The CR2 SE is used in both the numerator and the denominator.
    - UB = p̂ − q_δ(t*)·se.
    - Weights are drawn once per (realisation, draw, ε) and shared along the grid. Streams: `default_rng([20261128, 2, design code, cell, seed, draw, ε code])` in the simulation, and `[20261128, 1, site index, seed, draw, ε code]` on the development sites. Neither touches the draw stream.
  - Mancl–DeRouen is not added. CR2 is the variant asked for first, and one small-sample sandwich is enough for the comparison.
- **Settings.**
  - **Simulation:** all 34 main-study cells (E10 0–9, E12 0–23). Seeds 41–90, m = 22, 200 draws.
  - **Development sites:** the Table 5 draws (REV-1 m = 22, arm a2, seeds 41–45, 3 sites), via `m3_replay_lib.replay`.
- **Gate B (exit 3).**
  - Simulation:
    - (i) every S3 row of every task equals the stored main-study S3 row. That is 680,000 rows.
    - (ii) the anchor tasks (seed 41 of e10:0, e10:3, e12:5, e12:14, e12:19) also rerun CB, GEE and BB through `simw3_baselines.make_extra`, and every row equals simw3_rows.json (bed0438a).
  - Development sites:
    - (i) `m3_replay_lib.self_check`: S3 and B1u equal REV-1 within 1e-12.
    - (ii) all 36 CB, GEE and BB cells (c0 and c1s) equal w3.json (2f015508) within 1e-12, on 8 fields.
- **Reported, no prediction.** Coverage, MCSE, cannot-certify, ΔSens_ε and bound exceptions. A cell with exceptions above 1 % of its draws is flagged, as in W3.

### Part C: tuning constants (`p5a_s2_sim.py --part C`, `p5a_s2_dev.py --part C`)
- **One at a time around the defaults.** n_s ∈ {10, 40}, r_low ∈ {0.60, 0.80}, G ∈ {100, 400}, plus the default.
  - n_s and r_low are patched into `m3_replay_lib.cstat`. At the defaults, the original function object is kept.
  - G is set through `m3_replay_lib.G`.
  - m = 22, so every variant uses the main-study draws. Variants are paired.
- **Simulation:** e10:0, e12:19 and e12:10, with 50 realisations.
- **Development sites:** DE-Grid and single-candidate B1u on the Table 5 draws.
  - Split-Select is not run on the development sites. It needs the real-score Platt path of `m3_split_dev_eval.main`, which is outside S2.
- **Gate C (exit 3).**
  - The default variant equals the stored rows, as gate A in the simulation and as self_check on the development sites.
  - The patched centring at (20, 0.7) equals the frozen `cstat` on fixed test vectors.
- **Prediction C.**
  - Coverage of both rules stays ≥ 1 − δ − 2·MCSE in every simulation cell, at every variant, at ε 0.10.
  - Sensitivity and cannot-certify may move.
  - The development sites are descriptive. Five seeds give no usable MCSE for a judgement.

## 2. Run (VM only, after review and 实验管理's OK)
- **Wrapper.** `run_p5a_s2.sh <sim|dev> <A|B|C> <tag> [workers]`.
  - Package `p5s2_pkg_<tag>`, checked with sha256sum -c. It holds results_src copied from p5fifty_pkg_20261124T211616Z plus sizefpr, s2/ with the three scripts, and inputs/w3.json.
  - Busy gate on every <DATA>/*/RUN.log, on any <USER> GPU process and on MemAvailable ≥ 32 GB.
  - Job directory <DATA>/p5s2_<kind><part>_<tag>, with RUN.log START/END and PIDS.
  - nice 19, taskset 56–63, 8 workers, 1 BLAS thread.
  - Watchdog: tree RSS 24 GB, wall 4 h.
- **Order.** One job at a time: A sim, B sim, B dev, C sim, C dev.
- **Estimated time at 8 cores.** Pass-1 tasks took 15 s at m = 22 (fifty_b1 E10 replay log).
  - A: 1,200 tasks, ≈ 1–1.5 h. Time grows roughly with m, a mean factor of 1.4.
  - B sim: 1,700 tasks at ≈ 20 s, plus 5 anchor tasks at ≈ 720 s (BB). ≈ 1.5 h.
  - B dev: 15 tasks, dominated by BB. ≈ 0.5 h.
  - C sim: 1,050 tasks, ≈ 0.7 h.
  - C dev: 105 tasks, ≈ 0.2 h.
  - Total ≈ 4.5 h.
- **Write volume.** A ≈ 30 MB, B ≈ 60 MB, C ≈ 25 MB, dev < 5 MB. Under 150 MB in all.
- **Pull.** Through the scratchpad, checking SHA256SUMS.

## 3. Entry into the paper (fixed now)
- **Results master r2 for P5a.** An add-only builder adds key families `s2A/…`, `s2B/…` and `s2C/…` from the summaries. Each key records the source file sha and the path.
- **Main text.** One compact table or figure, plus one sentence per part.
  - A: coverage against m for both rules, in the 4 cells.
  - B: one row pair, CR2 and WCB, next to the existing baselines.
  - C: the range of coverage over the variants.
- **Details go to the SI.**
- A failed prediction is reported as a measured limit. The rules are not changed.

## 4. Known and not blocking
- WCB uses the unrestricted wild bootstrap-t with the CR2 SE. The restricted (null-imposed) form is not used, because the fixed-sequence test needs an upper bound at each grid point, not a test at one null.
- The m = 10 Split-Select certification group has only 4 fields. Exceptions there are counted and flagged, not hidden.


## r1 note (2026-11-28T01:52:39Z, written before the part A rerun; 实验管理 approved sim r1 1b9b52e6)
- **Sampler capacity.** In e12:19 the two-stage sampler draws s = 2 of S = 12 slides of exactly F // S = 25 fields. At m = 45 the two-stage design samples 45 of the 50 fields of the two drawn slides; at m = 60 it is not applicable. Prediction A(i) for e12:19 therefore covers m ∈ {22, 33, 45}. The other three cells keep m ∈ {22, 33, 45, 60}. The slide count is not scaled with m, because that would change the cell design.
- **Run log.**
  - Part A run 1, TAG 20261128T012737Z (package cbd71242): FAILED, rc 1, at task 543 of 1,200 (2026-11-28T01:29:16Z to 01:50:58Z). Cause: the sampler capacity above (ValueError in e12_stress.two_stage_sampler at e12:19, m = 60). No output was written and no gate was reached. The job directory <DATA>/p5s2_simA_20261128T012737Z is kept.
  - Part A run 2: under a new TAG, with p5a_s2_sim.py r1 (feasibility skip and capacity check only).


## r2 note (2026-11-28T05:56:03Z, written before the part A rerun; 实验管理 approved common r1 0b02fb98 and sim r2 fb28a71a)
- **Run log.**
  - Part A run 2, TAG 20261128T015239Z (package 6ab05e85): FAILED, rc 143. The watchdog killed it at the 4 h wall cap, in the gate step, after all 1,150 tasks had finished (the last task line was at 02:35:52Z; END at 05:53:35Z). No output was written and no gate failed. The job directory <DATA>/p5s2_simA_20261128T015239Z is kept.
  - Cause: stored_rows read the stored npz row by row, and NpzFile decompresses a whole array on every item access. Fix: common r1 reads each array once, and sim r2 loads the stored rows before the pool. Tasks, rules, gate conditions and outputs are unchanged.
  - Part A run 3: under a new TAG, with common r1 and sim r2.
