# P5 protocol freeze (2026-10-22)

Decision (user, 2026-10-22, via the P4 session): freeze the P5 protocol at **m = 22 fully annotated fields, rule S3**. Nothing below may change without a new dated freeze file. The reference implementation is `s11.py` (rule `S3`); the numbers in the paper are read from its result file.

## Rule S3 (as implemented in s11.py)
Input: a deployed classifier with scores s = logit(positive) − logit(negative); a target site with fields (images); an unlabeled score pool.

1. **Field sample.** Draw m = 22 fields uniformly at random without replacement from the site's fields (fields = images; units as in Methods). Annotate every cell of these fields: positives (clicks) and negatives (all remaining cells are labeled negative).
2. **Centering candidates (all unlabeled, per field i, n_i = total cells in the field, positives and negatives):**
   - c0: c_i = 0.
   - c1s: c_i = g + w_i (median_i − g).
   - c2s: c_i = g + w_i (c2_i − g), c2_i = median of the lowest ⌈0.7 n_i⌉ scores of the field.
   - w_i = n_i/(n_i + 20) (n0 = 20 fixed a priori); g = median of all cells of the site's evaluation pool.
   - Deployed rule for candidate k: flag cell j of field i as positive iff s_ij − c_i^k ≥ τ_k.
3. **Certification per candidate (B1u at level δ/3, δ = 0.10, K = 3).** Threshold grid: up to 200 distinct values, the pooled labeled negative centered scores at the j/200 quantiles, in descending order from the maximum. For each grid value τ: X_i(τ) = number of labeled negatives of field i with centered score ≥ τ; N = Σ n_i^neg. Compute
   - MSB = Σ n_i (p̂_i − p̂)² /(m−1), MSW = Σ n_i p̂_i (1 − p̂_i)/(N − m), p̂_i = X_i/n_i, p̂ = ΣX_i/N, n0 = (N − Σn_i²/N)/(m − 1); F = MSB/MSW (MSW = 0: ICC_U = 0 if MSB = 0 else 1);
   - θ_U = (F/F_lo − 1)/n0 with F_lo the δ/K quantile of F(m−1, N−m); ICC_U = θ_U/(1+θ_U) truncated to [0,1];
   - n̄_w = Σ n_i²/N, n_eff = N/(1 + (n̄_w − 1) ICC_U), x_eff = ΣX_i · n_eff/N;
   - upper bound UB = Beta^{-1}(1 − δ/K; x_eff + 1, n_eff − x_eff).
   Fixed-sequence (Learn-then-Test) over the descending grid: certify τ while UB ≤ ε; stop at the first grid value that fails; the candidate's threshold is the last certified (lowest) one. If the top grid value is not certified, the candidate is not certified.
4. **Selection.** Among certified candidates, pick the one with the highest sensitivity on the labeled positives of the same 22 fields (centered with the same candidate); ties go to c0, then c1s, then c2s. No labeled positives: ties. No certified candidate: output c0's highest grid value and report "cannot certify".
5. **Guarantee statement.** Each candidate is certified at δ/3, so if each B1u bound holds then P(site FPR of the chosen rule > ε) ≤ δ by the union bound, whatever the selection. B1u is an asymptotic/approximate bound (F-interval ICC upper limit for unequal cluster sizes, Clopper–Pearson on a continuous effective count); the guarantee is therefore **not finite-sample**; coverage is evaluated empirically.
6. **Targets.** (ε, δ) = (0.10, 0.10) main, (0.05, 0.10) secondary. FPR is cell-weighted over all true negatives of the site, computed on the chosen candidate's centered scores.

## Frozen evaluation setting
Classifier arms A2 (main) and A0-ctrl (secondary), seeds 41–45, 200 draws per seed; sites BBBC-ALL, MP-IDB, Lacuna-B; reference sensitivities at the site-level FPR = 0.05 and FPR = ε thresholds; baseline R1-2 (two random labeled negatives per field, R1 on 2m scores).

## Provenance (sha256)
Scripts on the A6000 VM, `~/P5_StainRef/iter1_20261021/code/`:
- s11.py c7d3ad6228254855dde16d08087b47d911558cd33869d12e769617daa36fca6f
- s10.py 6d60574edc60aea14f599d40ae291575c0f72c469a7c9bb8be29b2a5a535acdd; s10_bounds.py 6037e0fe08b86f9ef2be7ca42748ea2f36a2f7ff2d071e17648b3aefec295090
- s9.py b3030ac2cf4aee3d71d66623892ae7783f46bd9b37a0a3450b131cbbd530315d
- s8.py bbd2cf998cfcfbd6ef9c77bfa303ba52f23d714ea2f81d439578f9488dbe2247; s8_bounds.py 54a2af8d3704c16b0994dcdc584d4f407cd9564cc07d49ad74cdf7c4cad5ae7f
- s7.py 04be30707cc28eec6b678cf6f133bdb5ff672f06a60b438bc7de19e1d0860d24; s6.py 293da4c3d230d1a51f7b1a8616240a722b89fe41e0cba125b1b700b831508f79; s5.py 9c6f8e83eca6947f40695ce3aff924431c43d97f8fc352bc224a79e073a34ce3; s4.py 967aa81fc8be151fa9398686e06d233b09a26c06766c38aede0a1888427f000f
- s2_mb.py e9d11f67369a8823c2f62e777fa6d5ce9a3af5644b5dc6f7c1506bcef009328d; s2b.py ecfb64e3760a100934ea68bf6f9204c6ff2cb3e55e088d6892e774f6c1bd5422; s2c.py 0b7d74291ba58d3abf6a9361e199c3285a4facf4e911ea38d56813da0eb849e7
- p5_eval_phase1_iter1.py edf1b8fb70f97e8d6cce6e308f7736dabdae8f7cf80b8e0173b4ec1d460ff3dd
- Lacuna (`<DATA>/p5_e7_lacuna/logs/`): e7_crop_v2.py 5432a2fe4e78668388dff5143863ac2f8c2f38905700d6561e0242aabaec126a; e7_infer.py 9c666a9f08b42c69a7bd1f03779a08dcc52c8aa7c7ac07d683bda6b05d7f92f5; e7_analysis.py 3831fea92c1a6a4ce9d7ceab573f4238634dfc0a44dac1f2b6f87feb6f0193bc; e7_s2c.py 96f2e6a5bd7a16e854400fba6d042caceeb655d72f91649ab7da388b74dda8ad

Preregistrations (sha256 of the .md file): ITER1 143a07d47830a968a4955adbc5a5f42308ec40f977fb1d2302ea499cda5a4a4d; ITER1_ADDENDUM_S3S2b 266da7fb82b67a48f6f509f71434291c93d25fad9125728a7a18c99a9fb08e3a; ITER2_S2c 0dda7bb0a7bb14e9655547f8144e335bbe0b67066dd73ccd3a2b4afb76a00db9; ITER3_S2d ecd4f6ec2bb4d6510493d6873d75b232f23450cd5055e2d08730d3d324f6e934; ITER4_S4 81f95c590b889ceeb93eef0d0084b381fe557eedc0dec35596d667f0db672408; ITER5_S5 0bdfc41bce403b968b823b13abb110c21c8afa3ce8ef70f43553f77f444693c3; ITER5_DEFECT a505afb0a191028e7455e226a9ba663dda655235c5e59f32cdc955966d7d96ca; ITER6_S6 e0191de52d6ee154844cf5694a506ed9b869df808362298ee0fa89d62179942a; ITER7_S7 338e83229b684b299b41e4cb029d64ef6073a7a015f46e5ed2577c858dc2e04e; ITER8_S8 c8037e9a1764606f69205f7f5b16d69c5e7221529724c6c779f86096549126e0; ITER9_S9 e0aa6583479ffba5682fbdf373553eae9d4a69761d0350351608ac82f3325b80; ITER10_S10 986966c96aec82953a07058f240770953fa848c84766d1ce1bc2089f92a1a160; ITER11_S11 f51220c78630006f0d44ba58fd11dee7940d5ac479c5e30b723d8ceca2fea851; ITER11_DEFECT 83e96bba2708ca58a5efbc18ef2996aa2513fd35dca86c6910605e138c2e3629; E7_LACUNA da4297a54cfb1b32b821a826bf24c6cbeb8c4dcf5f42889663b4566f8a22feb8; E7 ADDENDUM_DETECTOR 66a048f47a3fbf96bb7c0e6a5e8c3d274166bc06971f13c4e4ed0e568790c8fb; E7 ADDENDUM_B_ONLY d0613641fc7edf03d254daf2f99854932a4687d5a3aab0b0432bdcb5396044c4.
