# 附录 S0′（封存）：前端 v3 的常数 s0′、r0′（PREREG_P5_A2PPP_SAMPLERBC_20261028_sealed.md，3f980c7a；勘误 E1a e2507f12、E2）

- **封存批准**：用户原话"批准 E2 并封存 S0′，继续"，2026-10-29T00:10Z（实验管理转达，date -u）。封存在新 M0 的任何结果之前。此后常数不再改动。
- **定义**（预注册 §2.1、§4.2，勘误 E1a）：
  - NIH 训练折的全部 case（split_cases_greedy(seed = 42) 的 case_id；不抽样，不用标签）；
  - 每个 case 用 v3（code/p5stage_a2ppp/sample_rbc_iqr_prod.py a7ec0817）的 prepare 和 sample_stats（TAU = 0.75）得到 ref_case、s_case；
  - s0′_c = s_case,c 按 case 的中位数；r0′_c = ref_case,c 按 case 的中位数。
  - 事先写定的附加规则（实验管理审定，台账第 205 行）：整样本回退的 case 不进入中位数并计数；若超过训练折 case 的 1%，标 stop_sample_fallback，不给常数。

## 1. 数值

| 通道 | s0′ | r0′ |
|---|---|---|
| R | 0.01861557364463806 | 0.11157477647066116 |
| G | 0.030933722853660583 | 0.20958289504051208 |
| B | 0.01961684226989746 | 0.20958289504051208 |

G、B 两个通道的 r0′ 相同，原因是 OD 由 uint8 换算而来、取值离散，两个中位数恰好落在同一个离散 OD 值上。

## 2. 来源

- 运行目录：runs/a2ppp_s0p_20261028T230649Z（VM 绝对路径 <VM_HOME>/P5_StainRef/runs/a2ppp_s0p_20261028T230649Z/）；旁文件在 runs/a2ppp_s0p_20261028T230649Z_side/。
- **s0p.json sha256 21e76a3f7a1b06720b94f59149d01c2710008d8274e3eecc51392df113192f17**（在 VM、scratchpad、OneDrive 三处一致）；SHA256SUMS 003877c5df2b2e89ff9f79d87e9883a13324d6a85fd606b05e5ce6878e7dea1d。
- 脚本 results_src/a2ppp/a2ppp_s0p_gate.py b58dd65c；包装 run_a2ppp_s0p.sh 5f41b1da；台账第 204–207 行。
- 运行 2026-10-28T23:06:51Z–23:07:22Z，exit 0，stderr 0 字节。
- NIH 训练折 19,338 张（阳 9,719），140 个 case；验证折 4,081 张、30 个 case；测试折 4,139 张、30 个 case（只计数，未打开）。三折 case 互不相交。
- **s0p_status = ok**：整样本回退 case 0 个，回退裁剪 0 张；全部 19,338 张走 NIH 分支；s_min 触发 0 次。

## 3. 只报告的读数（不影响常数）

- s_case 分位数（q10/q25/q50/q75/q90）：R 0.0103/0.0141/0.0186/0.0275/0.0350；G 0.0219/0.0268/0.0309/0.0382/0.0516；B 0.0133/0.0164/0.0196/0.0232/0.0285。
- ref_case 分位数：R 0.0797/0.0897/0.1116/0.1560/0.1855；G 0.1725/0.1901/0.2096/0.2264/0.2439；B 0.1387/0.1482/0.2096/0.2243/0.2321。
- 截尾后 s 与截尾前 s 之比的中位数：0.828 / 0.813 / 0.808；|R_S|/|F_S| 为 0.750–0.752。
- 每个 case 的图像数：中位数 86（q10 68，q90 310）。

## 4. 捷径关卡的结论

- **原关卡（预注册 §4.3）**：stop。感染 107 / 未感染 33 个 case；A* = 0.803（ref_B）。原样保留。
- **勘误 E2 下（逐细胞、相对 legacy、余量 +0.05）**：A*_cell(v3) = 0.5687，A*_cell(legacy) = 0.6094，Δ = −0.0406 ≤ +0.05 → **continue**（runs/a2ppp_ds2_20261028T232106Z，台账第 220 行）。

## 5. 代码中的封存

code/p5stage_a2ppp/sample_rbc_iqr_prod.py 中填入：
- S0P_JSON = "<VM_HOME>/P5_StainRef/runs/a2ppp_s0p_20261028T230649Z/s0p.json"
- S0P_JSON_SHA = "21e76a3f7a1b06720b94f59149d01c2710008d8274e3eecc51392df113192f17"

只改这两行。v3 的新 sha 和 VM 上 sealed_s0p() 的只读校验结果记在台账。此后所有脚本固定新的 v3 sha。
