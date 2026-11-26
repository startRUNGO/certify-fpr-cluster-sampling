# P5 E8 预注册（替代 Lacuna Ghana）：PlasmoCount 测试集一次性确认评估

写于 2026-10-23，早于下载任何图像、早于查看任何图像或模型得分。设计：P4 会话；执行口径补全：P5 会话（标“执行口径”）。用户在 P5 会话中的明确批准：只下载 Mendeley 记录 10.17632/j55fyhtxn4.2 中 test.json 引用、记录里可得的 116 张 TIFF 到服务器 <DATA>/p5_e8_candidate/images/，核对字节数与记录里的 sha256；先写本附录并记录 sha256，再跑确认评估；nice 10；不碰 ~/P4_rerun、~/P4_samemod。取代 PREREG_E8_GHANA_20261022.md（sha256 e64656e1…7902；Ghana 因类别 id 映射缺失按其停止规则停止，未裁剪、未打分）。
冻结物（不得改动）：PROTOCOL_FREEZE_20261022.md（sha256 f2565f309fff346cd8056ad2e543bea9bc3031041da7ce7be54b5dbb5aa0802b）中的规则 S3 与评估设置；s11.py（c7d3ad62…fca6f）；裁剪的 PAD=5 与缩放到 128×128 的推理流程（e7_infer.py 9c666a9f…2f5）；成功判据同 S11（m=22）。PlasmoCount 未参与规则开发（11 轮迭代只用了 BBBC-ALL、MP-IDB、Lacuna-B）。

## 1. 数据与评估集
- 来源：Mendeley Data 10.17632/j55fyhtxn4.2（Davidson & Baum，CC BY 4.0）。标注 test.json（sha256 984e26ccb72178e03ef5f332053dd2073089dd8631e5c7e6bc80f1b0248517fd，1,045,001 B）引用 121 张图；记录里有 116 张（5 张缺失，列于清单）。**评估集 = 这 116 张图**（视野 = 图像）；training.json 与 *_stages.json 不用。下载清单 e8_plasmo_manifest.json（sha256 c74ecb6d3cb6d3376671df538875dd5ce827fb0947e7d67eb08874e897bfbe8f；含每个文件的字节数与记录 sha256，合计 4,059,721,430 B）。逐文件核对字节数与 sha256，任一不符则重下，仍不符则停止报告。
- 标注（已在 JSON 里，打分前已知）：每个框有类别 uninfected / infected / not_sure；框为 {minimum:{r,c}, maximum:{r,c}}（r=行，c=列，像素，图像坐标）。**阳性 = infected；阴性 = uninfected（已标注，不用检测器）；not_sure 排除（既不当阳性也不当阴性）。**
- 裁剪：与 E7 相同，框每边外扩 PAD=5 像素并截断于图像边界；[r_min−5, r_max+5] × [c_min−5, c_max+5]。执行口径（读图）：用 PIL 打开 TIFF，若为 8 位 RGB（或可 convert("RGB") 的 8 位）则直接转 RGB；若图像为 16 位或其它非 8 位模式，则在看任何得分前停止并报告，不自行定缩放规则。若 TIFF 含多页只取第一页。
- 门槛（图像和标签层面，在 JSON 上即可判定）：(a) 每张评估图的 uninfected 数的中位数 ≥ 5（JSON 显示 test 全体最小 32，已满足）；(b) 至少有一张图含 infected。检测器召回门槛不适用（没有检测器）。
- 注意：标本类型（体外培养或临床）、物种、成像设备在记录页没有说明；视野=图像，没有玻片或患者 ID；这些写入局限。

## 2. 打分前报告（读图之后、推理之前）
报告：成功读到的图像数、图像尺寸/模式分布、各类框数（infected、uninfected、not_sure）、裁剪后阳性与阴性数、每图 uninfected 数的最小/中位/最大。若与 JSON 计数不符（框越界等）如实报告。

## 3. 推理与规则
推理 35 个检查点（批量 512；同时最多 1 个；显存总计 ≤30 GB，给 gaoqihuang 留 ≥8 GB；nice 10；CPU 进程 ≤4，且不碰 ~/P4_rerun 与 ~/P4_samemod），保存逐样本 logit。然后**一次性**运行冻结规则 S3（m=22；(ε,δ)=(0.10,0.10) 与 (0.05,0.10)；A2 为主、A0-ctrl 为次；5 种子 × 200 次抽样；随机数 numpy default_rng([20261022, 13, 臂序号, 种子, 22, 0])；候选 c0/c1s/c2s，δ/3 的 B1u，固定序列，选择规则同冻结文件；视野 = 图像；“视野内全部细胞”指该图的全部已标注细胞（含 not_sure？——**不含**：中心化统计 c_i 只用该图 infected 与 uninfected 细胞的得分，not_sure 不参与得分计算也不推理）。报告：覆盖率；灵敏度减 c0 的 FPR=0.05 参照；基线 R1-2（m=22，每视野 2 个阴性）；源阈值 FPR 调查（35 个检查点，源阈值 = NIH 验证折 Youden）与 AUROC–FPR 秩相关；cannot-certify 率；标注量（点击=infected 数、已标注阴性数）。标签噪声：阴性是人工标注，主口径 η=0；稳健性并列 η=0.10 下界 max(0, FPR−0.10)（与 Lacuna-B 对齐，稳健 = 两种口径都成立）。S3 对每个臂和设定只运行一次；看到数字后不改 S3。执行口径：评估代码只允许改数据加载（JSON→裁剪、路径、站点名），差异存档并记录 sha256。

## 4. 成功判据（与 S11 在 m=22 时相同；A2，(0.10,0.10)）
覆盖率 ≥ 0.88，且灵敏度代价（c0 的 FPR=0.05 参照灵敏度 − 实现灵敏度）不超过 0.05（即灵敏度差 ≥ −0.05），或不超过 R1-2 在 m=22 时代价的一半。两个条件都满足为“通过”，否则“未通过”，如实报告；η=0 为判定口径，η=0.10 并列。

## 5. 纪律
不覆盖已有结果；新结果写入 ~/P5_StainRef/runs/e8_plasmo/ 与 <DATA>/p5_e8_candidate/；异常先停下报告；每步完成后向 P4 会话简报。
