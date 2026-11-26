# P5 W1 预注册：M5 HCM-1000x 外部确认站点（一次性确认评估）——定稿稿 draft10（2026-10-31；相对 draft9 80e2950e 按勘误 E-M5-2、E-M5-3 改执行位置：裁剪拷到 VM 只读保存，VM 只推理、评估回 Mac，台账分推理与判读两行；规则、判据、数据、模型、抽样流不变；draft9：相对 draft8 6ccdfe7a 改推理路线为原始 png 字节流（P1）并写入更正；draft8：相对 draft7 9b2dbb32 只更新部署包（S2 定长批次后重建）与第 0.4 节一句；draft7：相对 draft6 0ca60852 改第 0.4、3、5、6 节的推理路线；draft6 相对 draft5 8dd26dc0 补第 3 节第 6 条的措辞；draft5 相对 draft4 e2b4f2e6 补第 0.4、3、6 节；待实验管理审查，再由用户批准封存；未运行任何分类器）

- **由来**：草案 3（7d21a950）+ 附注 A（57c6c59c）+ 勘误 E-M5-1（8f914577）+ 构造结果（check 阶段与 full 阶段），以及实验管理 2026-10-30 的两点补充（推理与判读路径写死；单站判读）。本稿取代草案 3 作为封存对象；草案 3 及其附注保留，不改。
- **性质**：对开发站上冻结的 SPLIT6 的**唯一**确认。PlasmoCount 和 IML 是按 S3 预注册并已读过的，对 SPLIT6 不算确认证据。
- **已知劣势（必须在本预注册和论文中如实写明）**：开发站上 SPLIT6 没有完全通过 draft7.1 的判据。旧版 MP-IDB 在 ε 0.05 时比 S3 低 0.0163，超过 0.01 的限值；ε 0.10 时为 −0.0003。其余 6 格都优于 S3，重建三站在 ε 0.10 的平均改善为 +0.0543（sp.json 364ebc08）。用户在知道这一劣势的前提下决定带 SPLIT6 进入确认。
- **封存须用户批准。推理脚本在用户封存之后才能运行，只运行一次。**
- 调用的 skill：light-research-plan、light-self-review。

## 0. 冻结物

### 0.1 规则（同一批抽样上同时评估）
- 与草案 3 第 0.1 节相同，不改：S3（REV-1 定义）、c1s（单候选 B1u-c1s，水平 δ）、SPLIT6（draft7.1，m_s = 6，(C1) + (C2)，每个候选 δ/3，三候选 union；K1 / K1-v 的 bootstrap-t B = 999；退化处理按 draft5 附注 A；Platt q 冻结、不重新拟合）。

### 0.2 脚本、文档与开发站选择记录（shasum 原始输出）
```
cbf1e7f4f33dd1c74cbfd0f3ac4d34a397e8a7ef44b82c9365ca117733bf6457  results_src/m3/split/m3_split_dev_eval.py
ed9d1f17701c9f30e92df0152a1c10bec6f7ee2d2bc392e65d617962ff0081de  results_src/m3/k2prep/m3_replay_lib.py
0091a3337074a227e404f0585966b87a9f5e9e55cfbb8dd766df7e4618893b67  results_src/m3/k2prep/m3_k12_dev_eval.py
75eb9afd3a130b00170757c4027f8b158bf95deb0ad6124229e9292c60e2d1ad  results_src/m3/stress/m3_stress.py
7b7e6ecda22b789c573451c757b6e197d747f77e10bcfc201ccdb85921e25abd  results_src/m3/k1vb/m3_k1vb_dev_eval.py
0e6aa9f8f423dbaed286db969afa8a14c9ce477e37bba962732d3b411002147c  results_src/m3/k1vb/m3_k1v_precheck.py
47f79f1df632e5efd1be5ea9485826d082eeaac38f89ee6f60a6aed2292116c5  results_src/mh/mh_rev1.py
d2c5a1827681066e6a4435f3929fa767773aeb35f7427267abc09af5813a4887  P5_机制创新_问题卡_先选后验_20261030_draft7.1.md
86a4c9e4b016f3387d3c33d38c6efef53d8a4d120074340f9d5bbe14499f5475  P5_理论草稿_K6_选择条件有效性_20261030_draft2.md
a876bf8f3793f21a05cd9672b85378d45626e517966036aad8c65c9872df3470  P5_理论草稿_K6_选择条件有效性_20261030_draft2_附注A.md
364ebc0819f1ac80d5db7a9ef373d08815591f592f9fed93b5676b7564288ec7  runs/m3_split_dev_20261030T170958Z/sp.json
a67ec5febadbaf1fa97514a6014c6ed86ead9699d4a5d368d0c67c871374ceef  runs/m3_splitr_dev_20261030T173014Z/spr.json
dca9ffa0b8bb8a4ab7647389edd517563ae711d1a552e60b59f7130aea8d3aee  P5_机制线决策材料_20261030.md
7d21a9508b995044d1b95b6121db5bd92c3cfacf645066e4bd164b5fcb64efb0  PREREG_P5_W1_M5_20261030_draft3.md
57c6c59cac7b518c701c840b96b1b7ab1812336495b0bd991d1261a349dfc792  PREREG_P5_W1_M5_20261030_draft3_附注A.md
8f9145773c33c04e7576aa053063afdc539bb91259b1b8685e61b9be6a5d0ddb  勘误_E-M5-1_M5构造尺度门槛_20261030.md
```
- Platt 参数（q = expit(a·s + b)；由 sp.json 的 platt 字段原样输出；推理脚本在 NIH 验证折上重算一次，必须与下列数值逐位相等，否则停止）：
```
seed 41: a = 1.3435275318401698, b = 0.7925587857445654
seed 42: a = 1.3285637311491154, b = 0.6609026438276017
seed 43: a = 1.2593773240565085, b = 0.7324933734980241
seed 44: a = 1.290789659587087, b = 0.7766941603495804
seed 45: a = 1.2591539840215842, b = 0.5389064779958056
```
- NIH 验证折文件（只读 logits 与 y 两个键）：
```
8e7fcd34d0dc119179dab406898ca9d08604e866885b6da83f5c262923b6cbd8  runs/c1a6_stage2_20261027T133543Z/a2ref/a2ref_s41_nih_validation_logits.npz
6d85e599802609f7a33dedc471a53f58d52bf4d1dbe52fb9b3ca94fa9fbab79f  runs/c1a6_stage2_20261027T133543Z/a2ref/a2ref_s42_nih_validation_logits.npz
71cd404368a600ad7a3a5470da93c0a7a83791e3173bdca14da8d06e6f3d886b  runs/c1a6_stage2_20261027T133543Z/a2ref/a2ref_s43_nih_validation_logits.npz
4b046baf7967608712d98f7fd30d3cdc1f6a0926c4014d6d8febdb43eaaa427f  runs/c1a6_stage2_20261027T133543Z/a2ref/a2ref_s44_nih_validation_logits.npz
7d71e4e47451fe32eb54530470176ddd05194330ac824df2eff458c77edaab66  runs/c1a6_stage2_20261027T133543Z/a2ref/a2ref_s45_nih_validation_logits.npz
7780177f854dc31b0122d4bba45dbf79c2c735d5f243e54f3954abf457039d62  runs/c1a6_stage2_20261027T133543Z/a2ref/ids_nih_validation.npz
```

### 0.3 构造产物（已完成；构造阶段不涉及任何模型输出）
```
dde8c9bad1fd7e5d9e5d87dbc78f4d217c43fd796ee314c0330fbf0e26c32142  results_src/m5/m5_construct.py
5432a2fe4e78668388dff5143863ac2f8c2f38905700d6561e0242aabaec126a  results_src/m5/vendor/e7_crop_v2.py
b5cdd8867b01e7d756fbd8ec13233365451e04e1af66744c169928453d024a40  <LOCAL_PATH>
a774649dbdb3cac65981e8f17f7fb984c57b72a0e457d3186d0bd75e0f26cc47  <LOCAL_PATH>
07a8204943f3d38bd23ad28cafc09d088ab9b2c9d14a4badd564f0398f2b874f  <LOCAL_PATH>
47619ccaf467d1fbd66a4f61ec15835af5acb0d558e12357ea7e7260039877e4  <LOCAL_PATH>
013f33cee606e7ea23cef8976ebb0f9ad6e256be7c446615cc40ce1b17d0b865  <LOCAL_PATH>
01659754186b31cbd27bfd15b9e95984e5e6b5453cb85e1f1b037021d9168b28  <LOCAL_PATH>
```

### 0.4 分类器与推理前端
- A2 的 5 个冻结检查点（种子 41–45）。文件 sha 固定为（来自 results_src/c1a6/f1_nihbg.py 的 CKPT_SHAS，即 C6s 运行的 c1a6_infer_meta.json）：
```
049aade5ce6481fc5f7c7bf63e09d33ceb4e2cb190f5f050451509f98578f1f1  a2_s41.pth
9d56a35ce40941d0c52467b0a3797301e39026e663b330af3c088c487e3e4014  a2_s42.pth
b0125a2a70282adacc4fa76e2804d90ccba7ad66591a9298091d4ec7b62b538b  a2_s43.pth
2da90395266ddc10c3e525387149a0fce4d5f7c3d83c754bc6fa2ac891c936ab  a2_s44.pth
f4cb22539b8ee452a8279abdd1fc3e03e2bb72466af7105210fb067c700cb7b2  a2_s45.pth
```
  （以上 5 行从 f1_nihbg.py 第 46–47 行原样转录；推理脚本以文件 sha 断言，不一致即停止。）
- 模型代码 malarianet.py：sha c51cc13deffd1a31cffd6af3dda1636dd5353f2f806baf3545d3fd709aff56f1（VM 上 C6s 运行的 model_py_sha256，同一来源）。
- 推理前端与 R8 / mh_infer.py（594e3549）相同：A2 臂用 p5stage.prod.legacy_float_pr_od（e638f61e），之后 PIL BILINEAR 缩放到 128 × 128，ImageNet 均值方差归一化，batch 512；分数 s = logit[1] − logit[0]。
- **推理路线（实验管理 2026-10-31 决定 P1）**：把 M5 裁剪的**原始 png 字节**按 manifest 顺序、每帧至多 512 个（PNG1 帧，每个文件带长度与 sha256），送到 VM 上的本地进程（标准输入；服务端为 bundle3c 里的 m5_stream.py serve，由推理步以精确路径启动）。VM 进程逐个核对 sha，在内存中解码并用**逐字复制的 IL.load_one**（参数为 (io.BytesIO(原始字节), True)）做前处理，补零到 512（每次推理调用的输入形状都断言为 (512, 3, 128, 128)），用存档管线的同一设置在 GPU 上推理——**torch 默认设置（TF32 开关不改动）**、cudnn.benchmark 为 False、确定性算法关闭——只回传前 k 行的 logits 和每个文件的 sha。**（E-M5-2 改）M5 裁剪事先拷到 VM 的 <DATA>/p5_m5crops_20261031T162417Z，只读，结束后保留到论文完成；拷贝逐文件 sha 在 Mac 与 VM 两端核对一致后才使用；VM 端除这份只读拷贝外只写运行日志与推理输出目录，服务端进程不写任何 M5 文件，也不写 /tmp。**客户端核对服务端回显的每个文件 sha，并把它们写入 crop_sha256.txt。这样解码、legacy_float_pr_od、缩放和前向计算全部与存档管线在同一台机器、同一环境中完成。
- **传输（E-M5-2 改）**：M5 HCM 评估位置的原始裁剪约 7.7 GB（183,105 个文件）由实验管理用 24 路 tar-over-ssh 事先拷到 VM（原单路 ssh 流式约 18–44 KB/s，估计要约 4 天，已放弃）；封存前由实验管理逐文件核对 Mac 与 VM 两端 sha，全部一致后才封存；Mac 端清单 <LOCAL_PATH> 893c53ace3290d85da1b7d825ae8fd89b86c28a40d6264a6186ed1b621952bed（183,107 行，含 manifest.csv 与 full_report.json；评估步 --mac_crop_sha 的输入）。推理步从 VM 本地读取这些文件。
- **更正与理由**：设备诊断（m5_devdiag_compare.json 4ec283ca）逐阶段只核查了 10 个裁剪，它们恰好逐位相同，当时据此认为差异来自 torch 的 CPU 前向——这一结论过强。随后的流式设备检查（m5_stream_devcheck.json 735c9fd7）让前向计算在 VM GPU 上与存档完全相同、只把前处理放在 Mac 上，差异仍为 0.018–0.037（与 Mac-CPU 几乎相同），补零已由 smoke 排除；所以差异的主要来源是 Mac 上的前处理（Pillow / numpy 版本在其余裁剪上的差别）。P1 让前处理也在 VM 上进行，预期差为 0。
- 部署包 bundle3c（VM 端运行 m5_stream.py serve；shasum 原始输出）：
```
45be135999b658fec2600f0b8949ebe582b35ced54d336b40c4c0da851fe3140  m5_bundle3c.tar
6d9dd322085a5ea349ea5eb5fb918136bbf9ab1e7eb047a63d25cef2c78d7dec  m5_bundle3c_SHA256SUMS
d3ac4d782568eceba4fe53455b200608ed4fd7579ae3a455882ea5ade72cfc90  results_src/m5/m5_infer_lib.py
9a0863360c701157eaf6959e52537f1d310d1fbcd7bc93557716d9c50b79c535  results_src/m5/m5_stream.py
049aade5ce6481fc5f7c7bf63e09d33ceb4e2cb190f5f050451509f98578f1f1  results_src/m5/vendor/a2_s41.pth
9d56a35ce40941d0c52467b0a3797301e39026e663b330af3c088c487e3e4014  results_src/m5/vendor/a2_s42.pth
b0125a2a70282adacc4fa76e2804d90ccba7ad66591a9298091d4ec7b62b538b  results_src/m5/vendor/a2_s43.pth
2da90395266ddc10c3e525387149a0fce4d5f7c3d83c754bc6fa2ac891c936ab  results_src/m5/vendor/a2_s44.pth
f4cb22539b8ee452a8279abdd1fc3e03e2bb72466af7105210fb067c700cb7b2  results_src/m5/vendor/a2_s45.pth
c51cc13deffd1a31cffd6af3dda1636dd5353f2f806baf3545d3fd709aff56f1  results_src/m5/vendor/malarianet.py
594e3549e3cdc6fb93b0d50e54439c36a91b48919598c328c287fa79bbdba48d  results_src/mh/mh_infer.py
30575c1b7b3314c7a422da68fc75177f954d63c27c8a851fe2b5427a4898db53  stage1/p5stage/__init__.py
e638f61e146e69c8f6993f3eb9c2ea42d2a7f5edcd7bce4e2c2ba889e078de68  stage1/p5stage/prod.py
```
- 原方案甲的本机 CPU 推理不再使用；其设备检查（m5_device_check.json 39b111c8）未通过，记录保留。
- 本机文件（shasum 原始输出）：
```
049aade5ce6481fc5f7c7bf63e09d33ceb4e2cb190f5f050451509f98578f1f1  results_src/m5/vendor/a2_s41.pth
9d56a35ce40941d0c52467b0a3797301e39026e663b330af3c088c487e3e4014  results_src/m5/vendor/a2_s42.pth
b0125a2a70282adacc4fa76e2804d90ccba7ad66591a9298091d4ec7b62b538b  results_src/m5/vendor/a2_s43.pth
2da90395266ddc10c3e525387149a0fce4d5f7c3d83c754bc6fa2ac891c936ab  results_src/m5/vendor/a2_s44.pth
f4cb22539b8ee452a8279abdd1fc3e03e2bb72466af7105210fb067c700cb7b2  results_src/m5/vendor/a2_s45.pth
c51cc13deffd1a31cffd6af3dda1636dd5353f2f806baf3545d3fd709aff56f1  results_src/m5/vendor/malarianet.py
e638f61e146e69c8f6993f3eb9c2ea42d2a7f5edcd7bce4e2c2ba889e078de68  stage1/p5stage/prod.py
594e3549e3cdc6fb93b0d50e54439c36a91b48919598c328c287fa79bbdba48d  results_src/mh/mh_infer.py
```
- 前端与模型构建的代码逐字复制自 mh_infer.py（594e3549）第 11、13–24、28–33 行（results_src/m5/m5_infer_lib.py，启动时逐行断言）。

### 0.5 参数
- m = 22；(ε, δ) = (0.10, 0.10) 和 (0.05, 0.10)；5 种子 × 200 次抽样。
- **新的抽样流**（与开发站不同）：REV-1 抽样流 default_rng([20261022, key8, ai, seed, m, sidx]) 中，M5 HCM 的 key8 = **20261104**、sidx = **5**。开发站用过的 key8 / sidx 均不是这两个值（推理脚本断言）。SPLIT 分组流 [20261102, site index, seed, draw] 与 bootstrap 流 [20261101, site index, …] 中的 site index 取 **4**（开发站为 0–3：bbbc041_all、mpidb_full、lacuna_B、旧版 MP-IDB），同样不同。开发站的 key8 为 8 或 31，sidx 为 0、1 或 2。

## 1. 数据（实验管理只读核查与下载后体检，均已备案）
- 来源：Sultani 等，CVPR 2022（10.1109/cvpr52688.2022.02003）；许可按用户决定，引用原文，不发邮件。
- 只用 **HCM-1000x**。下载 1,257 个位置；其中 40 个为构造检查集（check_set.json），不进入评估。**评估位置 = 1,217 个**，m = 22 的抽样框就是这些位置。
- LCM-1000x **不参与确认，也不作探索性判读**（E-M5-1 第 5 节）：低价显微镜图像上检测器召回 0.39，标注对位误差中位数 0.72 个细胞半径，构造失败。这一事实写入论文讨论。
- 物种：*P. vivax*（论文数据集一节原文："we collected data from P.vivax infected patients"）。
- 视野 = 位置；没有病人或玻片编号。
- 标签：Ring 与 ring 合并；所有寄生阶段都算阳性。出界框先截断，面积为 0 或不到原框 50% 的丢弃（HCM 无此类框）。

## 2. 构造（已完成，结果见 0.3）
- 冻结的 Lacuna-B Hough 检测器（逐字复制，逐行断言），阴性与阳性在同一检测器上取景（Lacuna-B 重建规则）；r0 = 104.25。
- 构造门槛（HCM）：召回 0.790（≥ 0.70，过）；central_frac 中位差 −0.0040（≤ 0.05，过）；阳性 ≥ 300：PASS；目检 23 / 24 框住整个细胞（≥ 20，过）。
- 尺度门槛 r0 / R_cell = 0.443，经勘误 E-M5-1 判为方法无效，不调常数、不重跑。
- 全集（1,217 个位置）：阴性 180,344 个，每视野阴性 最小 4 / 中位 154 / 最大 226；阴性为 0 的视野 0 个，所以 F 与 REV-1 的抽样框相同。

## 3. 推理与评估（只运行一次；路径写死）
1. **读 manifest**：<LOCAL_PATH> 见 0.3），列 file、subsite、field_id、kind、idx、label、x0、y0、x1、y1、matched_iou；断言 subsite 全为 HCM、field_id 共 1,217 个且与 check_set.json 无交集、label 与 kind 一致、每个 file 存在。行序即输出行序。
2. **推理（VM，推理步，E-M5-3）**：5 个检查点按种子 41 → 45 的顺序，每个只推理一次（前端见 0.4）；每个种子写 logits 文件。推理步（m5_confirm_e3.py --stage infer）只写 ids_m5hcm.npz、5 个 logits npz、crop_sha256.txt、stream_record.json、SHA256SUMS，然后停；不加载 Platt、sp.json、NIH 验证集，不调用评估代码。
3. **NIH Platt 核对**：按 m3_k1v_precheck.load_nihval / platt 在 NIH 验证折上重算 Platt 参数，必须与 0.2 的数值逐位相等。
3b. **推理步与评估步分开（E-M5-3，两步分机）**：推理步在 VM 上按上述路线（必须是 PNG1，代码断言并把路线名写入 stream_record.json）得到 5 个种子的 logits 与 crop_sha256.txt，以 "x" 写入并写 SHA256SUMS。结果拉回 Mac（先到 scratchpad，核 SHA256SUMS，再放入 <LOCAL_PATH> 环境，即开发站所用环境；m5_confirm_e3.py --stage eval）先核：SHA256SUMS 恰好覆盖目录内文件且文件集合固定；crop_sha256.txt 行序等于 manifest 且逐行等于 Mac 清单；ids 等于 manifest；stream_record 的路线、总数与各 pin；logits 形状；再做第 3 步的 Platt（逐位相等，断言原样）与第 4 步的评估。任一核对失败，不创建输出目录。流式记录（每批的输入与输出 sha、服务端计数）写入 m5_confirm.json。**封存后、推理前，必须在 Mac 上先跑一次 `--stage eval --preflight`（pin、冻结代码、抽样键、manifest、Platt 含 NIH 验证折；不写任何东西）并通过，作为推理的门槛。**
4. **评估**：以导入方式调用 0.2 所列冻结代码（不复制、不修改；只新增 M5 的站点读取）：先按 REV-1 抽样流生成 5 × 200 次抽样，在同一批抽样上依次计算 S3 → c1s → SPLIT6；覆盖口径为全站 p（判定），另报 p_R。
   - **按设计丢弃的两个臂**：冻结循环（m3_split_dev_eval.run_site_seed）在同一批抽样上还会在内存中计算 K1 与 K1vB。它们的行在汇总之前就在 evaluate() 内丢弃（代码断言只有 S3、B1u_c1s、SPLIT6 离开该函数），不写盘、不打印、不进入任何报告或日志，**未读取**；它们不是 M5 的第二个或第三个读取对象。
   - p_R 所需的阈值取自 m3_splitr_dev_eval.py（8ca6b605）中复制的同一循环（SPLITS_R = ()）；两个循环在 S3、B1u_c1s、SPLIT6 上的行必须逐行相同（断言）。
5. **输出（E-M5-3，两个目录）**：推理步写 <VM_HOME>/P5_StainRef/runs/m5_infer_<UTC>/（ids_m5hcm.npz、a2_s{41..45}_m5hcm_logits.npz、crop_sha256.txt、stream_record.json、SHA256SUMS，每个文件以 "x" 打开），拉回 Mac 后放入 <LOCAL_PATH> <LOCAL_PATH> sha；SHA256SUMS；每个文件以 "x" 打开）。
6. **台账（E-M5-2、E-M5-3）**：先记一行旧运行「中止，未读结果，原因：传输」（Mac 上 m5_confirm_20261031T110535Z，只有 ids_m5hcm.npz，没有 logits；保留不删，标"中止"；推理 0 次、判读 0 次）。之后两行：推理步前写「开始推理（VM，第 1 次，也是唯一一次）」，完成后写「推理完成（VM）」；评估步前写「开始判读（Mac，第 1 次，也是唯一一次）」，完成后写「判读完成（Mac）」。读取计数：推理 1 次、判读 1 次。任何失败只能按勘误报告，不重跑、不改规则。
   - **输出目录创建之前的断言失败**（数据 pin、逐字段、冻结代码 sha、抽样键、manifest 检查、Platt 逐位比较，均在第 (3) 步之前、输出目录创建之前）：此时没有计算任何分数，**不算一次推理读取**。用勘误写明原因并修正后，可以用同一份代码、在新的输出目录中重新开始，并记入台账。
   - **输出目录创建之后的任何失败**：照上一条"只按勘误报告，不重跑、不改规则"处理。

## 4. 成功判据（单站：M5 HCM；η = 0 与全站 p 为判定口径；ε = 0.10 为主）
- **确认层**（ε 0.10）：
  1. SPLIT6 的覆盖率 ≥ 0.88；
  2. ΔSens_own(SPLIT6) − ΔSens_own(S3) ≥ +0.02；
  3. 没有标记（异常 > 1%）。
- **非劣层**：第 2 条未达 +0.02 但 ≥ −0.01，且第 1、3 条满足。写作"保证守住，灵敏度非劣，改善未确认"，稿中只能主张非劣。
- **未通过**：第 1 或第 3 条不满足，或第 2 条 < −0.01：写成"开发站上的改进，外部确认未通过"。
- **另报（不作判据）**：ε = 0.05 档的全部量；c1s 的覆盖率是否 ≥ 0.88；p_R 口径的全部量；类别条件覆盖；SPLIT6 选界比例与 (C2) 瓶颈次数；逐种子 AUROC。
- 无论结果如何，都同时写明本稿开头的已知劣势。

## 5. 风险栏（写于任何模型输出之前）
- 草案 3 第 5 节第 1–10 条全部保留。补充：
  11. **物种偏移**：M5 为 *P. vivax*，开发站以 *P. falciparum* 为主，与 IML 的情形相同。
  12. **检查集上阳性稀疏**：HCM 检查集每视野阳性占比中位 0.0075、均值 0.019、sd 0.033（无阳性视野 4 / 40）。草案 3 风险栏第 2 条担心的"阳性偏密、接近旧版 MP-IDB"在检查集上没有出现；这只是检查集上的描述，第 2 条照旧保留。
  13. **LCM 构造失败**（E-M5-1）：确认只剩一个子站，没有第二个站可以对照。
  14. **构造尺度门槛方法无效**（E-M5-1）：框是否框住整个细胞以目检为准（23 / 24）。
  15. **推理设备**：解码、前处理与推理都按第 0.4 节的 PNG1 路线在 VM 上进行（裁剪由 VM 本地只读副本提供），设置与存档管线相同。PNG1 路线上的设备复现检查（第 6 节第 1 条）必须通过。

  16. **开发站 logits 带有 TF32 误差（只报告）**：全部开发站与 NIH 验证折的 logits 都是用 TF32 默认设置在 GPU 上算的，约带 1e-3 的误差（设备诊断：关闭 TF32 后与 VM CPU 相差 ≤ 5e-6）。排序几乎不受影响（AUROC 差 ≤ 4.4e-5）。M5 用同一设置，所以可比。
  17. **评估环境（E-M5-3）**：评估在 Mac 的 deeplearning 环境里进行，与开发站评估相同。VM 环境（scipy 1.15.3，Mac 为 1.16.3）里冻结的 Platt 逐位断言失败（相差 ≤ 1e-7，跨环境核对 (a) 未通过），所以评估不在 VM 上做；M5 的 logits 在 VM 上算出，所以推理与评估之间多一道传输，由 SHA256SUMS 与 Mac 清单把住。核对 (b)（VM 上用冻结的 m3_split_dev_eval 重跑四个开发站，与 sp.json 364ebc08 逐格比，只作诊断）：56 格、852 个浮点字段、1029 个精确字段，精确字段 0 处不同、无判定翻转；超过 1e-9 的只有 6 个 Platt 参数（s41–s43，最大 9.2e-8）。即 scipy 版本差只落在 Platt 拟合上，没有传到格内判定；评估回 Mac 是为了保证 Platt 逐位相同。

## 6. 封存前必须完成
1. **设备复现检查（在 PNG1 路线上重做）**：m5_stream_devcheck.py 用同一份 1,000 个开发站裁剪（清单 ec006f11，Lacuna-B，前 500 个阴性加前 500 个阳性），把原始 png 字节经 PNG1 路线送到 VM，在 VM 上解码、前处理并在 GPU 上推理（服务端回显的文件 sha 必须等于 SHA256SUMS_vm），与存档的 CUDA 参照逐行比较。只读开发站数据，不碰 M5。
   - **判据不变**：5 个种子都满足「最大 |Δlogit| ≤ 1e-3，且 AUROC 差 ≤ 1e-4」（这条路线上预期为 0）；不满足就停下报告。必须通过之后才请用户封存。
2. 推理与判读脚本（m5_infer_lib.py、m5_stream.py、m5_stream_devcheck.py、m5_confirm_e3.py）、部署包 bundle3c / bundle4 与 smoke（m5_confirm_e3_smoke_local.py）经实验管理审查；smoke 全部通过（smoke 不读 M5 的任何文件）。脚本 shasum 原始输出（实验管理 2026-10-31 已部署核对）：
```
03ac6a67104b79719996af5b37c73c924d8df3a51444ffdb81fcc2b48290d8b5  results_src/m5/m5_confirm_e3.py
503a0e0d2667702198ce45602967b38dc2a273656d48572403d5856959d9af6f  results_src/m5/m5_confirm_e3_smoke_local.py
9a0863360c701157eaf6959e52537f1d310d1fbcd7bc93557716d9c50b79c535  results_src/m5/m5_stream.py
d3ac4d782568eceba4fe53455b200608ed4fd7579ae3a455882ea5ade72cfc90  results_src/m5/m5_infer_lib.py
```
   e3 smoke 输出 2db865ea10aec59044bc8d8dbe3d41507feffd94858fe8fe7f72fae576cecda7（(1)–(10) ALL PASS）。勘误：E-M5-2（裁剪拷到 VM、本地 serve、路径钉死）、E-M5-3（VM 只推理、评估回 Mac）。
3. 用户批准封存。封存后先在 Mac 上跑 preflight，通过后才运行推理步（VM，一次）与评估步（Mac，一次），即 m5_confirm_e3.py 的两个阶段各一次。
