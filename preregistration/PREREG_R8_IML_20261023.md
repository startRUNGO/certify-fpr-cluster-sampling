# P5 R8 预注册：IML Malaria 第二确认站点（一次性确认评估）

写于 2026-10-24（CEST），早于对 IML 图像的任何裁剪、推理、模型得分。设计：P4 会话；执行口径补全：P5 会话（标“执行口径”）。模板：PREREG_E8_PLASMOCOUNT_20261023.md（sha256 2ce88dfde044d508e8a4ae258086e25b67c343d4e6f9df8fbf914d84f05aabda）。

用户批准（本会话内亲自批准下载，经 AskUserQuestion 回答“批准下载”）：仅在 Mac 上下载 GitHub 仓库 QaziAmmar/A-Dataset-and-Benchmark-for-Malaria-Life-Cycle-Classification-in-Thin-Blood-Smear-Images 的 IML_Malaria/ 文件夹与根目录 annotations.json，存到 Paper5_StainRef/data_r8/iml/；不进 Git、不进投稿包、不上传；论文只报告汇总结果，不再分发图像。用户经 P4 转述的原话：“批准下载 IML Malaria 做 R8：先写预注册，再向实验管理申请虚拟机时段，批准后再跑；论文里只报告汇总结果，不再分发图像。”虚拟机上的任何运行须先由 P4 向实验管理会话申请，获“可以开始”后才开始。

## 0. 冻结物（不得改动）
- PROTOCOL_FREEZE_20261022.md（sha256 f2565f309fff346cd8056ad2e543bea9bc3031041da7ce7be54b5dbb5aa0802b）：规则 S3 及评估设置；
- s11.py（sha256 c7d3ad62…fca6f，母文件 PREREG_ITER11_S11）；
- 裁剪 PAD=5 与缩放到 128×128 的推理流程（e7_infer.py，sha256 9c666a9f…2f5；E8 的 e8_infer.py 为其同流程副本，sha256 18b25467…78d0e 的数据加载路径版）；
- 35 个检查点（7 臂 × 5 种子，~/P5_StainRef/runs/main/ckpt/*.pth），A2 为主、A0-ctrl 为次；
- IML 没有参与 S3 的任何开发（11 轮迭代只用了 BBBC-ALL、MP-IDB、Lacuna-B），也没有参与 E8；本预注册之前没有看过任何 IML 的模型得分。

## 1. 数据（已在 Mac 上读取元数据，未打分、未推理）
- 来源：见上；论文 Arshad et al., Neural Computing and Applications 34:4473–4485 (2022), doi 10.1007/s00521-021-06602-6。**许可证：仓库无 LICENSE，数据集许可证未查到**；因此只在内部评估、论文只报汇总数字。
- 下载：IML_Malaria/ 共 346 个文件（345 张 JPG + mean_image.png，后者不用）与根目录 annotations.json；git blob sha1 与仓库树逐个核对，0 个不符。逐文件 sha256 见 data_r8/iml/MANIFEST_sha256.txt（347 行）。
- **annotations.json：sha256 c9b180b43fb65c18fa62fb4f752dc089dca4fb23659a2073e5eac2c4fd2a5db8，3,242,388 字节。**
- 标注读取结果（JSON，345 条记录，全部有对应图像）：图像 1280×960（sips 读取 345/345），RGB；框总数 38,449。原始类名及框数：red blood cell 37,899；gametocyte 261；ring 164；trophozoite 77；schizont 27；difficult 21。每图 red blood cell 数最小 14、中位 110、最大 164；含感染细胞（四个期别任一）的图 318 张，感染细胞合计 529（ring 164 + trophozoite 77 + schizont 27 + gametocyte 261）。
- **标签：阳性 = gametocyte、ring、trophozoite、schizont；阴性 = red blood cell（已标注，不用检测器）；difficult（21 个框）排除，既不当阳性也不当阴性、不推理、不进入中心化统计。**
- 视野 = 图像（无玻片或患者 ID；这是局限）。
- 门槛（在 JSON 上即可判定，已满足）：(a) 每图 red blood cell 数的中位数 110 ≥ 5；(b) 至少一张图含感染细胞（318 张）。
- 框几何：bbox 为 {x,y,h,w}（字符串），x=列、y=行、w=宽、h=高；全部 x,y≥0；1,250 个 red blood cell 框的右边或下边越过图像边界，且全部恰好超出 5 像素（其他类别无越界）。**执行口径：裁剪时按冻结流程 [y−5, y+h+5]×[x−5, x+w+5] 并截断于图像边界；不另行修正。**（此观察的含义，即标注框是否已含 5 像素外扩，留给 P4 判断；不改变本预注册。）
- 标本：P. vivax、Giemsa、XSZ-107 显微镜 100×、拉合尔当地医院（来源论文）；与 BBBC041、MP-IDB、PlasmoCount、Lacuna、NIH 的重叠：来源未声明；Mac 上文件名与 md5 核对结果见 data_r8/iml/OVERLAP_CHECK.md（由 P5 在完成后写入；Lacuna 与 NIH 待在虚拟机上核对，本次不上虚拟机去查）。

## 2. 打分前报告（读图裁剪之后、推理之前）
成功读取的图像数、图像尺寸/模式分布、各类框数、裁剪后阳性与阴性数、每图阴性数最小/中位/最大、退化框数；与上面 JSON 计数不符如实报告。若任何图像为非 8 位或不能转 RGB，在看得分之前停止并报告。

## 3. 推理与规则
- 裁剪：新脚本 r8_crop.py（改自 e8_crop.py，只改 JSON 读取、路径和类名映射）；输出 crops/ 与 manifest.csv（含 field_id=图像名）。
- 推理：r8_infer.py（改自 e8_infer.py，只改 crops 路径和输出名 ids_iml.npz / *_iml_logits.npz），35 个检查点，批量 512，保存逐样本 logit。
- 规则：一次性运行冻结规则 S3（m=22；(ε,δ)=(0.10,0.10) 与 (0.05,0.10)；A2 为主、A0-ctrl 为次；5 种子 × 200 次抽样；随机数 numpy default_rng([20261022, 14, 臂序号, 种子, 22, 0])；候选 c0/c1s/c2s，δ/3 的 B1u，固定序列，选择规则同冻结文件）。S3 对每个臂和设定只运行一次；看到数字后不改 S3。
- 分析脚本：r8_rev1.py = rev1.py（sha256 7bcb1635464757458968be79c2a33cd4bd2b1c77e7a4d31bd9fa870d49fc3ae0）加站点 “iml”（key8=14，m 只取 22），只改数据加载；r8_rev7.py = rev7.py（sha256 e55829156204597b4c069ae0a2884b1e989cbad95862acae6f561419a6167499）同样只加站点，用于 S3-G。评估代码只允许改数据加载（JSON→裁剪、路径、站点名），差异存档并记录 sha256。
- 报告（均为 REV-1 的修订定义）：覆盖率 P(FPR ≤ ε)（200 次抽样 × 5 种子）；ΔSens_ε = 灵敏度 − oracle_ε 灵敏度（oracle = 全站点阴性在未中心化分数上 FPR=ε 处的灵敏度）；ΔSens^own（oracle 在规则自己的中心化分数上）；cannot-certify 率；实现 FPR 的均值、中位数、90% 分位；点击数（= 所抽图内感染数）与已标注阴性数。η：阴性是人工标注，主口径 η=0，并列 η=0.10 下界 max(0, FPR−0.10)。
- **灵敏度的不确定性**：感染细胞仅 529 个，所以所有灵敏度和 ΔSens 一律附 90% 区间，分别给出种子层面（5 个种子的区间）与抽样层面（1,000 次抽样的 5%–95% 分位）。**主判据落在覆盖率上**；灵敏度条件在区间下如实说明。
- 同批配对抽样同时报告：R1-2（每图 2 个阴性）、R1-field、R1-all、B1、B1u，(0.05,0.10) 一并报告。
- 次要判据：S3-G（源阈值固定网格，母文件 PREREG_REV_FIXEDGRID_20261023.md，sha256 06f6b5b5…），同样 m=22、(0.10,0.10)。
- REV-4 未触发，因此没有 R4 修正候选；本评估不含 R4 的任何变体。

## 4. 成功判据（A2，(ε,δ)=(0.10,0.10)，m=22；η=0 为判定口径，η=0.10 并列）
- 覆盖率 ≥ 0.88；
- 且 ΔSens_0.10(S3) ≥ −0.05，或 ΔSens_0.10(S3) ≥ ½·ΔSens_0.10(R1-2)（代价不超过 R1-2 代价的一半）。
- 两个条件都满足为“通过”，否则“未通过”，如实报告。
- 以上判据换用修订后的参照 oracle_ε（PlasmoCount 的 E8 判据用的是旧参照），两者不可混用；PlasmoCount 的结果不受本预注册影响。

## 5. 资源（供实验管理审批）
- 机器：A6000 虚拟机 <VM_HOST>，python 环境 <VM_HOME>/miniconda3/envs/pyproject/bin/python；不碰 ~/P4_rerun、~/P4_samemod、~/P4_nih。
- GPU 推理：35 个检查点，批量 512，同一时间最多 1 个作业；显存 ≤ 30 GB（全部 38,449−21 = 38,428 个裁剪片，每个检查点显存很小，估计峰值 ≪ 30 GB）；nice -n 19 ionice -c3；`ulimit -v` 内存上限按实验管理当前规定（上次实验：`ulimit -v 3000000` 与 OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=1；如规定已变，以规定为准）。
- CPU 分析：并行进程 ≤ 4，同样 nice -n 19 ionice -c3，同样内存上限。
- 上传：**预注册默认：向 <DATA>/p5_e8_candidate/r8_iml/ 上传 345 张原图（约 66 MB）和 annotations.json，在服务器上裁剪**（裁剪片约 38,428 张小 PNG，估计几百 MB 以内）。图像不进 Git、不进投稿包。
- 估计用时（估计，不是测量）：裁剪 < 10 分钟；推理 35 个检查点合计 < 30 分钟；CPU 分析（rev1 风格，2 臂 × 5 种子）数小时以内；总计 < 4 小时。
- 输出：~/P5_StainRef/runs/r8_iml/，不覆盖任何已有结果。
- 能否中途停止：可以。推理逐检查点保存，分析逐 (臂, 种子) 任务保存；中途停止后可从已保存部分继续，不影响已有结果。
- 起止 UTC 与峰值 RSS 在运行后报告。

## 6. 时间顺序与纪律
本文件的 sha256 与创建时间必须早于任何 IML 的裁剪输出、推理输出与得分；sha256 写入同目录 PREREG_R8_IML_20261023.md.sha256。脚本 r8_crop.py、r8_infer.py、r8_rev1.py、r8_rev7.py 在本预注册之后写成并记录 sha256（results_src/r8/r8_scripts.sha256），仍早于任何推理。异常先停下报告；发现实现缺陷先登记缺陷文件再改。运行后向 P4 会话报告。
