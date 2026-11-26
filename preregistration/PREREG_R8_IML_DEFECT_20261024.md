# R8 IML 实现缺陷登记（2026-10-24）

母文件 PREREG_R8_IML_20261023.md（sha256 591d155de9162884208f93f25bbbae04758a2f5b908348d370befa8adf35a57f）及勘误（sha256 c29d33440eb614521f76cf31a773fa4dce5280d5c3f874825ece91cb79819320）。

**发现**：run_r8.sh（sha256 0c4f7221daa140447085d1eee5fe785e5d7a7d14719d956e24a17872fd1ab045，PID 150984，2026-10-24 07:47:03Z 启动）的第 3 步（推理）在加载模型时崩溃：脚本里写的 `P5_MODEL_PY=<VM_HOME>/P5_StainRef/code/malarianet_asd.py` 是我在上一段没有核实的猜测；该文件内部要读 `<VM_HOME>/MedPrune/models/malarianet.py`，那个路径在虚拟机上不存在（FileNotFoundError），进程在加载任何检查点之前退出，07:48:09Z 之前结束，watchdog 随之正常退出。正确值取自已成功运行过的 iter1_20261021/code/run_s3.sh：`P5_MODEL_PY=<DATA>/medprune/MedPrune/models/malarianet.py`。

**状态**：第 0 步（重叠核对：无命中）、第 1 步（裁剪：38,428 张）、第 2 步（打分前报告：与预注册计数无不符）已正常完成，输出保留在 <DATA>/p5_e8_candidate/r8_iml/crops 与 ~/P5_StainRef/runs/r8_iml/。**没有产生、也没有查看任何模型得分、logit 或覆盖率数字**；infer/、rev1/、rev7/ 目录都不存在。

**处理（提案，待 P4 与实验管理确认后才执行）**：新脚本 run_r8b.sh 从第 3 步续跑，只把 P5_MODEL_PY 改为 <DATA>/medprune/MedPrune/models/malarianet.py；开头重新校验脚本、annotations.json、裁剪计数（r8_check_counts.py）；第 3–5 步命令、参数、资源限制、watchdog 不变；第一次运行的日志保持原名，新日志与耗时文件带 _b 后缀。S3、规则、随机数、判据、预注册均未动。
