# 勘误：PREREG_P5_MPIDB_HOSTCELL_REEVAL_20261025.md（封存 sha256 ecd3d6ac…）

写于 2026-10-25，在任何脚本对新集运行之前。原预注册不改，本勘误只澄清一处自相矛盾的数据边界表述。

## 矛盾在哪

- 第 3 节规定：S3-G 照 `r8_rev7.py`，Table 3 的无标签校正照 `rev5.py`，"定义一律照各自原脚本"，来源阈值"读取已存档的值"。
- 但 Q6 里写的是"断言不读任何 NIH 文件"，第 5 节写的是"不读 NIH 的任何图像（包括测试折）"。
- 照原定义：
  - S3-G 的固定网格取自**源验证折分数**的分位数；
  - UPA 用到源分布和源流行率，这些都来自已存档的 **NIH 验证折分数与标签**（`INF/ids_iter1.npz` 中验证折的标签、`INF/<stem>_scores.npz` 中的 `nih_validation_logits`、`INF/<stem>_eval.json` 中的来源阈值）。
- 所以"不读任何 NIH 文件"与第 3 节的定义相矛盾，是主会话起草时表述不严谨。脚本改编（fork，mh_rev7.py 和 mh_rev5.py）照第 3 节的原定义执行，这个矛盾就暴露了出来。

## 澄清后的边界（以此为准）

- **允许读：** 已存档的 NIH **验证折**分数、标签和来源阈值（只读数组，不是图像）。它们是 S3-G 和 UPA 定义所需的输入，与 PlasmoCount、IML 确认评估时使用的完全相同。
- **不读：** 任何 NIH 图像；NIH **锁定测试折**的任何数组（`nih_locked_test_*` 键），脚本里只按键名取 `nih_validation_*`；任何确认站点；BBBC041 和 Lacuna 的数据。
- 判据、方法和其余边界都不变。

## 实现要求

- mh_rev5.py 和 mh_rev7.py 读取 npz 时只访问 `nih_validation_*` 键，并断言没有访问任何以 `nih_locked_test` 开头的键（在加载处检查键名）。
- 包装脚本的日志写明本次读了哪些 INF 文件。
