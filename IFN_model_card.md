# Frozen IFN model card (2026-09-12)

## 目的
在接触任何 S:CORT 外部验证数据之前, 固化"IFN 分数"的完整定义与训练数据口径,
使 external validation 成为真正的盲测 (pre-registered definition)。

## 冻结的模型定义

**特征**: IFN score = ssGSEA(HALLMARK_INTERFERON_ALPHA_RESPONSE) 与
ssGSEA(HALLMARK_INTERFERON_GAMMA_RESPONSE) 两分数的算术平均,
实现 = GSVA::gsva(ssgseaParam(expr_symbol_matrix, ifn_sets)), 表达矩阵为
队列内 symbol 化(首 token 别名 + 最大方差去重)后的**队列内原始刻度**
(不做跨队列标准化; ssGSEA 为样本内秩次方法, 对单调变换不变)。

**方向约定**: 正 g = 缓解组更高。验证时的检验假设: external cohort 中
responders 的 IFN 分数高于 non-responders (单侧方向检验 + 双侧报告)。

**无阈值**: 本模型不定义 cut-off; 验证指标 = 方向一致性 (g>0)、标准化均数差、
AUC(其自然尺度), 不做切点优化(遵循导师意见#2)。

**训练数据**: 9 个 nCRT 队列 (features_GSE*.csv, 见哈希表; 不含 GSE213331
化疗队列、GSE56699 放疗队列、GSE46862 原代培养、任何 S:CORT 数据)。

**冻结的合并估计**: IFN g=0.384 (95% CI 0.161–0.606), P=0.0041, I²=1%
(REML + modified Hartung-Knapp, meta_pooled_results.csv 行 primary_ncrt9/IFN)。

## 输入文件 SHA256 (model_inputs_sha256.txt)

| 文件 | 角色 |
|---|---|
| h.all.symbols.gmt (8ff1f036…) | 基因集定义 (MSigDB Hallmark) |
| mcpcounter_genes.txt (408f6c5d…) | MCP 标志基因集 (NK/Cyto 支持特征) |
| project_config.R (01b6fa37…) | 队列注册表与分析集定义 |
| cohort_utils.R (d17cf36b…) | symbolize/build_cohort_features 实现 |
| meta_utils.R (4dc8f8a2…) | hedges_g/pool_random_effects 实现 |
| meta_pooled_results.csv (fe7167b7…) | 冻结的合并估计 |
| features_GSE<9 nCRT>.csv | 9 队列样本级分数与标签 |

## External validation 预案 (⑤, 数据到手即执行, 不再改动)

对每个 S:CORT 队列 (Grampian / ARISTOTLE; Copernicus/TREC 若标签可得):
1. 表达矩阵 → symbol 化(与 cohort_utils.R 同一函数) → 同一 ifn_sets ssGSEA → 取均值
2. 按队列内存活/疗效标签二分 (定义取原文, 不重编码)
3. 报告: Hedges g (responders vs poor) + 95%CI、Wilcoxon P、AUC (poor/good 方向固定)
4. 汇总: 不与训练集重合并; 仅做"9 队列估计 g=0.384 是否落于外部 g 的 95% CI 兼容范围"与
   方向一致性计数
5. 通过标准 (预先声明): 方向一致 + P<0.05 (至少一个队列) 或随机效应合并后 P<0.05;
   失败则如实报告, 不得改定义、换基因集或找切点

## 变更纪律
冻结后任何对基因集/打分/分组/方向的改动都必须在本文件记录日期与理由,
并使验证结果标记为 post-hoc。
