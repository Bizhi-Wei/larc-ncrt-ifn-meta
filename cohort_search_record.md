# 队列系统检索与筛选记录 (PRISMA-style)

**检索日期**: 2026-09-01
**数据库**: NCBI GEO DataSets (gds)
**检索式**: `(rectal[Title]) AND (neoadjuvant OR chemoradiotherapy OR chemoradiation OR radiochemotherapy) AND ("expression profiling by array"[Filter] OR "expression profiling by high throughput sequencing"[Filter]) AND "Homo sapiens"[Organism]`
**命中**: 38 个数据集；另经引文追踪补充 2 个 (GSE109057, GSE94104)。2026-09-10 的检索后引文/基准队列复核再发现 1 个原检索式漏网队列 (GSE56699，其标题不以 rectal 开头) → 共评估 41 个

## 纳入标准与分析范围

1. 人类直肠癌 (rectal) 肿瘤组织
2. 新辅助治疗前活检/手术标本；nCRT 队列进入主分析，单纯放疗或单纯新辅助化疗队列仅进入各自方案敏感性分析
3. bulk 转录组 (芯片或 RNA-seq)
4. 可公开获取的逐样本疗效标签 (pCR / TRG / response 二分)
5. 探针→基因注释可程序化获得

## 纳入 (9 个 nCRT 队列 577 例；另有放疗 51 例作方案敏感性分析)

| 队列 | n | 平台 | 疗效定义 | 备注 |
|---|---|---|---|---|
| GSE209746 | 93 | RNA-seq | 2年持续完全缓解 CR/iCR | 发现队列 (PMID 35970919) |
| GSE35452 | 46 | Affymetrix GPL570 | JSCCR 退缩分级二分 | UFT/LV 方案 (PMID 24316942) |
| GSE87211 | 202 | Agilent GPL13497 | pCR=ypT0N0 | 203 个肿瘤中 1 例 ypT/ypN 均缺失，排除于疗效分析；DFS/OS (PMID 29119627) |
| GSE150082 | 39 | Agilent GPL13497 双色 | TRG 二分 | PMID 32784964 |
| GSE45404 | 42 | Affymetrix GPL570 | Mandard TRG 二分 | PMID 26023803 |
| GSE119409 | 56 | Affymetrix GPL570 | sensitive/resistant (capecitabine+50.6Gy) | PMID 33106387 |
| GSE133057 | 33 | Illumina GPL6102 | AJCC/CAP 评分 0-1 vs 2-3 | PMID 31704889; 注释覆盖 51% |
| GSE53781 | 26 | CodeLink GPL18134 | TRG1-2 vs 3-5 (Mandard) | PMID 25380052; 注释覆盖 38% |
| GSE94104 | 40 | Illumina HT-12 (FFPE) | TRG3 vs TRG1-2 (升序制, 依 PMID 28621318 校正) | 仅治疗前活检亚组 |
| GSE56699 | 51 | Illumina HumanWG-DASL v4 (FFPE) | Mandard TRG1-2 vs TRG3-5 | 术前放疗、非 nCRT；58 张治疗前芯片来自 52 名患者，去除 1 名缺失 Mandard 分级者后 51 人（23/28）；重复活检先按患者平均后评分；仅作放疗方案敏感性；PMID 25706627 |

## 排除 (31 个) 及原因

### 疗效标签未公开寄存 (4)
| 队列 | n | 说明 |
|---|---|---|
| GSE216616 | 298 | JAMA Netw Open 2023 (PMID 36662520; JFCR 日本); GEO 仅有 tissue/age/sex, TRG/W&W 标签未寄存; 该研究独立报道了与本研究一致的结论 (细胞毒性淋巴细胞 MCP-counter 分数预测疗效, 并在 GSE87211/GSE45404 验证), 作为趋同证据引用 |
| GSE109057 | 91 | 同上 (Cancer Institute Hospital; 标签未寄存) |
| GSE318268 | 73 | Galunisertib+CRT 试验 (PMID 41653703); GEO 无疗效字段, 试验补充材料无法程序化提取 |
| GSE242786 | 15 | 计数表列名 (sample60-80) 与 GEO 样本标签无法映射 |

### 平台/注释不可行 (2)
| GSE123390 | 33 | GPL17586 (HTA 2.0) 平台注释在分析时点无法获取 (NCBI 服务异常), 仅有原始 CEL |
| GSE139255 | 156 | 定制低密度面板 (784 基因), 转录组覆盖不足, 无法支持 Hallmark/ssGSEA 打分 |

### 重复/重叠 (1)
| GSE40492 | 245 | 德国 GRCSG 试验队列, 与 GSE87211 源自同一 ~300 例采集 (治疗字段完全相同), 判为近重复, 保留注释更全的 GSE87211 |

### 设计不适用 (24)
- 原代癌细胞培养而非直接 bulk 肿瘤组织: GSE46862 ("primary cancer cell")
- 治疗后切除组织、无可用治疗前疗效标签: GSE213331 (2026-09-12 撤出; 带标签 18 例均为 nCRT 后切除标本)
- 显微切割且治疗前肿瘤上皮样本过少: GSE93375 (~10)
- 实验药物联合 (贝伐珠单抗) 且治疗前肿瘤 ~20: GSE60331
- 10 例配对前后设计: GSE15781
- 单细胞测序: GSE190826, GSE254249, GSE278405, GSE278406
- 血浆转录组: GSE333629
- 外周血单核细胞: GSE44172
- 临床前药物筛选/细胞系: GSE205035, GSE207521, GSE304710, GSE294953, GSE316027, GSE319208
- CAF 临床前模型为主: GSE156281
- 治疗前后配对免疫标志 (22 例, 设计不符): GSE233517
- 结肠为主/单核细胞: GSE221926, GSE221925, GSE222300
- 甲基化 (非表达): GSE145037
- 无新辅助注释 (误用记录): GSE39582 (源自 Marisa 结肠癌分型队列, PMID 23700391; PMC9558072 所称 "200 例 nCRT LARC" 无法核实)

## 流程计数 (供 PRISMA 图)

- 鉴定: 38 (预设检索) + 2 (初始引文追踪) + 1 (2026-09-10 检索后复核发现 GSE56699) = 41
- 排除: 31 (标签未寄存 4 / 平台注释 2 / 重复 1 / 设计不适用 24，含 GSE213331)
- 纳入定量综合: 9 个 nCRT 队列、577 例疗效可分类肿瘤 (93+46+202+39+42+56+33+26+40)；另有 GSE56699 放疗队列 51 名独立患者，仅进入放疗方案敏感性分析。合计 10 个数据集 628 人。TCGA-READ 164 例仅作预后旁证。

> 注: GSE209746 93 例为 94 例中去除 1 例缺 RNA (P-0044433)；GSE87211 的 203 个肿瘤样本中，1 例因 ypT/ypN 均缺失不进入疗效分析；GSE94104 为 80 例中的 40 例治疗前活检；GSE56699 的 58 张治疗前芯片不等于 58 名独立患者（52 名患者，其中 51 名有 Mandard 分级），分析在患者层面完成。非 nCRT 队列不进入主分析。

## 更正 (2026-09-12): GSE213331 撤出全部数据集

GEO 官方记录 (!Series_overall_design / !Sample_description) 证实: 27 样本 = 9 pCR(切除标本, nCRT后) + 9 npCR(切除标本, nCRT后) + 9 Pre(治疗前活检, **无疗效标签**)。此前将其作为"新辅助化疗疗效队列"(n=18)纳入敏感性分析是**错误**: 带疗效标签的 18 例为治疗后组织, 不能用于治疗前生物标志物关联; 且其治疗实为 nCRT 而非单纯化疗(两次误判)。

- 处置: 从全部分析集移除; features 表改名 REMOVED_features_GSE213331_posttreatment.csv 留档; project_config 版本号 2026-09-12-gse213331-removed-v1
- 影响: 主分析 9 队列 n=577 完全不变; "化疗敏感性(10 数据集)"与"全 11 数据集"分析组删除; 样本总量 646→628
- 教训: 疗效标签必须与样本时间点(treatment-naive)同时核验, 不能只看标签字段
