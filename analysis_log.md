# 多队列转录组分析日志

**课题**: 基于多队列转录组分析的局部晚期直肠癌新辅助治疗反应相关分子特征研究
**分析日期**: 2026-08-31 ~ 2026-09-01
**环境**: R 4.6.1 (Windows), Bioconductor 3.23, GEOquery / limma / DESeq2 / GSVA / pROC / survival

---

## 一、数据集

| 队列 | 平台 | 总样本 | 分析样本 | 表型定义 | 分组 |
|---|---|---|---|---|---|
| GSE209746(发现队列) | RNA-seq (STAR→GRCh38 原始 count) | 114 | **93** | 论文 Table S1A：两年维持完全缓解(CR)/未完全缓解(iCR) | CR=26, iCR=67 |
| GSE35452 | Affymetrix GPL570 | 46 | 46 | TRG 类 Responder/Non-responder | R=24, NR=22 |
| GSE87211 | Agilent GPL13497 | 363(含160黏膜) | 203(肿瘤) | 由 ypT/ypN 推导：ypT0N0=pCR | pCR=35, non-pCR=168 |
| GSE150082 | Agilent GPL13497(双色) | 39 | 39 | TRG 归纳 Good/Poor | Good=16, Poor=23 |
| GSE45404 | Affymetrix GPL570 | 42 | 42 | Mandard TRG 归纳(原文 PMID 26023803) | R=19, NR=23 |
| TCGA-READ(旁证) | RNA-seq (Xena star counts) | 177 | 164(原发、排除1例新辅助) | 非新辅助队列, 仅用于预后旁证 | OS事件25, PFI事件36 |

### 关键决策记录

1. **GSE209746 的疗效标签没有用 GEO 自带的 `response:ch1` 列**(27 CR/70 iCR),而是用论文 Nature Medicine (PMID 35970919) 补充表 Table S1A 的官方定义(`Use_in_outcome_analysis=1 & Has_RNASeq=1 & Use_for_Response_2year=1` → CR=26/iCR=68)。两者有差异：有 1 例 pCR 后两年内复发被论文改判为 iCR。论文 94 例中 1 例(P-0044433)的 RNA 数据未上传 GEO,实际可分析 93 例。
2. **GSE87211 没有直接的 pCR 列**,pCR 定义为 `depth of invasion after rct=0 且 lymph node metastasis after rct=0`(ypT0N0)。
3. **GSE122288 误下载后排除**(实际是脐带血甲基化数据，与主题无关)。
4. GEOquery 2.80 解析空表达表的 series matrix(RNA-seq 数据集)有 bug,绕过 `getGEO()` 直接调用 `GEOquery:::parseGSEMatrix()`。
5. **GSE39582 不可用**:文献(PMC9558072)声称从中选取 200 例 nCRT LARC,但 GEO 全部注释(含 family XML)无任何新辅助/放疗字段,tumor.location 仅 distal/proximal(结肠亚部位,源自 Marisa 结肠癌分型队列)——该说法无法核实,弃用。改用 TCGA-READ 做预后旁证。
6. GSE45404 的基因注释直接复用 GSE35452 的 fData(同平台 GPL570),覆盖率 100%。

---

## 二、分析流程与结果

### 1. 发现队列复现(GSE209746, DESeq2)

- count 矩阵 19662 基因 × 93 例，全部非负整数;STAR 原始计数直接可用
- 宽松低表达过滤(总 count≥10)后 18798 基因;DESeq2(Wald, BH 校正),iCR vs CR
- **复现论文 15 个蛋白编码差异基因中的 14 个**(|log2FC|>1 且 FDR<0.05),方向 15/15 一致
- IGF2: log2FC=+3.83, FDR=2.3e-08(iCR 高表达);L1CAM: log2FC=+2.36, FDR=4.9e-04
- 唯一未复现:CPS1(方向一致 -1.60 vs 论文 -3.05,但 padj=0.24)
- 注意:严过滤(≥26例 count≥10)会误删 VAX1/NHLH2/SOX2(少数 iCR 样本特异高表达),复现已发表论文时应用宽松过滤
- **校正 MSI 后**(design=~msi+response)IFN 通路及 IGF2/L1CAM 仍显著:IFN-α P 2.9e-3→2.3e-2;IGF2 padj=4.3e-07;L1CAM padj=7.8e-04

### 2. 单基因跨队列验证(15 基因方向一致性)

| 队列 | 一致/总数 | 备注 |
|---|---|---|
| GSE35452 | 8/14 | 仅 INHBB 方向一致且 P<0.05;IGF2 方向一致不显著,L1CAM 方向相反 |
| GSE87211 | 6/14 | 无显著;INHBB 再次一致且 P<0.05 |
| GSE150082 | **10/14** | IGF2 方向一致 P=0.063(趋势) |

结论:15 基因特征(尤其 IGF2/L1CAM)是严格"完全缓解"表型特异的,在 TRG 宽松表型队列中稀释;**INHBB 是唯一在 3 个芯片队列都方向一致的基因**。

### 3. 通路层面(Hallmark 50, limma camera,方向统一为"好反应组 vs 差反应组")

- 17/50 通路四队列方向一致(前三个队列时;第四队列加入后核心结论不变)
- **最稳定主题:干扰素/免疫激活**
  - IFN-α: 四队列全部 Up;P = 7.2e-03 / 0.065 / 7.0e-13 / **9.0e-18**
  - IFN-γ: 四队列全部 Up;P = 4.3e-03 / 0.145 / 4.9e-11 / **1.6e-15**
  - IL6-JAK-STAT3、炎性反应同向;代谢(糖酵解/OXPHOS)在好反应组方向一致下调
- **排除的主题**:缺氧(方向矛盾)、增殖 E2F/MYC/G2M(GSE35452 反向)、EMT/间质(GSE87211 反向显著)、DNA 修复(不一致)

### 4. 免疫浸润(MCP-counter ssGSEA,10 类细胞)

- **NK 细胞、细胞毒性淋巴细胞四队列方向全部一致(好反应组高)**;GSE87211 显著(P=0.009/0.016),GSE150082 显著(NK P=0.021,另 B lineage P=0.025)
- 方法学验证:我们的 ssGSEA 与论文 BINDEA 分数高度相关(T cells r=0.96, B lineage r=0.96, Neutrophils r=0.89, Cytotoxic r=0.87)
- 发现队列中单一细胞类型 CR vs iCR 均无 FDR 显著,但**免疫热聚类 IG3/IG4 的 CR 率 64% vs IG1/IG2 的 23%(Fisher P=0.043)**——IFN 信号来自免疫热亚群
- 注意:MCP-counter 的 CD8 签名仅 1 个基因,CD8/NK 绝对分数需谨慎,方向结论不受影响

### 5. 预测效能

| 模型 | GSE209746 | GSE35452 | GSE87211 | GSE150082 |
|---|---|---|---|---|
| IFN 评分单独 | 0.556 | 0.578 | 0.599 | **0.764** |
| IFN+IGF2+L1CAM(系数固定迁移) | 0.638(训练) | 0.500 | 0.618 | 0.761 |

结论:IFN 评分在表型定义清晰的队列中有中等偏上预测力,但跨表型迁移不可靠;**定位为"关联特征"而非独立预测标志物**。

### 6. 生存分析(GSE87211,DFS/OS)

- DFS(175例,34事件):IFN 连续评分 HR=0.89(0.64–1.25),P=0.51;KM 中位切分 P=0.33
- OS(190例):HR=0.79(0.50–1.25),P=0.32
- 方向均为保护性,事件数不足导致不显著

### 7. 四队列 meta 分析(IFN 评分,Hedges g)

| 队列 | g (95%CI) |
|---|---|
| GSE209746 | 0.28 (-0.18, 0.73) |
| GSE35452 | 0.24 (-0.33, 0.81) |
| GSE87211 | 0.48 (0.11, 0.84) |
| GSE150082 | 0.96 (0.30, 1.62) |
| **随机效应合并** | **0.44 (0.19, 0.70), P=6.9e-04** |

异质性极低:I²=11.6%,Q P=0.335。**效应方向跨平台、跨表型定义完全一致**。

---

## 三、总结论

1. **干扰素/免疫激活是 LARC 新辅助治疗反应最稳健的跨队列分子主题**(4 队列、3 平台、3 种表型定义方向全部一致,meta P=6.9e-04,I²=12%)
2. 机制链:治疗前 IFN 信号高 ← NK/细胞毒性淋巴细胞(B 细胞)浸润高 ← 免疫热亚群 → 完全缓解率约 3 倍(64% vs 23%)
3. IGF2/L1CAM 为严格完全缓解表型特异的耐药相关基因,校正 MSI 后仍成立
4. 建议标题收窄方向:干扰素/免疫激活主题

*(2026-09-01 深化分析后的修订版结论见第六节约 7 小节)*

## 四、主要输出文件

| 文件 | 内容 |
|---|---|
| **Fig1_main.png**(3800×4400, 300dpi) | **主图终版:A=IFN九队列 nCRT 森林图,B=NK/Cytotoxic九队列森林图,C=IG分组CR率,D=pCR-KM** |
| fig_A_forest_IFN_all.png / fig_B_forest_immune_all.png | 主图 A/B 单面板（九队列 nCRT） |
| fig_C_ig_crrate.png / fig_D_km_pcr.png | 主图 C/D 单面板 |
| ~~IFN_forest_4cohorts.png / meta_forest_5cohorts.png / immune_meta_forest.png~~ | 旧版 base-R 森林图(文字重叠,已被 fig_A/fig_B 取代) |
| hallmark_3cohorts_barplot.png / immune_cells_3cohorts.png | 通路/细胞跨队列对照图 |
| IFN_score_ROC_3cohorts.png / combined_model_ROC.png | ROC 曲线 |
| GSE87211_IFN_DFS_KM.png | 生存曲线 |
| GSE209746_DESeq2_iCR_vs_CR_loose_filter.csv | 发现队列完整 DESeq2 结果 |
| hallmark_3cohorts_camera.csv / GSE150082_hallmark_camera.csv | 通路富集结果 |
| immune_ssgsea_3cohorts.csv | 免疫细胞评分对照 |
| IFN_meta_4cohorts.csv | meta 分析效应量 |
| paper_94_labels.csv / paper_S8A_results.csv / paper_S1D_bindea_scores.csv | 论文官方标签/结果/免疫评分 |

## 五、遗留与后续

- [x] ~~GSE35452 的 Responder 定义未逐例核对原文~~ → 已核对,见六.5
- [ ] GSE209746 的 94→93 例差异(P-0044433 缺 RNA)建议在论文方法中注明
- [x] ~~生存分析事件数不足~~ → 已通过多变量/分层/非线性加深,见六.3;pCR 亚组零复发为新发现
- [ ] 联合模型的跨平台迁移需要更稳健的归一化策略(如分位数匹配)或限于同平台验证

---

## 六、深化分析(2026-09-01 下午)

### 1. MSS-only 敏感性分析(GSE209746, 剔除 5 例 MSI + 3 例 MSI-NA → n=85)

- MSI 与疗效:MSI 患者 CR 3/5(60%) vs MSS 23/85(27%),Fisher P=0.14, OR=3.97(趋势一致,样本少不显著)
- **IFN 关联不依赖 MSI**:MSS-only 中 IFN g=0.349(全队列 0.275),P=0.199,AUC=0.567——效应量不降反升
- IGF2(g=-0.29, P=0.077)、L1CAM(g=-0.50, **P=0.0054**)在 MSS-only 中与全队列基本一致
- **DESeq2 MSS-only 复现 15 基因:14/15 FDR<0.05 且方向一致,15/15 方向一致**;此前全队列唯一未复现的 CPS1 在 MSS-only 中显著(log2FC=-3.07, FDR=1.8e-4);NHLH2 方向一致(padj=NA,独立过滤)
- 脚本 `GSE209746_MSS_only.R`,输出 `GSE209746_MSS_only_features.csv` / `GSE209746_MSS_only_15genes.csv`

### 2. 免疫分组 IG1-IG4 方法学澄清 + MSS-only 验证

- **日志原"免疫热聚类 IG3/IG4 64% vs IG1/IG2 23%"实为论文(Nature Med, PMID 35970919)官方分组体系**:IG4=全部 dMMR/MSI(n=5)单独成组,MSS 内层次聚类得 IG3(免疫热)/IG2(中)/IG1(冷)
- 用论文 Table S1D 的 22 个 Bindea 免疫 ssGSEA 分数按论文方法复现:**IG3/4 CR=64%(7/11) vs IG1/2 24%(19/79),P=0.012**——与日志值一致
- **MSS-only(IG3 vs IG1/2):CR 67%(4/6) vs 24%(19/79),Fisher P=0.043**——剔除 MSI 后免疫热→完全缓解的富集仍然显著
- kmeans 种子稳定性(50 种子):IG3 CR 率中位 67%,IG1/2 中位 24%,中位 P=0.043,32/50 种子 P<0.05——结论稳健
- 脚本 `GSE209746_immune_groups.R`,输出 `GSE209746_ig_groups.csv` / `GSE209746_ig_groups_summary.csv`

### 3. 免疫细胞四队列 meta 分析(Hedges g, DL 随机效应)

| 细胞 | 合并 g (95%CI) | P | I² |
|---|---|---|---|
| NK cells | 0.39 (0.12, 0.66) | 5.1e-03 | 19.9% |
| Cytotoxic lymphocytes | 0.33 (0.08, 0.58) | 9.4e-03 | 9.0% |

- 四队列方向全部一致(好反应组高),与 IFN meta(g=0.44)构成完整机制链定量证据
- 脚本 `immune_meta_4cohorts.R`,输出 `immune_meta_4cohorts.csv` / `immune_meta_forest.png`

### 4. 生存分析加深(GSE87211)

- **多变量 Cox(DFS, n=175, 34 事件)**:校正 ypT/ypN 后 IFN HR=1.05(P=0.80);ypT 是最强预后因子(HR≈2.0/级, P=0.001)。NK(HR=1.16)、Cyto(HR=0.90)同样不独立显著
- **非线性检验**:自然样条 ns(df=3) vs 线性,LR P=0.44 → 线性假设成立(附图 `GSE87211_IFN_DFS_spline.png`)
- **pCR 分层:pCR 亚组 32 例 DFS 事件为 0**(随访期内零复发);non-pCR 亚组(n=143, 34 事件)IFN HR=1.01, P=0.97
- OS 多变量(n=175, 21 事件):IFN HR=0.94, P=0.80
- **结论修正:IFN 评分的弱预后信号经由"更高的病理缓解率"介导,它不是独立于 yp 分期的预后因子**;pCR→极好 DFS 是该队列的强事实
- 脚本 `GSE87211_survival_v2.R`,输出 `GSE87211_survival_v2_summary.csv`

### 5. GSE35452 表型定义核对(已解决)

- 原文:**Watanabe T, Kobunai T, et al. Dis Colon Rectum 2014;57(1):23-31 (PMID 24316942)**,东京大学;46 例芯片训练集 + 16 例 qPCR 独立验证集(后者不在 GEO)
- 方案:50.4 Gy/28 次 + **UFT(300-500mg/d)/LV(75mg/d)**(日本方案,非 5-FU 输注/卡培他滨),放疗结束后手术
- 反应定义:切除标本组织病理学检查,按**日本大肠癌学会(JSCCR)半定量退缩分级**;**逐例分级未公开,GEO 仅有作者二分标签(R=24/NR=22),逐例核对不可行,只能引用定义体系**
- 提示:该队列化疗方案(UFT/LV)与其余队列(5-FU 基础)不同,可能是 15 基因信号稀释的原因之一,建议讨论中注明

### 6. 肿瘤纯度/间质校正(MCP Fibroblasts + Endothelial)

| 队列 | IFN~response 未校正 | 校正间质后 | response~IFN 校正后 OR (P) |
|---|---|---|---|
| GSE209746 | β=0.28, P=0.23 | β=0.15, P=0.51 | 1.18 (0.52) |
| GSE35452 | β=0.24, P=0.42 | β=0.22, P=0.43 | 1.34 (0.42) |
| GSE87211 | β=0.47, P=0.011 | β=0.49, **P=0.005** | 1.68 (0.010) |
| GSE150082 | β=0.89, P=0.005 | β=0.88, **P=0.006** | 3.39 (0.013) |

- 间质本身与疗效无显著关联(四队列 P 均 >0.09);IFN 与间质相关中等(r=0.15-0.46)
- **结论:在 IFN 显著的两个大队列中,校正间质含量后关联不减反略增——IFN-疗效关联不是肿瘤纯度/间质比例的伪影**
- 脚本 `tumor_purity_adjust.R`,输出 `tumor_purity_adjust.csv`

### 7. 深化后总结论(修订版)

1. 干扰素/免疫激活主题进一步加固:MSS-only 成立(排除 MSI 混杂)+ 间质校正后成立(排除纯度伪影)+ NK/细胞毒性 meta 定量合并(P=5e-03/9e-03)
2. 机制链更新:免疫热亚群(IG3/4)→ CR 率约 2.7 倍(64% vs 24%),**MSS-only 仍 67% vs 24%(P=0.043)**
3. **预后层面定位明确:IFN 经由病理缓解间接关联预后,非独立预后因子;pCR 患者 DFS 零复发(GSE87211, 32 例)**
4. 15 基因 MSS-only 复现度 14/15 显著 + 15/15 方向一致,优于全队列

### 8. 新增队列 GSE45404(第5队列, GPL570, Mandard TRG)

- 42 例(Responder 19 / Non-responder 23),与文献计数一致(PMC10134325)
- IFN g=+0.39(P=0.21, 方向一致)、Cyto g=+0.46(P=0.13, 方向一致)、NK g≈0(本队列无信号)
- 15 基因方向:9/14 一致(INHBB 显著 P=0.010),符合 TRG 宽松表型队列的稀释规律
- 脚本 `GSE45404_analysis.R`,输出 `features_GSE45404.csv` / `GSE45404_stats.csv`

### 9. 五队列 meta 分析(IFN/NK/Cyto, Hedges g, DL 随机效应)

| 特征 | 合并 g (95%CI) | P | I² | 方向一致队列数 |
|---|---|---|---|---|
| IFN signature | 0.44 (0.21, 0.66) | **1.1e-04** | 0% | 5/5 |
| NK cells | 0.32 (0.06, 0.58) | 1.4e-02 | 24.5% | 4/5 |
| Cytotoxic lymphocytes | 0.35 (0.13, 0.57) | **1.6e-03** | 0% | 5/5 |

- 加入 GSE45404 后 IFN 合并 P 从 6.9e-04 强化到 1.1e-04,异质性降为 0
- 脚本 `meta_5cohorts.R`,输出 `meta_5cohorts_all.csv` / `meta_forest_5cohorts.png`(三面板森林图,主图候选)

### 10. TCGA-READ 预后旁证(非新辅助, 164 例治疗初治原发肿瘤)

- OS:IFN HR=0.92(P=0.70),校正分期后 HR=0.93(P=0.76);分期强预后(HR=2.1, P=0.003)
- PFI:IFN HR=0.90(P=0.53);**NK 细胞 HR=0.70 (0.51-0.97), P=0.031——NK 浸润对无进展生存有独立保护趋势**
- 结论与 GSE87211 一致:IFN 方向保护但不显著;NK 的 PFI 信号为机制链提供独立队列的预后旁证
- 脚本 `TCGA_READ_survival.R`,输出 `TCGA_READ_ifn_survival.csv` / `TCGA_READ_IFN_OS_KM.png`
- 方法备注:ssGSEA 基于样本内秩次,count→log2(CPM) 等单调变换不改变结果

### 11. 主图重绘(2026-09-01 晚, ggplot2, 修复文字重叠)

- 问题:旧版 base-R 森林图队列标签与置信区间须线/数值文字互相重叠,不宜投稿
- 重绘方案(ggplot2 + gridExtra,脚本 `figures_main.R`):
  - 森林图改用显式数值 y 坐标,队列名放 y 轴,g(95%CI) 与 P 值放入图外右侧专用文字列(`coord_cartesian(clip="off")` + 右 margin),合并效应用橙色菱形单独一行
  - KM 曲线按 `rep(names(km$strata), km$strata)` 正确映射分层,补删失点(shape=3)
  - 组合图改 3 行布局(A 全宽 → B 全宽双森林 → C+D 双栏),避免宽面板挤压
- 输出:**Fig1_main.png**(3800×4000, 300dpi, 主图终版)+ fig_A/B/C/D 单面板
- 已逐张目检:无重叠、无截断、无乱码(箭头 "→" 在 Windows png 默认字体下会乱码,统一写作 "to")

### 12. 论文初稿(2026-09-01 晚, 英文)

- 文件:`manuscript_draft.md`(完整 Original Article:标题/摘要/前言/方法/结果/讨论/声明/表1-2/图注/参考文献19条+编辑备注)
- 标题方向:Pre-treatment interferon activation and an immune-hot microenvironment characterize complete response to nCRT in LARC: a five-cohort transcriptomic meta-analysis
- 所有数据均取自本日志六节各小节精确数值;参考文献全部经 eutils 核实 PMID(Bindea 24138885、MCP-counter 27765066、GSVA 23323831、camera 22638577、OPRA 35483010、pROC 21414208、GSE87211 原文 29119627、GSE150082 原文 32784964、PMC9558072→PMID 36248969)
- 待作者补:署名/单位/基金/贡献声明;投稿前按目标期刊调参考文献格式

---

## 七、审稿防御分析 + 系统性队列扩充(2026-09-01 晚)

### 1. meta 分析诊断(脚本 `meta_diagnostics.R`)

- **Leave-one-cohort-out**(10 队列版):IFN 合并 g 范围 0.29–0.38,**P 恒 ≤0.027**——无任何单队列驱动;NK g 0.17–0.30(剔除 GSE87211 时 P=0.11,权重使然);Cyto g 0.18–0.30(同)
- **Egger 检验**:IFN P=0.21,NK P=0.53,Cyto P=0.11——无小样本效应迹象(k=10 功效有限,如实注明)
- 输出 `meta_leave_one_out.csv` / `meta_egger.csv` / `fig_S_funnel.png`

### 2. 跨队列免疫热二分验证(脚本 `immune_hot_xcohort.R`)

- 方法:各队列内按特征分数**上三分位**定义免疫热,Fisher 检验反应率富集,log-OR 做 DL 随机效应合并
- **10 队列结果**:NK-热 OR=2.00 (1.33–3.01) **P=8.6e-04**, I²=12%;Cyto-热 OR=1.61 (1.10–2.34) P=0.013, I²=1%;IFN-热 OR=1.53 (0.98–2.37) P=0.060
- 发现队列 IG3/4 的 CR 富集由此升级为**跨队列可复现**的关联
- 输出 `immune_hot_xcohort.csv` / `fig_S_immune_hot.png`

### 3. GSE87211 临床增益 + KRAS 校正 + 随访(脚本 `GSE87211_clinical_utility.R`)

- **IFN 在临床分期之外有独立增量**:pCR ~ cT+cN AUC=0.600;+IFN 后 AUC=0.655,似然比检验 **P=0.016**(DeLong P=0.21,n=197 功效有限);全模型中 ifn_z z=2.37, P=0.018
- **KRAS 校正后 IFN 仍显著**(pCR ~ IFN+KRAS: IFN P=0.0088);KRAS 突变者 pCR 率偏低趋势(12% vs 20%, P=0.17)
- **中位随访 71.2 个月**(reverse KM, DFS 队列 n=175)——pCR 零复发不是短随访伪影
- 输出 `GSE87211_clinical_utility.csv`

### 4. GSE209746 基因组学混杂校正(脚本 `GSE209746_genomic_adjust.R`,论文 Table S1C WES)

- WES+RNA 交集 n=57:**TMB、WES 纯度、fraction_cna、genome_doubled 与疗效均无关联**(P=0.16/0.52/0.95/0.16)
- 校正纯度+TMB+CNA 后 IFN~response β=0.26(未校正 0.35),response~IFN OR=1.40——方向不变、幅度略降(子集功效有限)
- **WES 真纯度**(非间质代理)本身与疗效无关 → 纯度伪影论证闭环
- 输出 `GSE209746_genomic_adjust.csv`

### 5. 15 基因逐基因 meta(脚本 `genes15_meta.R`)

- 14/15 基因可评(VAX1 在 vsd 中因低计数缺失 + GPL570 无探针);**方向一致 11/14**
- **INHBB 合并 g=-0.52, P=3.2e-4, FDR=0.0045, I²=33%**;KLK6 g=-0.32, P=0.004, FDR=0.03——两基因 FDR 显著
- 输出 `genes15_meta_percohort.csv` / `genes15_meta_pooled.csv` / `fig_S_genes15_meta.png`

### 6. 系统性队列检索与扩充(PRISMA 记录: `cohort_search_record.md`)

- GEO DataSets 系统检索(rectal + neoadjuvant/chemoradiotherapy + 表达谱 + 人)**38 命中 + 引文追踪 2 → 评估 40**
- **新纳入 5 队列**:GSE119409 (n=56, IFN g=+0.74)、GSE133057 (n=33, -0.16)、GSE53781 (n=26, +0.07)、GSE94104 (n=40 活检, +0.19)、GSE213331 (n=18, -0.58)
- **GSE94104 TRG 方向陷阱**:GEO 标签 1/2/3,初判 TRG1=好导致方向全反;原文(Alderdice, Mod Pathol 2017, PMID 28621318)明确 **TRG4=完全缓解(升序制)**,校正为 TRG3 vs TRG1-2 后方向转正——**改判须有原文依据,已记录**
- **排除 30 个并逐条记录原因**,其中:GSE216616 (n=298, JAMA Netw Open 2023, Akiyoshi) 与 GSE109057 (n=91) **疗效标签未寄存**——该研究独立报道了细胞毒性淋巴细胞→疗效的同一结论并在 GSE87211/GSE45404 验证,作趋同证据引用;GSE40492 (n=245) 与 GSE87211 同源(GRCSG)判为近重复;GSE318268 (Galunisertib 试验) 与 GSE242786 标签不可得;GSE123390 平台注释不可获取;GSE139255 仅 784 基因面板
- 平台注释攻坚战:GPL6102 annot ✔;GPL18134 (CodeLink) 无 symbol → org.Hs.eg.db ACCNUM 映射(覆盖 38%);GPL14951 → illuminaHumanv3.db(90%);GPL17586 不可得
- NCBI 触发了 reCAPTCHA 限流(并行下载所致)——后续下载一律串行+延时

### 7. 十队列 meta 最终结果(脚本 `meta_all_cohorts.R`, 输出 `meta_all_cohorts.csv` / `fig_S_forest_allcohorts.png`)

| 特征 | 合并 g (95%CI) | P | I² | 方向一致 |
|---|---|---|---|---|
| IFN | 0.33 (0.11–0.55) | 3.2e-03 | 25% | 8/10 |
| NK | 0.26 (0.08–0.44) | 5.0e-03 | 0% | 8/10 |
| Cyto | 0.27 (0.09–0.46) | 3.5e-03 | 0% | 8/10 |

- **高注释覆盖敏感性分析**(剔除 GSE133057 51%/GSE53781 38% 两个低覆盖队列, 8 队列):IFN g=0.40 (0.16–0.63) **P=8.0e-04**;NK 0.28 P=0.009;Cyto 0.32 (0.13–0.52) **P=1.1e-03, 8/8 方向一致**
- 解读:两个反向队列均为小样本+低覆盖/小 n(GSE133057 n=33 覆盖 51%;GSE213331 n=18),随机效应模型已按 SE 降权;主结论在 10 队列全集与高覆盖子集同向成立

### 8. 十队列版主图

- `figures_main.R` 已更新读 `meta_all_cohorts.csv`,输出 `Fig1_main.png`(3800×4400, 300dpi)+ `fig_A_forest_IFN_all.png` / `fig_B_forest_immune_all.png`
- 面板 A/B 为 10 队列森林图,C(IG 分组 CR 率)与 D(pCR-KM)不变

### 9. 修订后总结论(十队列版)

1. IFN/免疫激活主题在 **10 队列 596 例**中稳健成立(g=0.33, P=0.003, I²=25%, 8/10 同向;高覆盖子集 g=0.40, P=8e-04)
2. NK/细胞毒性淋巴细胞同向成立(NK P=0.005, Cyto P=0.003, I²=0),并有**二分免疫热的跨队列复现**(NK-热 OR=2.0, P=8.6e-04)
3. 混杂控制链完整:MSI(MSS-only)× 间质(MCP)× WES 真纯度 × TMB/CNA × KRAS(GSE87211)全部不能解释该关联
4. 临床定位:IFN 在临床分期之上有独立增量(LR P=0.016);预后经由缓解介导(pCR 零复发,中位随访 71 月);独立 298 例研究(JAMA Netw Open 2023)与本结论一致
5. 论文题目改为 ten-cohort 版,手稿 `manuscript_draft.md` 已整体改写为十队列版(2026-09-02):标题/摘要/方法(系统检索+5 新队列+注释路线+新分析)/结果(10 队列 meta、LOO/Egger、免疫热 OR 合并、临床效用、KRAS/基因组校正、15 基因 meta)/讨论(三条独立收敛证据:Akiyoshi n=298、Alderdice NK 签名、Rajamanickam CMS1↔CR)/Table 1(10 队列)/Table 2(新估计+高覆盖敏感性+免疫热 OR)/图例(Fig1 10 队列 + S0–S6)/参考文献 26 条(新增 8 条,PMID 均经 eutils 核实)

## 八、投稿前优化包(2026-09-02 上午)

### 1. 95% 预测区间 (IntHout 2016, BMJ Open; PMID 27406637, 已核实)
- 新脚本 `meta_prediction_intervals.R`: PI = mu ± t(k-2)·√(τ²+SE²); 输出 `meta_prediction_intervals.csv`
- 连续特征 (g):
  - IFN 全集 PI = **-0.15 ~ 0.81**(跨 0,由两个小样本/低覆盖反向队列驱动); 高覆盖子集 -0.09 ~ 0.88
  - NK 全集 PI = **0.05 ~ 0.48**(不跨 0); Cyto 全集 PI = **0.06 ~ 0.49**(不跨 0); 高覆盖 Cyto 0.08 ~ 0.57
- 免疫热 OR: NK PI 0.98~4.11(临界不跨 1), Cyto PI 1.01~2.56(不跨 1), IFN PI 0.61~3.84
- 论文叙事相应调整: NK/细胞毒性分数定位为"最可移植的定量信号", IFN 为通路层面总括; PI 跨 0 主动报告并归因
- `figures_main.R` 的 dl_pool 加 PI 计算, forest_gg 在合并菱形下画 PI 横杠+端帽+标注; I² 文本下移至 y=0.02 防压线; Fig1 已重绘并目验通过

### 2. Figure S0 PRISMA 流程图
- 新脚本 `fig_S0_prisma.R` (grid 绘制), 输出 `fig_S0_prisma.png` (3000×2100, 300dpi), 已目验
- 数字与 cohort_search_record.md 一致: 38+2=40 评估 → 排除 30 (标签未寄存 4 / 平台注释 2 / 重复 1 / 设计不适用 23) → 纳入 10 队列 596 例; 脚注注明 TCGA-READ 仅旁证

### 3. 手稿更新 (manuscript_draft.md)
- 摘要 353 → **299 词**(不含关键词): 方法段去掉队列罗列, 结果段压缩混杂清单, 结论段精简; PI 写入摘要 (NK 0.05-0.48 / Cyto 0.06-0.49)
- 方法 2.5 加 PI 计算方法 + IntHout 引文 [27]
- 结果 3.4 加 PI 段; 3.5 免疫热加 PI
- Table 2 加 "95% prediction interval" 列 + 敏感性/免疫热脚注补 PI
- 讨论 4.1 加 PI 解读段
- 图 1 图例注明菱形下橙杠 = 95% PI
- 参考文献 27 条 (新增 IntHout 2016, PMID 27406637 经 eutils 核实)
- 编辑备注更新 (S0 已绘 / PI 已加)

### 4. 可复现性打包
- `run_all.R`: 30 个正式脚本的全流程主控 (debug_*.R 除外), 分 5 阶段, 带计时; 阶段 0 数据获取为手工 (README 记录)
- `session_info.txt`: 加载全部 12 个关键包后的 sessionInfo 落盘 (R 4.6.1 / Win11 x64)
- `README.md`: 项目概览 + 队列表 + 复现方法 + 数据获取说明 + 脚本→产物映射 + 关键方法决策

### 5. 遗留事项(需用户行动)
- 给 GSE216616 通讯作者 (Akiyoshi 组) 发邮件索取疗效标签: 拿到后 10→11 队列, n=894, 且"收敛证据"升级为"纳入证据"
- GitHub/Zenodo 代码托管 (DOI) 需用户账号操作; 本地打包已就绪

## 九、投稿前技术审计与修正（2026-09-02）

> 本节替代第七至第八节中的“10 队列 596 例”和 DL 主分析口径。历史记录保留用于追溯，不再作为当前结果引用。

### 1. 队列与疗效标签修正

- **GSE213331**：原文和 GEO 对应研究为新辅助化疗，不是 nCRT。因此从 nCRT 主分析中排除，仅保留为全新辅助治疗敏感性分析。
- **GSE87211**：203 个直肠肿瘤中，GSM2325248 的 ypT/ypN 均缺失，不能自动归为 non-pCR；疗效分析改为 202 例（pCR 35、non-pCR 167）。GSM2325371 虽 ypT 缺失，但 ypN=2，可明确归为 non-pCR。
- **GSE94104**：处理脚本统一为该论文的升序 TRG 定义，TRG3=good，TRG1–2=poor；加入 11/29 的断言，避免未来重跑反转标签。
- 当前样本口径：nCRT 主分析 9 队列、577 例；加入 GSE213331 后的敏感性分析 10 数据集、595 例；高注释覆盖 nCRT 敏感性分析 7 队列、518 例。

### 2. 统计模型升级

- 主模型由 DerSimonian-Laird + 正态近似改为 **REML τ² + modified Hartung-Knapp**，方差缩放因子下限为 1；DL 结果保留为敏感性分析。
- 预测区间使用 t(k−2)；诊断限定于 9 个 nCRT 主分析队列。
- 漏斗图边界已改为围绕合并效应居中，而不是错误地围绕 0 居中。
- top-tertile 分析改称 **score-high sensitivity analysis**，不再称为独立验证；仅在存在零单元格时应用 0.5 连续性校正。

### 3. 审计后核心结果（独立数值复核；待 R 全流程重跑确认输出文件）

| 分析 | IFN | NK | Cytotoxic |
|---|---|---|---|
| nCRT 主分析（9 队列，n=577） | g=0.38 (0.16–0.61), P=0.0041, PI 0.16–0.61 | g=0.23 (−0.005–0.47), P=0.054, PI −0.12–0.59 | g=0.26 (0.04–0.49), P=0.026, PI 0.003–0.52 |
| 全新辅助治疗敏感性（10 数据集，n=595） | g=0.35 (0.10–0.59), P=0.011 | g=0.25 (0.03–0.48), P=0.032 | g=0.27 (0.06–0.48), P=0.018 |
| 高注释覆盖 nCRT（7 队列，n=518） | g=0.45 (0.20–0.70), P=0.0046 | g=0.26 (−0.03–0.54), P=0.072 | g=0.32 (0.07–0.57), P=0.020 |

- nCRT 主分析方向一致：IFN 8/9，NK 7/9，Cytotoxic 7/9。
- Leave-one-cohort-out：IFN g=0.34–0.43，所有删一结果 P=0.0036–0.020；NK 和 Cytotoxic 在部分删一分析中不再显著。
- Egger P：IFN 0.58、NK 0.24、Cytotoxic 0.050。队列数少、功效低；Cytotoxic 的名义信号需谨慎。
- nCRT top-tertile score-high：NK OR=1.98 (1.13–3.44), P=0.022, PI 0.70–5.54；IFN P=0.056；Cytotoxic P=0.090。
- 15 基因 REML-mKH 复核：无基因通过 BH 校正；INHBB P=0.024、KLK6 P=0.046，二者 FDR 均约 0.324，只能视为探索性结果。

### 4. 方法命名与论文措辞修正

- 实际实现为“使用 MCP-counter marker sets 的 ssGSEA scores”，不是 MCPcounter 包的丰度去卷积。
- “pre-specified high-coverage analysis”改为“sensitivity analysis”，因为没有带时间戳的预注册方案。
- 临床 T/N + IFN 结果改称同队列内增量关联；AUC 未作 bootstrap 或交叉验证乐观偏倚校正。
- 删除“独立于所有混杂”“预后由疗效介导”“carried by”等因果化表述；生存结果只描述为 response-linked association。
- 补充材料清单已在 `supplementary_manifest.md` 中逐项映射，并明确哪些图表需重跑或组版。

### 5. 复现性状态

- 新增 `meta_utils.R`、`validate_outputs.R`，并将 `run_all.R` 改为可移植路径、带时间戳日志和失败中止。
- 主图脚本增加 SVG/PDF 矢量输出；PRISMA 脚本增加 SVG 输出。
- 当前系统找不到 `Rscript`。因此本节所述代码已静态核对、核心结果已独立数值复核，但审计后的 R 全流程、图件重绘和视觉 QA 尚未完成。工作区内 2026-09-02 上午生成的旧 PNG/CSV 不应作为最终投稿文件。

## 七、工程可复现性优化（2026-09-03）

- 计数口径以 `project_config.R` 为准：9 个 nCRT 队列 577 例，加 GSE213331 18 例，共 595 例；早期段落中出现的 596 是旧版计数，不能用于最终报告。
- 新增 `project_config.R`，集中维护 10 个队列的样本数、good/poor 计数、治疗范围和三类分析集；`validate_outputs.R`、`meta_all_cohorts.R` 与 `meta_diagnostics.R` 统一读取该注册表。
- `validate_outputs.R` 现在额外检查样本 ID 的空值/重复、IFN/NK/Cyto 的有限数值和完整 schema，并将旧版额外列记录为警告；IGF2/L1CAM 允许因平台注释覆盖不足而缺失。
- `run_all.R` 增加 preflight：在运行任何分析前检查 R 包、输入文件和流水线脚本；`--preflight-only` 可只做环境检查。
- `run_all.R` 以独立子环境执行每个步骤；步骤间只通过文件传递结果，避免前一步残留对象掩盖缺失输入。
- `combined_model.R` 将训练专用特征写入 `combined_model_features_*.csv`，不再把训练结局 `y` 混入标准 `features_*.csv`。
- `new_cohorts_process3.R` 统一处理 GSE53781、GSE94104 和 GSE133057；`new_cohorts_process1.R` 不再先行计算 GSE133057，`run_all.R` 也不再调用早期的 `new_cohorts_process2.R`。

## 九、GSE216616 标签获取途径核查(2026-09-02)

### 公开途径逐一验证——均不可行
1. **GEO series matrix**(重新下载,11.5KB 元数据):逐样本字段仅 tissue/age/sex 三项,确认无 TRG/疗效/生存
2. **论文补充材料**(NCBI 被反爬挡 → Europe PMC supplementaryFiles 成功,s001.pdf 13 页):仅聚合 TRG 分布 (TRG1 99 / TRG2 113 / TRG3 31 / TRG4 50, 合计 293) 与散点图,**无 GSM→TRG 逐样本映射**
3. **BioProject PRJNA894509 / SRA**:仅测序元数据
4. **散点图数字化**:无样本 ID 映射,不可行
5. **样本文件名内码**(如 6C886T1):内部编码,无公开对照表

### 结论
唯一途径 = 联系通讯作者。通讯作者: Takashi Akiyoshi, MD, PhD (Cancer Institute Hospital, JFCR, Tokyo), takashi.akiyoshi@jfcr.or.jp
→ 邮件草稿见 `email_GSE216616_request.md`(英文正文 + 中文说明 + 无回复备选方案)

### 附带收获
- 论文补充材料已下载存档: `new_cohorts/akiyoshi_suppl/` (eTable 临床特征 + eFigure 1-8)
- 补充材料证实其 TRG 口径为升序制 (TRG3/4 = good responder, 81/293 = 27.7%), 与本研究二分框架兼容
- GSE216616 表达矩阵为逐样本 RSEM 文件 (GSE216616_RAW.tar), 拿到标签后处理管线与 GSE213331 类似 (counts → log2CPM → ssGSEA)

### 后续进展
- 2026-09-09: 索取邮件已由用户发送 (Bizhi Wei, 15619056250wbz@gmail.com → takashi.akiyoshi@jfcr.or.jp)
- 拿到标签后的接入预案: GSM+TRG 两列表 → 按 GSE213331 管线处理 RSEM (GSE216616_RAW.tar) → 生成 features_GSE216616.csv (TRG3/4=good) → meta_all_cohorts.R 自动扫描纳入 → 重跑 meta_prediction_intervals.R / immune_hot_xcohort.R / figures_main.R → 手稿数字与 Table 1/2 更新为 11 队列 n=894

## 十、署名落定与零版面费投稿策略(2026-09-09)

- **一致性核查**: 9月7日重跑的 meta_pooled_results.csv / score_high_pooled.csv / validation_summary.csv 与手稿(技术审计版)Table 2 数字逐项一致(9队列 577 例, IFN g=0.38 P=0.0041 PI 0.16-0.61; NK P=0.054; Cyto P=0.026; NK-high OR=1.98 P=0.022; 全队列 PASS)——旧 PNG/CSV 已全部被审计后版本覆盖
- **手稿署名定稿**: Bizhi Wei 独立作者(一作+通讯), Puai Medical College, Shaoyang University; Funding = None; ICMJE 贡献声明已填
- **导师确认的约束**: 无版面费预算、无 IF/分区要求 → 期刊名单重做, 纯 OA 刊(BMC/MDPI/Frontiers/PLoS/Sci Rep)全部排除, 只留订阅制发表免费的 hybrid 刊
- **投稿策略 v2** (`投稿策略讨论稿.docx`): 首投 Radiotherapy and Oncology (Elsevier, IF≈4.9) → Clinical Colorectal Cancer / Dis Colon Rectum → J Cancer Res Clin Oncol / Int J Colorectal Dis → 保底 Strahlenther Onkol / Cancer Biomarkers
- 建议: 注册免费 ORCID; 留意 hybrid 刊选订阅制路线即零费用

## 十一、导师意见落实轮(2026-09-09 下午)

导师 10 条意见逐条处理; 手稿为本节口径。

### 1. 效应方向统一(意见#1)
- 固定句 "**positive Hedges g indicates higher scores in responders**" 贯穿: 摘要开头 / 方法 2.5(加粗定义,并指明 IFN=primary feature, CTL/NK=supporting) / 结果 3.4 / Table 2 表头下 / 图 1 图例 / 森林图 x 轴(figures_main.R xlab "> 0: higher in responders") / 稳健性图轴 / 两张 harmonization 表 direction 列
- 森林图文本列 -0.00 → 0.00(fmt_g 函数防 sprintf 负零)

### 2. IFN 主线叙事(意见#2)
- 方法 2.5 明确 IFN=primary, CTL/NK=supporting; NK 二分明确写"single pre-chosen threshold (tertile) without any cut-off optimization"
- 主图面板 B 改为 CTL 左(标 B)/NK 右; 图例顺序同步; 讨论 4.1 明确 NK 只作 supportive sensitivity evidence, 不做 primary claim
- 未做任何 NK 最佳切点搜索(按导师意见拒绝)

### 3. IFN 签名稳健性分析(意见#3) — 新脚本 ifn_robustness.R
- 替代签名: Ayers IFN-γ 6 基因(JCI 2017, PMID 28650338, eutils 核实) + Reactome IFN-α/β (R-HSA-909733) + Reactome IFN-γ (R-HSA-877300), 基因集来自 Enrichr Reactome_2022 (reactome2022.gmt)
- 单基因: CXCL9/CXCL10/IDO1/STAT1/IFNG
- 逐队列加载器与各原脚本逐行一致(含 GSE87211 三态判定与 35/167 断言、GSE94104 升序 TRG 与 11/29 断言)
- 结果(9 nCRT 队列, REML-mKH): Ayers g=0.29 (0.03-0.55) P=0.034; Reactome α/β g=0.35 (0.12-0.58) P=0.008; Reactome γ g=0.33 (0.07-0.60) P=0.019; STAT1 g=0.39 P=0.006; CXCL10 g=0.32 P=0.017; IFNG g=0.28 P=0.019; CXCL9/IDO1 正向不显著 → **主发现不依赖基因集选择**, 支持稳定 IFN 激活状态
- 与 primary IFN 分数 Spearman ρ 0.58-0.92
- 产物: ifn_robustness_percohort/pooled/correlations.csv + fig_S_ifn_robustness.png(手稿 Fig S7 / Table S8)

### 4. 跨平台透明化(意见#4) — table_S2_platform_harmonization.csv
- 逐队列: 平台/n/表达来源/标准化/探针→基因路线/注释覆盖率/打分方法/方向约定; 总行声明"raw matrices never merged across cohorts, effect-size level 合成"; 方法 2.3 加同义句并引 Table S2

### 5. 疗效终点 harmonization 表(意见#5) — table_S1_response_harmonization.csv
- 逐队列: 原论文定义 / GEO 元数据来源 / 最终编码 / good/poor n / 人工核查备注(含 94→93、pCR-then-relapse、GSM2325248 排除、GSE94104 升序陷阱及代码断言、GSE213331 化疗类型); 方法 2.2 引 Table S1

### 6. 影响诊断补全(意见#6) — 新脚本 meta_baujat.R
- Baujat 图(x=Q 贡献, y=对合并效应影响)三特征并列, 图内标注 Q/df/τ²/I²; IFN 最大 Q 贡献 GSE150082, NK/Cyto 最大影响 GSE87211, 均无翻案性影响
- Table 2 增 τ² 和 Q(df),P 两列; 结果 3.4 原文报告 τ²=0/Q=8.12/df=8
- **PI≈CI 疑点解释(导师点名的统计审稿风险)**: τ²=0 时 PI=mu±t(k-2)·SE_mu 与 CI=mu±t(k-1)·sqrt(max(1,HK))·SE_mu 数值上几乎重合; 手稿明示三位小数 CI 0.161-0.606 vs PI 0.156-0.611

### 7. Egger 措辞(意见#7)
- 改为 "no strong evidence of small-study effects, although power was limited because fewer than 10 cohorts were available"; 全文无"无发表偏倚"式断言

### 8. 介导措辞降级(意见#8)
- 摘要: "consistent with, not proof of, response as an intermediate clinical correlate"; 结果 3.9 保持 response-linked patterns 框架并报告 pCR 0/32 vs 34/143; 全文无 formal mediation 声明

### 9. 预测模型定位(意见#9)
- 维持 in-sample incremental association 表述(审计版已做), 不新增任何预测模型

### 10. GSE216616 独立验证叙事(意见#10)
- 讨论 4.5(ii) 明确: 若拿到标签, 首选作为"大型独立验证队列"检验 9 队列效应的方向与大小是否样本外复现, 而非简单并入; 日志第九节接入预案改为独立验证设计

### 署名
- ORCID 0009-0008-9481-3024 已入作者块与通讯作者行; 邮箱按用户原文 gmai.com 填写(见"遗留")

### 其他
- 摘要重写后 305 词(含方向定义句与稳健性结果)
- 补充材料重编号: Table S1/S2=两张 harmonization 表, 原 S1-S5 顺移为 S3-S7, 新增 Fig S7/Table S8; supplementary_manifest.md 同步
- 新增参考文献 [28] Ayers 2017 (PMID 28650338, eutils 核实卷期页 127(8):2930-2940)
- run_all.R 尚未追加 ifn_robustness.R / meta_baujat.R(待下一轮工程整合)

### 遗留
- **邮箱拼写待用户确认**: 用户两次提供 15619056250wbz@gmai.com; 手稿已按原文填 gmai.com, 但致 Akiyoshi 的信函 docx 用的是 gmail.com —— 若实际为 gmail.com 需改手稿通讯作者行(投稿通讯地址错误会导致编辑部失联)
- Fig S2 组版(LOO+漏斗+Baujat 三合一)待做

## 十二、审计轮②: GSE87211 方案 + GSE46862/GSE56699 复核(2026-09-10)

### ① GSE87211 治疗方案审计 — 方案异质性明确但可解释
- 363 样本(含正常粘膜),直肠肿瘤 203: 5-FU+RT 111 / 5-FU+Oxaliplatin+RT 87 / 5-FU+Oxa+Cetuximab+RT 5
- **方案与 IFN 无关**(IFN~方案 P=0.53 未调 pCR / 联合模型 pCR 项 P=0.007 而方案项 NS) → 骨干不同(±奥沙利铂)不构成 IFN-疗效关联的混杂; 5 例 Cetuximab 亚组太小单列无意义
- 手稿 Table 1/方法需补一句: GSE87211 内部方案异质(5-FU±Oxa±Cetux), 方案与 IFN 无关
- 注意: 早期"GRCSG 300例 5-FU 单药"印象不准确——实际 GEO 记录明确为三方案

### ② GSE46862 — 维持排除, 但排除理由需改
- 原记录"原代细胞培养(primary cancer cell)"**属实**: source_name 全部为 "primary cancer cell", 69 样本
- 但重新审计发现: **该系列其实有逐样本 Dworak TRG 四分类**(MI=grade1 n=9 / MO=grade2 n=32 / NT=grade3 n=11 / TO=grade4 n=17), Shafi et al. Radiat Oncol 2016 (PMID 27005571), 治疗前活检 77 例建 163 基因多分类模型
- 平台 GPL6244 (Affymetrix Human Gene 1.0 ST) 注释可得 → **若接受"原代细胞(活检组织原代培养)表达谱"为肿瘤转录组, 该队列可纳入**; 保守处理仍可排除, 但记录须改为"样本为原代培养的 primary cancer cell(非组织块直接提取), 与其余队列的 bulk 转录组不同质, 排除"; 是否纳入待导师定夺
- 若纳入: TO(17) = good vs MI+MO(41) = poor(Dworak 4 vs 1-2, NT=3 居中或并入 poor), 方向: Dworak 升序 TRG4=完全缓解, 与 GSE94104 同口径

### ② GSE56699 — 维持排除, 理由需改且实质不变
- 72 样本: PB(治疗前活检) 58 + SS(手术后) 14; Isella et al. "Stromal contribution"系列(都灵/布鲁塞尔, FFPE, GPL14951 Illumina HT-12 v3, 与 GSE94104 同平台同型号)
- **有逐样本疗效**: 3 分类 CR 34 / RES 28 / PR 10 + Mandard TRG 1-5(活检内: TRG1 n=19, 2 n=8, 3 n=8, 4 n=17, 5 n=4, NA 2)
- Mandard 降序制: TRG1=完全缓解。活检内可二分: TRG1-2=good(27) vs TRG3-5=poor(29)
- 原记录未收录该队列(检索式中 rectal[Title] 不命中——标题是 "colorectal cancer (rectal samples)", 引文追踪亦未及) → **这是检索漏网, 须在 PRISMA 记录补上**
- 可纳入性: 平台注释路线与 GSE94104 相同(illuminaHumanv3.db 90% 覆盖), 58 例治疗前活检 → **技术上可纳入为第 10 个 nCRT 队列**; 需导师确认后执行

### 处置决定(待导师)
- GSE46862: 默认排除(原代培养不同质), 记录已改写; 若导师要纳入, 半天可出
- GSE56699: 建议纳入(PB 58 例, Mandard TRG1-2 vs 3-5, 与 GSE94104 同注释管线) → 9→10 nCRT 队列; 46862 不动
- 无论纳入与否, cohort_search_record.md 与 PRISMA 计数都要更新(56699 是漏网, 需透明披露)

---

## 十三、GSE56699 官方文件复核、患者层面重建与治疗范围校正（2026-09-10）

> 本节以官方 ArrayExpress/BioStudies SDRF/IDF 和源论文核验结果，**取代第十二节把 GSE56699 暂拟为 nCRT 队列的判断**。GSE56699 是术前放疗队列，不进入 9 队列 nCRT 主分析。

### 1. 官方数据与可审计输入

- 官方条目：E-GEOD-56699 / GEO GSE56699；源论文 Isella et al., *Nat Genet* 2015，PMID 25706627。
- 下载并纳入预检的输入：`new_cohorts/E-GEOD-56699.processed.1.zip`、`new_cohorts/E-GEOD-56699.sdrf.txt`、`new_cohorts/A-GEOD-14951.adf.txt`。
- SHA256：processed ZIP `9789A2A1A2DFC9D282EC8BBE83B06DD6A69DF7DEF1E65A4A61C4EEA3D8C77798`；SDRF `45AF7A5987B5447BC6B600F77C7F60DB0F5A50B3F61753690DA584F2C4AB37B2`；ADF `E50B8B4375ED0996FF530C832F14DBF5F5255C2667DFDFE8396616D0D460D553`。
- 旧归档中的 `GPL14951.annot.gz` 只有 990 bytes，内容为 HTML 404，不是平台注释；当前工作区的无效副本已删除，正式流程改用官方 ADF。

### 2. 独立样本口径

- SDRF 共 72 张芯片：治疗前活检（PB）58 张，手术标本（SS）14 张。
- 58 张 PB 芯片只对应 52 名患者；有 5 名患者存在重复活检（其中 1 名 3 张、4 名各 2 张）。
- PB Mandard 分级：TRG1=19 张、TRG2=8、TRG3=8、TRG4=18、TRG5=4、缺失=1。
- 患者层面排除唯一 Mandard 缺失者后 n=51：TRG1-2（good）23 人，TRG3-5（poor）28 人。
- 正式方案先在表达矩阵上按患者平均重复活检，再做 ssGSEA；替代方案先逐芯片评分再按患者平均。两种方案的患者分数 Spearman 相关：IFN 0.994、NK 0.998、Cyto 0.999。

### 3. 单队列与合并结果

- GSE56699 单队列 Hedges g：IFN −0.317、NK −0.186、Cyto −0.157，均与 9 队列 nCRT 主方向相反。
- 9 队列 nCRT 主分析完全不变：IFN g=0.384（95% CI 0.161–0.606，P=0.0041）；NK g=0.234（P=0.054）；Cyto g=0.264（P=0.026）。
- 加化疗 GSE213331（10 队列、n=595）：IFN g=0.347（P=0.0107），NK g=0.253（P=0.0321），Cyto g=0.270（P=0.0178）。
- 加放疗 GSE56699（10 队列、n=628）：IFN g=0.300（95% CI 0.037–0.564，P=0.0299，I²=33.9%，PI −0.241–0.841）；NK P=0.106；Cyto P=0.058。
- 全部 11 数据集（n=646）：IFN g=0.256（95% CI −0.021–0.532，P=0.0662，I²=41.0%）；NK g=0.204（P=0.0681）；Cyto g=0.220（P=0.0418，PI 跨 0）。
- 解释：GSE56699 提示治疗方案可能是效应修饰因素，也可能混有平台/终点/队列差异；不能据单个放疗队列下因果结论。主张必须限定为 nCRT。

### 4. 注册表、复现与 PRISMA 修订

- `project_config.R` 现定义 5 套分析集：`primary_ncrt9`、`chemotherapy_containing10`、`radiotherapy_containing10`、`all_neoadjuvant11`、`high_coverage_ncrt7`。
- `new_cohorts_process4.R` 从官方 ZIP/SDRF/ADF 重建 `features_GSE56699.csv`、患者标签映射和重复活检顺序敏感性表；注释覆盖 97.5%。
- PRISMA 透明记为 38 个预设检索命中 + 2 个初始引文追踪 + 1 个检索后复核发现 = 41 个评估；排除 30；纳入 11（9 nCRT 主分析 + 1 放疗敏感性 + 1 化疗敏感性）。
- GSE46862 仍排除：材料是原代癌细胞培养，不是直接 bulk 肿瘤组织；其标签存在与否不改变设计不合格结论。

### 5. 最终端到端重建与审计（2026-09-11）

- 使用 R 4.6.1 执行 `run_all.R`：2026-09-11 10:40:37 开始，10:48:45 结束，退出码 0；完整日志为 `run_all_20260911_104037.log`。
- preflight：14 个 R 包、28 项必需输入、36 个流水线脚本、3 个支持脚本全部通过；`validation_summary.csv` 十一队列全部 PASS。
- 流水线内部 warning 计数为 0；R 启动阶段另有 4 条 `LC_*` locale 设置警告，不影响分析执行和输出生成。
- 全流程结束后独立执行 `python audit_meta.py`，结果为 `AUDIT_ASSERTIONS_PASS`；`git diff --check` 退出码 0。
- Figure S1、S2、S4、S6 正式组版已目视检查，无裁切、重叠或不可读标签；其余流水线图表亦已刷新。Table S6/S7 的底层分析 CSV 已生成，但投稿版汇总表仍需人工组装与标签化。

### 6. Table S6/S7 投稿版自动汇总与最终重跑（2026-09-11）

- 相关模型脚本现额外保存完整、机器可读的模型摘要，包括模型实际样本量/事件数、效应量、95% CI、P 值、校正变量、适用总体和限制说明；原有简版输出保持不变。
- `assemble_supplementary_tables.R` 自动生成 `table_S6_confounder_adjustment.csv`（30 行：MSI、间质、WES、KRAS 四部分）和 `table_S7_survival_models.csv`（17 行：GSE87211/TCGA-READ，DFS/OS/PFI）。pCR 亚组 0 个 DFS 事件，明确标为 Cox 模型不可估计，未伪填 HR 或 P 值。
- 双工作表投稿审阅版为 `outputs/01a08bc2-2052-7d60-8183-6fe862fb032d/supplementary_tables_S6_S7.xlsx`；两张表均通过数值类型、公式错误扫描和目视检查。
- 更新后的 `run_all.R` 共 37 步：2026-09-11 12:11:36 开始、12:17:52 结束，退出码 0；日志 `run_all_20260911_121136.log`。preflight 14 个包、28 项输入、37 个流水线脚本、3 个支持脚本全部通过；流水线内部 warning 计数为 0。

## 十四、导师十项清单执行(2026-09-12)

### 已完成
- **① GSE87211 方案审计**: 直肠肿瘤 203 = 5-FU+RT 111 / 5-FU+Oxa+RT 87 / 5-FU+Oxa+Cetux+RT 5; **方案与 IFN 无关**(lm IFN~方案 P=0.53, 联合模型中 pCR 项 P=0.007 而方案项 NS)。方案异质性写入手稿方法
- **② GSE46862/56699**: 46862=Shafi Radiat Oncol 2016 (PMID 27005571), Dworak 4分类, 但样本为原代培养细胞 → 维持排除, 记录改写; 56699 已在日志十三节正式纳入(放疗敏感性, 反向 g=-0.32)
- **④ IFN 模型冻结**: `IFN_model_card.md` + `model_inputs_sha256.txt`(17 文件 SHA256: 基因集/工具函数/9队列 features/合并估计), 定义了 external validation 盲测协议与失败处置(不得改定义/换基因集/找切点)
- **⑥ 真实 MCPcounter 包敏感性** (`mcp_package_sensitivity.R`, 本地缓存官方 genes 签名): Cyto g=+0.31 (0.05-0.57) P=0.025 同向显著; NK g=+0.21 P=0.072(与 ssGSEA 版 0.054 相当); 两实现每队列 spearman 0.93-0.95 → **结论不依赖 ssGSEA 实现**
- **⑦ IFN 替代打分** (`ifn_scoring_sensitivity.R`, 队列内 mean-z): g=+0.38 (0.14-0.61) P=0.0066 vs ssGSEA 0.38 P=0.0041; 每队列 spearman 中位 0.94 → **打分方法不敏感**
- **⑧ bootstrap + 序贯** (`gse87211_bootstrap_sequential.R`): ΔAUC(表观)=0.055, **optimism 校正后 0.049 (95% CI -0.058~0.136)** — 增量真实但 CI 含 0, 与"关联非临床效用"定位一致; DFS 序贯(ypT,ypN→+IFN→+age,sex): IFN HR 1.05-1.07 全程 NS, c-index 0.675→0.713 主要由年龄贡献
- **⑨ 文献更新**: 2025-26 直肠 nCRT 标志物文献 1017 篇扫描; 关键重叠论文 = **iScience 2026 (PMID 41550766) 6 数据集 ML 186 基因签名** — 讨论需引用并区分(本研究=系统检索+免疫主题 meta, 非单基因集 ML); 新增引用候选: JITC 2026 (42547264) 基质-免疫交互预测, JCO PO 2026 (42566730), Med Oncol 2026 (42565918)
- **⑤-a Grampian 表达侧完成**: 223 CEL(与 S:CORT 样本表 223/223 匹配) → oligo RMA(pd.clariom.s.human 强制几何, 27189 探针集×223) → clariomshumantranscriptcluster.db 映射(TC#####.hg.1→SYMBOL, 71.7% 覆盖) → symbol 矩阵 19287×223 → **冻结 IFN 分数已落盘 scort/grampian_ifn_score.csv**(min 1.40 / med 1.67 / max 1.97)

### 进行中/阻塞
- **⑤-b Grampian 疗效标签**: pCR/TRG **不在** S:CORT 公开导出内(只有表达 CEL + 队列名样本表)。BJC 2026 (PMID 42613364, Roxburgh CS 通讯, Aberdeen) 用了 Grampian 表达+疗效 → 标签在 Aberdeen 团队手中。需邮件申请(可并入给 Akiyoshi 的批次)。ARISTOTLE 试验标签在 Lancet Oncol 2026 (42508427) 论文/试验注册处, 也非逐样本公开
- **③ ARISTOTLE 1.48GB** 后台重下中(约 70MB/1.48GB)
- S:CORT 技术要点(记录防再踩): Xcel 卡 CEL 用 oligo::read.celfiles(filenames=..., pkgname="pd.clariom.s.human") 强制读入; RMA 对 ExpressionFeatureSet 无 target 参数; 特征名已是 TC#####.hg.1 样式直接走 clariomshumantranscriptcluster.db

### 下一步
- Grampian/Roxburgh 邮件草稿(索取 Grampian pCR + ARISTOTLE 对照臂标签, 引 BJC 2026 与 Lancet Oncol 2026)
- ARISTOTLE CEL 到位后同管线处理(同为 Xcel 平台)
- ⑩ 手稿 narrative: 加 6/7/8 三个稳健性结果 + S:CORT 外验框架(表达侧就绪, 标签为申请中)

### 十四(续): ⑩ narrative 更新落稿(2026-09-12)
- 摘要 Results 重写: 稳健性从"3 签名+3 基因"扩为 "3 签名+3 基因+mean-z 打分+官方 MCPcounter 实现"(g=0.21-0.39 全部主比较 P<0.035); 加入 bootstrap 校正 ΔAUC 0.049 (CI -0.058~0.136)
- 方法 2.3 末尾加两项实现检查(mean-z / MCPcounter 包)与 2000 次 Efron-Gong bootstrap、生存序贯框架描述
- 结果 3.5 改题为 "...across alternative signatures, scoring methods and deconvolution implementations", 加入 6/7 两项数值
- 结果 3.8 加 bootstrap 校正句; 3.9 加 GSE87211 序贯调整(c-index 0.675→0.680(IFN)→0.713(年龄), IFN HR 全程 NS)
- 讨论 4.1 加"two scoring methods and two deconvolution implementations"; 4.2 加 iScience 2026 [29] 与 JITC 2026 [30] 区分句; 4.5 加 S:CORT Grampian/ARISTOTLE 冻结打分+预注册外验声明
- 补充清单加 Table S9(两项实现敏感性); 参考文献 30 条([29] Corrò iScience 2026, [30] Wang JITC 2026, 作者/卷期页经 Europe PMC 核实)
- 待 ARISTOTLE zip 下完(1.48GB)后: 同管线(oligo+pd.clariom.s.human+tc db)处理 → 与 Grampian 分数一起等标签

### 十五、等待中的外部依赖(截至 2026-09-12)
1. **Akiyoshi 组** (GSE216616 n=298, 已发邮件 2026-09-09)
2. **Roxburgh 组/Aberdeen** (Grampian pCR + ARISTOTLE 对照臂标签; 草稿 email_Grampian_ARISTOTLE_request.md, 发送前需从 BJC 2026 论文或 Aberdeen 官网核实邮箱)
3. 两批标签任一到达 → 按 IFN_model_card.md 的预注册协议出外验结果(方向/标准化差/AUC, 无切点优化)

### 十六、LaTeX 投稿版手稿(2026-09-12)
- `manuscript.tex`(pdfLaTeX, TeX Live 2026 编译) + `manuscript.pdf`(14 页)
- 版式: article 单栏 11pt / a4 / 2.5cm 对称边距 / lineno 行号(审稿用) / natbib 数字引用 / booktabs
- 参考文献 31 条按首次出现自动重编号 —— **修复了 markdown 版 [29] 撞号**(Isella 与 Corrò 均为 29)
- Table 1: adjustbox 包 tabular(p{6.8cm} 响应定义列), 宽高均受控; Table 2: adjustbox tabular + 表下脚注块(minipage, 因合并超页拆出 float)
- 质检: 0 错误 / 0 Overfull≥10pt / pdf_qa 仅剩 cover 全出血提示(期刊投稿稿无需封面, 忽略)
- 目验: 标题页、Table 1/2 页、Figure 1 页全部正常; 图 1 直接嵌入 Fig1_main.png
- 重编译命令: `pdflatex manuscript.tex` ×2 (D:\GEO2R 下)

## 十七、GSE213331 撤出全部数据集(2026-09-12, 导师审查发现)

### 证据(GEO 官方记录, 重新核验)
- `!Series_overall_design`: "resected tissues from nine LARC patients with complete pathological response (pCR), nine without pathological response (npCR) and **biopsy tissues from nine LARC patients prior nCRT**"
- `!Sample_description` 分布: 9 Biopsy prior to nCRT(无疗效标签) + 9 Resected after nCRT with pCR + 9 Resected after nCRT without pCR
- 结论: 18 例带标签样本 = **nCRT 术后切除组织**; 唯一的治疗前样本(9活检)反而无标签。此前两轮审计均误判("pCR vs npCR 后新辅助化疗")——既错在时间点, 也错在治疗类型(实为 nCRT 非化疗)

### 处置
1. `project_config.R` 移除该队列, 版本号 `2026-09-12-gse213331-removed-v1`; 总样本断言 646→628
2. `features_GSE213331.csv` → `REMOVED_features_GSE213331_posttreatment.csv` 留档
3. `meta_all_cohorts.R` 标签表清理; 重跑: meta_pooled_results / PI / immune_hot(score-high) / LOO / Egger / Baujat / 主图 —— 全部成功
4. "化疗敏感性(10数据集)" 与 "全11数据集(646)" 分析组永久删除

### 影响评估(重跑后核验)
- **主分析 9 队列 n=577 完全不变**: IFN g=0.384 (0.161–0.606) P=0.0041 I²=1.5% τ²=0; NK P=0.054; Cyto P=0.026 —— 导师"好消息"判断证实
- 放疗敏感性(10数据集 n=628): IFN g=0.300 (0.037–0.564) P=0.0299 I²=33.9% —— 与十三节数字一致(该组本就不含 GSE213331)
- 高覆盖7队列: 不变
- score-high OR: 主分析不变(IFN 1.66/NK 1.98/Cyto 1.57)
- LOO/Egger/Baujat: 全部重算, 数值同前(仅9队列)
- bootstrap/序贯/稳健性(⑥⑦⑧/签名稳健性): 本来就只用9队列或GSE87211, 无需重跑

### 文档同步
- manuscript.tex: 16 处修改(摘要2/引言1/方法4/结果3/讨论1/数据可用性/表1行/表2脚注2行/S0图例/排除计数30→31); 重编译 14 页, 0 错误 0 超宽
- manuscript_draft.md: 同步同样修改
- cohort_search_record.md: 追加更正段; 排除计数 30→31(评估41/排除31/纳入10: 9 nCRT+1放疗敏感性)
- README/manifest: 队列表与状态同步
- 教训(写入记录): **疗效标签必须与样本时间点(treatment-naive)同时核验**, 不能只看标签字段; GEO 的 !Sample_description 常带组织时点信息, 应成为标准检查项

## 十八、导师逐句核对修订(2026-09-12 第二轮)

1. **41/30 矛盾**: Results 3.1 残留 "excluded 30" → 31(tex+md)。41评估/31排除/10纳入口径全线统一
2. **摘要 robustness 拆两层**: IFN 关联(3替代签名+3核心基因+mean-z, g=0.28-0.39 全 P<0.035) 与 cytotoxic 关联(MCPcounter g=0.31 P=0.025; NK 0.21 P=0.072 不确定) 分开表述——原文把 MCPcounter 的 NK 塞进 "all reproduced IFN" 统计上不成立
3. **bootstrap 旧句清理**: Methods 2.7 "lacked bootstrap correction" 与 Discussion "no internal validation corrected optimism" 均与 2000-resample Efron-Gong 校正结果矛盾 → 改为 "Bootstrap optimism correction yielded a similar point estimate ... CI included zero; not validated predictive performance"
4. **GSE213331 残留单数化**: "non-nCRT cohorts were combined"→"this cohort was combined with the nCRT set"; "radiotherapy or chemotherapy alone"→"radiotherapy alone"(已无化疗队列)
5. **命名统一**: "radiotherapy-containing" → "treatment-scope sensitivity (9 nCRT + RT-only)"(nCRT 本身含放疗, 原名不当); tex表2脚注/md表2/两处摘要
6. **deconvolution 拆分**: Results 3.5 标题改 "Robustness across IFN definitions, scoring methods and immune-scoring implementations"; 结尾拆两句(IFN 对基因集/打分法稳健; cytotoxic 被 MCPcounter 丰度估计复现); Discussion 4.1 同步——原文说 "IFN 不依赖 deconvolution" 逻辑不通(MCPcounter 检验没重算 IFN)
7. **ROC 4 vs 3**: 核实=3张ROC图(209746 0.556/35452 0.578/87211 0.601)+150082 仅数值0.764; 正文 "four cohorts" 数值正确, S4 图例加注 "GSE150082 AUC reported in text without plotted curve"
8. **pre-registered → locked+timestamped**: 无公开注册号, 措辞改为 "locked and timestamped ... with input-file hashes"(诚实)
9. **Statistics 补全**: 声明 REML-mKH 为项目自实现(meta_utils.R pool_random_effects, 依 IntHout 公式), 未用 meta 包; MCPcounter 包(ebecht/GitHub master)+官方签名
10. **文献审计(eutils核实)**: Chatila 补卷页 Nat Med 2022;28:1646-1655; Xue 补 Front Oncol 2022;12:993726; Corrò 题名去 "locally advanced"(正式题为 "...in rectal cancer"); Wang 题名补 "plus immunotherapy"

### 十九、ARISTOTLE 表达侧完成(2026-09-12)
- 1.48GB zip 下载完整(300 CEL, 分段续传多轮); 样本表内对照臂 121/121 全匹配(另 179 CEL 为伊立替康试验臂, 导出多给, 不用)
- 同 Grampian 管线: pd.clariom.s.human RMA (27189×121) → clariomshumantranscriptcluster.db (71.7%) → symbol 19287×121 → 冻结 IFN 分数落盘 `scort/aristotle_ifn_score.csv`
- 分布与 Grampian 相似(中位 1.68 vs 1.67), 平台/管线一致性旁证
- scort_aristotle_prep.R: 硬编码只取对照臂(样本表过滤), 试验臂 CEL 留在解包目录不参与
- 两个 S:CORT 队列的表达侧均已就绪, 仅等疗效标签(Roxburgh 邮件)

### 二十、英文 Word 版手稿(2026-09-12)
- `manuscript.docx` (Times New Roman 11pt / 1.3 行距 / A4 / 编号引用与 tex 同源): 由 make_manuscript_docx.py 从 manuscript.tex 直接转写, 保证 tex/docx 内容单源一致
- 结构: 标题块→Abstract(四段)→Introduction→Methods(9小节)→Results(9小节)→Discussion(7小节)→Statements→Table 1/2→Figure 1(嵌入 Fig1_main.png)→Supplementary legends(S0-S7+补充表)→References(31条)
- 生成器踩坑记录(供复用): tex 转写时正则必须 re.escape(chr(92)+cmd) 构造, 直接 r"\section" 会因 \s=\s空白类 失配; "\t" 类有效转义会破坏 minipage 匹配; \vspace 前缀需剥离才能命中 S1-S7 图例; \posg 宏需展开; Discussion 提取边界需止于 Statements 注释
- postcheck: 0 错误; 2 个无害警告(表格内 1.15 行距属排版惯例; 分页符空段为刻意分页)
- LibreOffice 渲染目验: 标题页/表 1/表 2 正常; word 版 17 页(单栏行距更大)
- 文件族: manuscript.tex + manuscript.pdf(LaTeX 14页) + manuscript.docx(Word) + manuscript_word.pdf(Word 版预览)

## 二十一、终点定义与 GSE87211 方案敏感性(2026-09-12)

### 设计(导师审定口径)
1. **Endpoint-class 亚组 meta**(仅 primary_ncrt9, 不改主分析)
   - complete-response endpoints: GSE209746(sustained CR) + GSE87211(pCR=ypT0N0), k=2, n=295
   - broader response-classification endpoints: 其余 7 队列(TRG 或 sensitivity/resistance), k=7, n=282
   - 输出: 各亚组 pooled g+95%CI+I²/τ²; Q_between(固定效应 Q 分解); 每队列方向
   - **严格终点亚组不报 prediction interval**(t 分布 df=k−2=0 不可估); k=2 的 mKH CI 仅作探索性解读
   - 措辞锁定: 不写 "no difference between endpoint types", 写 "positive under both classes, no evidence of effect modification; comparison underpowered"
2. **GSE87211 队列内 regimen sensitivity**
   - 字段 `preoperative radiochemotherapy (rct):ch1`
   - 分层: fluoropyrimidine_CRT (5-FU+RT, n=111) vs oxaliplatin_containing_CRT (5-FU+Oxa±Cetux, n=91; Cetux 并入)
   - 终点仍为 pCR=ypT0N0(与主分析一致)
   - 分层 Hedges g + z(score)~response×regimen 交互

### 实现
- `meta_utils.R`: `pool_random_effects(..., prediction=FALSE)` 允许 k≥2; prediction=TRUE 仍要求 k≥3
- `meta_endpoint_subgroups.R` → `meta_endpoint_effects_classified.csv`, `meta_endpoint_subgroups.csv`, `meta_endpoint_q_between.csv`
- `gse87211_regimen_sensitivity.R` → `gse87211_regimen_counts.csv`, `gse87211_regimen_effects.csv`, `gse87211_regimen_interaction.csv`
- 两脚本已挂入 `run_all.R` 阶段 4; 单独重跑退出码 0

### 结果
**Endpoint class (IFN)**
- complete-response: 两队列 g 均正 (0.275, 0.480); 合并点估计 g=0.399; mKH CI 极宽 (−1.45–2.25), P=0.22; I²=0; **无 PI**
- broader: g=0.368 (0.016–0.720) P=0.043, I²=21.4%, 6/7 方向为正
- Q_between: IFN P=0.889 / NK P=0.189 / Cyto P=0.264 → 无效应修饰证据, 但功效有限

**GSE87211 regimen**
- fluoro n=111, pCR 22 (19.8%); oxa-containing n=91, pCR 13 (14.3%)
- IFN g: fluoro 0.473 / oxa 0.540; 交互 P=0.70
- NK g: 0.472 / 0.580; 交互 P=0.79
- Cyto g: 0.333 / 0.905; 交互 P=0.19
- 结论: 方案分层内 IFN 方向一致, 交互不显著 → ±奥沙利铂不驱动主关联

### 手稿同步
- Methods 2.5: 补 endpoint-class 与 regimen 敏感性方法段(含 k=2 无 PI 声明)
- Results 新 3.5 "Response-definition and regimen sensitivity"; Discussion 4.1 补一句
- 补充表图例加 S10; Table S9 仍为 mean-z/MCPcounter 实现敏感性
- Discussion 4.3 残留 "no internal validation corrected optimism" → 改为 bootstrap CI 含 0 的表述(与十八节第3条一致)
