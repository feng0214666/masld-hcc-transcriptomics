# 候选队列资料审查

核对日期：2026-11-18。已获取三个公开队列的表达和样本资料，复现入口及文件校验见[三队列分析](reproduce_three_cohort.md)。

## GSE126848 · discovery

[GEO source](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126848)

- 57 个人肝脏样本，记录包含正常体重、肥胖、NAFL、NASH 分组。
- 提供 `GSE126848_Gene_counts_raw.txt.gz`。
- 已读取 19,786 个唯一基因 ID，使用 SOFT description 对应表达列的数字编号完成 GSM 映射。匹配时显式处理前导零，保留原始编号；未做基因 ID 转换。
- 正常体重 14、肥胖 12、NAFL 15、NASH 16。四组共同拟合 `~ sex + group`，估计 NASH 相对三组的差异；[复现方法](reproduce_three_cohort.md)。

## GSE225740 · histology validation

[GEO source](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225740)

- 93 例肝活检，记录有组织学活动、纤维化和脂肪性肝炎相关信息。
- 提供 `GSE225740_NAFLD_RNA_counts.csv.gz` 与 `GSE225740_NAFLD_phenotypes.csv.gz`。
- 样本标题中存在缺失的组织学字段。保留为 NA；不能归入 0 级。
- GEO SOFT description、计数列、表型患者编号逐一匹配，表型数值与 SOFT 双向核对。NAS 为 93 例，纤维化为 87 例；排查未发现完全相同的表达向量。公开表型文件缺少可用于调整的年龄、性别和批次字段。

## GSE164441 · paired HCC

[GEO source](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164441)

- 20 个样本，即 10 对 NAFLD 相关 HCC 肿瘤和邻近组织。
- 提供 Cuffdiff 的 FPKM 和 count 文本摘要；GEO 注明原始测序数据因隐私原因未提供。
- 根据公开样本标题整理 [GSM 与配对表](../data/metadata/GSE164441_pairs.tsv)。表中 patient_id 是记录标题的样本编号，配对关系需在导入表达矩阵时再次核对。
- FPKM 表包含可对应到 10 对患者的 20 个逐样本字段。count 表首批数值明确包含小数，属于连续估计。分析输入使用逐样本 FPKM，保留配对设计与伪计数敏感性结果。

## 标签与研究解释

原始研究使用 NAFLD/NASH 术语。没有逐例代谢条件等资料时，不能声称已按现代 MASLD/MASH 标准重新诊断。跨队列差异还可能反映采样、人群、技术和组织组成等差异，应先在队列内估计效应，再设计综合策略。
