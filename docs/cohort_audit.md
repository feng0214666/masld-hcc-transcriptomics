# 候选队列资料审查

核对日期：2026-11-17。依据为 NCBI GEO Series 页面，尚未完成所有表达文件和临床协变量的逐列审计。下述规模来自公开记录，不是本项目分析后的有效样本量。

## GSE126848 · discovery

[GEO source](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126848)

- 57 个人肝脏样本，记录包含正常体重、肥胖、NAFL、NASH 分组。
- 提供 `GSE126848_Gene_counts_raw.txt.gz`。
- 下一步：检查基因 ID、表达文件列名和 GSM 映射，记录组别、性别和可用批次字段。对照组的正常体重／肥胖标签分别保留。

## GSE225740 · histology validation

[GEO source](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225740)

- 93 例肝活检，记录有组织学活动、纤维化和脂肪性肝炎相关信息。
- 提供 `GSE225740_NAFLD_RNA_counts.csv.gz` 与 `GSE225740_NAFLD_phenotypes.csv.gz`。
- 样本标题中存在缺失的组织学字段。保留为 NA；不能归入 0 级。
- 下一步：核对表型文件编码、计数列顺序与可用协变量，再预设组织学相关分析。相关／回归模块尚未实现。

## GSE164441 · paired HCC

[GEO source](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164441)

- 20 个样本，即 10 对 NAFLD 相关 HCC 肿瘤和邻近组织。
- 提供 Cuffdiff 的 FPKM 和 count 文本摘要；GEO 注明原始测序数据因隐私原因未提供。
- 根据公开样本标题整理 [GSM 与配对表](../data/metadata/GSE164441_pairs.tsv)。表中 patient_id 是记录标题的样本编号，配对关系需在导入表达矩阵时再次核对。
- 下一步：检查 Cuffdiff 表是否具有逐样本表达、数值单位和估计量定义。如果只包含组均值／差异摘要，就不能还原个体矩阵或拟合配对模型；如为连续估计值，需采用相应表达量方法。当前通用整数计数入口不会尝试转换它。

## 标签与研究解释

原始研究使用 NAFLD/NASH 术语。没有逐例代谢条件等资料时，不能声称已按现代 MASLD/MASH 标准重新诊断。跨队列差异还可能反映采样、人群、技术和组织组成等差异，应先在队列内估计效应，再设计综合策略。
