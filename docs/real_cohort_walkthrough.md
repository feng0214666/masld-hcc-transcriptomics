# GSE126848：从公开计数到真实差异分析

本项目已在本地执行导入与 QC，并在 GitHub Actions 的 Ubuntu 24.04 环境完成 DESeq2 分析。原始公开计数不随代码重复分发，仓库保留派生结果和来源校验和。

## 输入与样本映射

[GEO Series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126848) 提供基因计数和 family SOFT 元数据。表达矩阵的列如 `0869` 对应 SOFT description 中的 `869`。适配器只在匹配时去掉数字编号的前导零，验证一一对应且无冲突，并保存原始列名和 description；最终分析列名使用 GSM accession。

导入得到 19,786 个基因、57 个样本：NAFL 15、NASH 16、Normal-weight 14、Obese 12。分组保留原始研究标签。QC 的全队列过滤保留 14,190 个基因；PC1/PC2 解释方差分别约 24.12% 和 18.59%。

## 复现命令

在仓库根目录下载两个 NCBI 官方文件。Windows 可使用 `curl.exe`：

```sh
curl --fail --location --retry 3 --output data/raw/GSE126848_Gene_counts_raw.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126848/suppl/GSE126848_Gene_counts_raw.txt.gz
curl --fail --location --retry 3 --output data/raw/GSE126848_family.soft.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126848/soft/GSE126848_family.soft.gz
Rscript --vanilla analysis/import_gse126848.R data/raw/GSE126848_Gene_counts_raw.txt.gz data/raw/GSE126848_family.soft.gz data/processed/GSE126848
Rscript --vanilla analysis/run.R --dataset GSE126848 --counts data/processed/GSE126848/counts.tsv --metadata data/processed/GSE126848/samples.tsv --deseq2 --reference Normal-weight --case NASH --covariates sex --out results/GSE126848-de
```

输出目录必须不存在。若之前导入过，复用已有规范输入，或选择新的导入目录。DESeq2 安装说明见首页；本次快照使用 R 4.3.3、DESeq2 1.42.0，完整依赖版本见结果目录中的 sessionInfo.txt。

## 已完成的对比

- 模型为 `~ sex + group`，纳入 14 例 Normal-weight 和 16 例 NASH；NAFL 和 Obese 样本仅保留在全队列 QC 中。
- 对比内按 count≥10 且至少 3 个样本满足进行过滤，13,979 个基因进入模型，13,969 个得到非缺失调整 p 值。
- 按 BH 调整 p＜0.05 且 |log2FC|≥1，上调 981 个、下调 569 个。方向为 NASH / Normal-weight。
- 使用 DESeq2 默认 Wald 检验、独立过滤及离群值处理；LFC 未收缩。

此对比估计这些人群标签之间、调整性别后的表达差异。年龄、BMI、组织组成和其他未纳入因素仍可能影响差异；不能由这一次对比单独归因于某个分子机制。当前结果未进行 YAP 或脂质代谢签名的正式富集验证。

## 结果与来源

[真实 QC 报告](../examples/GSE126848-qc/report.html) 和 [真实 DE 报告](../examples/GSE126848-de/report.html) 可下载后离线打开。结果目录包括完整基因结果、模型矩阵、样本清单、来源映射、输入与代码校验和。

DE 快照来自 [GitHub Actions run 44106406192](https://github.com/feng0214666/masld-hcc-transcriptomics/actions/runs/44106406192)，代码提交为 `2af3857fbd3ea75cd62ea97e15bad4ee46e13ee3`。快照的 code_checksums.tsv 对应该提交的代码；后续报告界面改进不会追溯修改这份分析快照。CI 会继续用最新代码重新运行同一队列，并上传最新产物。
