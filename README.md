# MASLD–HCC Transcriptomics

**从脂质代谢到 Hippo–YAP：面向人肝脏公开队列的可复现转录组分析工具。**

Cross-cohort investigation of lipid-metabolic and Hippo–YAP transcriptional programs across MASLD severity and MASLD-associated hepatocellular carcinoma.

[![R workflow checks](https://github.com/feng0214666/masld-hcc-transcriptomics/actions/workflows/checks.yml/badge.svg)](https://github.com/feng0214666/masld-hcc-transcriptomics/actions/workflows/checks.yml)
[![R](https://img.shields.io/badge/R-%E2%89%A54.3-276DC3?logo=r)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-17876f.svg)](LICENSE)

本项目研究脂质代谢与 YAP/TEAD 相关转录程序在独立人肝脏队列中的变化。已完成 GSE126848 真实计数接入、57 个样本的 QC，以及调整性别的 NASH 对正常体重组差异分析；同时提供 GEO 导入、配对建模、协变量调整、程序评分、组织学关联和离线报告工具。

## 已完成的真实分析

| 数据与分析 | 已保存结果 |
|---|---|
| GSE126848：19,786 个基因、57 个样本 | [真实 QC 报告](examples/GSE126848-qc/report.html)、PCA 坐标与样本映射 |
| NASH 16 例 vs 正常体重组 14 例，设计 `~ sex + group` | [完整 DESeq2 结果](examples/GSE126848-de/differential_expression.tsv)、[模型矩阵](examples/GSE126848-de/model_matrix.tsv) |
| FDR＜0.05 且绝对 log2FC≥1 | 上调 981、下调 569；[结果摘要](examples/GSE126848-de/de_summary.tsv) |

这些数值来自指定队列和模型。[复现步骤与结果说明](docs/real_cohort_walkthrough.md)记录了输入来源、样本规模、软件版本与代码提交。程序签名的生物学验证及其他队列的独立复现仍是下一阶段工作。

![真实队列报告预览](docs/assets/report-preview.png)

## 立即运行

在仓库根目录执行；基础示例只用 R 自带包，无需下载表达数据。

```sh
Rscript --vanilla tests/run_tests.R
Rscript --vanilla analysis/run.R --demo --out results/my-demo
```

用浏览器打开 `results/my-demo/report.html`，即可查看样本组成、文库统计、PCA 方差与所有结果表的下载入口。每次运行请使用新的输出目录。

**示例包含 12 个模拟样本、6 对配对标签和 120 个模拟基因，仅用于验证软件。** 已运行的示例位于 [examples/demo](examples/demo)，其 [报告源码](examples/demo/report.html) 可下载后在浏览器打开。

## 现在能做什么

| 能力 | 实际实现 |
|---|---|
| GEO 元数据获取 | 下载 family SOFT，导出 GSM、标题、组织、物种及原始特征；保存来源、时间与 MD5 |
| 输入质量检查 | 检查计数类型、缺失值、重复基因与样本、样本对应关系及空文库 |
| 探索性 QC | 文库大小、检出基因数、表达过滤、可悬停查询的 PCA 散点图、解释方差、样本相关矩阵 |
| 基因程序评分 | 读取自备 GMT，计算队列内平均基因 z 分数，并输出覆盖率与跳过原因对应指标 |
| 差异表达 | DESeq2；显式指定 case/reference，支持配对、年龄／性别／批次等协变量及秩检查 |
| 组织学关联 | GMT 程序分数与数值表型的 Spearman 相关或协变量调整 OLS；记录缺失值与 BH 校正 |
| 可追溯结果 | 独立 HTML 报告、TSV 表格、输入校验和、参数和 `sessionInfo()` |
| 持续检查 | GitHub Actions 运行测试与离线示例，上传报告作为运行产物 |

```mermaid
flowchart LR
    A[GEO 元数据与来源记录] --> B[人工确认分组与表达列映射]
    B --> C[整数计数与样本检查]
    C --> D[队列内 QC 与 PCA]
    D --> E[GMT 程序评分]
    D --> F[可选 DESeq2 配对或非配对分析]
    E --> G[HTML 报告与可追溯表格]
    F --> G
```

## 公开队列与研究设计

核心问题：**脂质代谢和 Hippo–YAP 相关转录程序的变化，能否在不同 MASLD 严重程度及 MASLD 相关 HCC 的独立队列中重复观察到？**

| 队列 | GEO 记录样本数 | 用途 | 导入前需要确认 |
|---|---:|---|---|
| [GSE126848](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126848) | 57 | 正常体重、肥胖、NAFL 与 NASH 的比较 | 基因注释、表达列与 GSM 对应、原始分组 |
| [GSE225740](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225740) | 93 | NAS 与纤维化等组织学指标 | 保留缺失值，检查协变量与表型编码 |
| [GSE164441](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164441) | 20 | 10 对 NAFLD 相关 HCC 肿瘤／邻近组织 | Cuffdiff 文件结构与计数定义；原始测序数据未公开 |

样本数是 GEO 记录中的规模，并非本项目最终纳入数。旧研究保留原始 NAFLD/NASH 标签；MASLD/MASH 的现代临床分类需要相应个体资料支持。详见 [队列审查记录](docs/cohort_audit.md)。

## 接入自己的数据

先获取元数据，随后人工核对标签和表达列，避免按列顺序猜测样本身份：

```sh
Rscript --vanilla analysis/fetch_geo.R GSE126848 data/metadata/GSE126848_run1
```

准备以 `gene_id` 开头的计数 TSV，以及含 `sample_id`、`group` 的元数据 TSV；[输入规范和示例](docs/data_dictionary.md) 列出了完整要求。

```sh
Rscript --vanilla analysis/run.R --counts data/processed/counts.tsv --metadata data/metadata/samples.tsv --out results/cohort-qc
```

可选地通过 `--gmt config/programs.gmt` 加入自行核定来源与版本的基因集。演示基因集只包含 `DEMO_GENE_*`；项目尚未提供经过文献核验的 YAP 活性签名。

差异分析前先安装 Bioconductor 的 DESeq2：

```r
install.packages("BiocManager")
BiocManager::install("DESeq2")
```

```sh
Rscript --vanilla analysis/run.R --counts data/processed/counts.tsv --metadata data/metadata/samples.tsv --deseq2 --reference adjacent --case tumor --paired --out results/hcc-de
```

`log2FoldChange > 0` 表示 case 相对 reference 升高。配对分析要求每个受试者两组各一个样本，至少三个完整配对。其他组保留在 QC 中，差异分析仅使用指定两组。通过 `--covariates age,sex,batch` 显式加入协变量；缺失、常量、完全混杂和配对项共线会报错。模型公式、纳入样本、设计矩阵和对比专用过滤表随结果保存。

程序分数与组织学指标的关联示例：

```sh
Rscript --vanilla analysis/run.R --counts counts.tsv --metadata samples.tsv --gmt programs.gmt --trait NAS --covariates age,sex --out results/nas-association
```

不提供 `--covariates` 时计算 Spearman 相关；提供时拟合 OLS 并输出每单位表型对应的分数变化、95% CI 和 BH 校正 p 值。同一患者重复采样会被拒绝，需另行采用适当的重复测量模型。

## 方法与复现

- QC 使用 count ≥10 且至少 3 个样本满足的过滤规则，基于全部输入基因的文库大小计算 log2(CPM+1)，再进行中心化 PCA。
- 程序评分为队列内各基因标准化后的均值；至少匹配 3 个非恒定基因且覆盖率 ≥50%。这是描述性评分，不等同于 GSVA、ssGSEA 或富集显著性检验。
- DESeq2 使用选定对比的整数计数和相同过滤阈值，输出 Wald 检验及 BH 调整结果；独立过滤等因素产生的 NA 保留。未进行 LFC 收缩。
- 每个队列独立分析。程序关联的 BH 校正范围是该次运行、该表型成功检验的程序集合；多个表型的整体检验家族需另行预设。跨队列程序评分不能直接当作统一标尺比较；综合效应和异质性分析尚待实现。
- 基础工作流使用 `--vanilla` 避免首次启动自动下载 renv。原有 `renv.lock` 仅记录 renv，**尚未锁定可选 DESeq2 环境**。安装并验证完整分析依赖后，使用 `renv::snapshot()` 更新锁文件；每次运行同时保存实际 R 会话版本。

## 文件导航

```text
R/workflow.R                 输入检查、QC、评分、DESeq2 与报告函数
R/associations.R             数值表型关联、协变量调整与检验状态
R/import_geo.R               GSE126848 公开计数的显式编号映射
R/report_widgets.R           离线 PCA 交互视图
R/geo.R                      GEO SOFT 下载及解析
analysis/run.R               分析命令行入口
analysis/fetch_geo.R         GEO 元数据命令行入口
config/cohorts.tsv           已核对的候选队列概览
data/metadata/GSE164441_pairs.tsv  公开的 10 对样本对应表
docs/                        输入说明、方法和队列审查
examples/demo/               已运行的模拟示例与输出
examples/GSE126848-qc/       真实队列 QC 与来源映射
examples/GSE126848-de/       已完成的真实 DESeq2 分析快照
tests/run_tests.R            输入边界、设计、评分和解析检查
.github/workflows/checks.yml CI 与示例报告产物
```

[分析日志](analysis_log.md)记录研究决策；[研究方案](docs/study_protocol.md)说明后续验证顺序。代码采用 [MIT License](LICENSE)。

本项目输出为研究分析结果。横断面关联不直接证明纵向进展或因果关系；转录签名也不直接测量 YAP 核定位、蛋白活性或棕榈酰化。
