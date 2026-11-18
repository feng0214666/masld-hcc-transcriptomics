# MASLD–HCC Transcriptomics

**在三个人肝脏公开队列中检验脂质代谢和 YAP 相关转录程序。**

Cross-cohort investigation of lipid-metabolic and Hippo–YAP transcriptional programs across MASLD severity and MASLD-associated hepatocellular carcinoma.

[![R workflow checks](https://github.com/feng0214666/masld-hcc-transcriptomics/actions/workflows/checks.yml/badge.svg)](https://github.com/feng0214666/masld-hcc-transcriptomics/actions/workflows/checks.yml)
[![Three-cohort reproducibility](https://github.com/feng0214666/masld-hcc-transcriptomics/actions/workflows/three-cohort.yml/badge.svg)](https://github.com/feng0214666/masld-hcc-transcriptomics/actions/workflows/three-cohort.yml)
[![R](https://img.shields.io/badge/R-%E2%89%A54.3-276DC3?logo=r)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-17876f.svg)](LICENSE)

三个队列覆盖四组肝活检比较、疾病活动度与纤维化关联，以及肝癌肿瘤与邻近组织的配对比较。统一使用四个固定来源的表达程序。［[程序总表](examples/three-cohort/cross_cohort_programs.tsv) · [完整运行产物](https://github.com/feng0214666/masld-hcc-transcriptomics/actions/runs/44291619378/artifacts/14420874689) · [复现指南](docs/reproduce_three_cohort.md)］

## 已完成的真实分析

| 队列与问题 | 主要结果 |
|---|---|
| GSE126848 · 57 例，四组统一性别调整模型 | NASH vs 肥胖组有 823 个上调、494 个下调基因，NASH vs NAFL 有 63 个上调、7 个下调基因（FDR < 0.05，绝对 log₂FC ≥ 1）。[三个对比](examples/three-cohort/discovery_contrasts.tsv) |
| GSE225740 · NAS 93 例，纤维化 87 例 | YAP 签名与纤维化正相关，Spearman ρ = 0.46，95% CI 0.27–0.62，FDR = 5.96 × 10⁻⁵。[组织学结果](examples/three-cohort/histology_programs.tsv) |
| GSE164441 · 10 对肿瘤与邻近组织 | 脂肪酸代谢基因集在肿瘤中下向富集（FDR = 6.79 × 10⁻⁴），YAP 签名的 FDR 为 0.987。[配对结果](examples/three-cohort/hcc_programs.tsv) |

YAP 签名在 NASH vs NAFL 的主设置中上向富集（FDR = 0.00117），提高假定基因间相关性后 FDR 为 0.131。组织学和 HCC 对比的含义不同，三队列结果按各自问题解释。[方法与敏感性结果](docs/three_cohort_results.md)

## 复现完整分析

```bash
bash analysis/reproduce_study.sh results/three-cohort
```

此入口适用于已安装 R、DESeq2 和 limma 的环境。使用[锁定的 Docker 环境](docs/reproduce_three_cohort.md)可复现 GitHub Actions 的软件版本和三个分析。所有 11 个公开输入均核对 SHA-256。

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
| [GSE126848](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126848) | 57 | 正常体重、肥胖、NAFL 与 NASH 的比较 | 计数样本编号与 GEO SOFT 一一对应 |
| [GSE225740](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225740) | 93 | NAS 与纤维化等组织学指标 | 表型 CSV 与 GEO SOFT 相互核对，6 例缺少纤维化分期 |
| [GSE164441](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164441) | 20 | 10 对 NAFLD 相关 HCC 肿瘤与邻近组织 | 使用逐样本 FPKM，依据公开患者编号配对 |

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

可通过 `--gmt programs.gmt` 加入与输入使用相同 gene ID 的基因集。三队列分析另用[固定版本的四个程序](config/programs.tsv)和相应 ID 映射。

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

## 通用计数入口的方法

- QC 使用 count ≥10 且至少 3 个样本满足的过滤规则，基于全部输入基因的文库大小计算 log2(CPM+1)，再进行中心化 PCA。
- 程序评分为队列内各基因标准化后的均值；至少匹配 3 个非恒定基因且覆盖率 ≥50%。这是描述性评分，不等同于 GSVA、ssGSEA 或富集显著性检验。
- DESeq2 使用选定对比的整数计数和相同过滤阈值，输出 Wald 检验及 BH 调整结果；独立过滤等因素产生的 NA 保留。未进行 LFC 收缩。
- 通用程序关联的 BH 校正范围是单次运行的单个表型。三队列研究的校正家族和具体方法见[复现指南](docs/reproduce_three_cohort.md)。
- 基础工作流使用 `--vanilla`，适合独立运行通用脚本。三队列环境由 `environment/` 内的 Docker 文件锁定，程序同时保存实际 R 会话和安装包版本。

## 文件导航

```text
R/workflow.R                 输入检查、QC、评分、DESeq2 与报告函数
R/associations.R             数值表型关联、协变量调整与检验状态
R/import_geo.R               GSE126848 公开计数的显式编号映射
R/report_widgets.R           离线 PCA 交互视图
R/geo.R                      GEO SOFT 下载及解析
R/discovery.R                四组 DESeq2 对比
R/histology.R                组织学样本映射与程序相关
R/hcc_paired.R               Cuffdiff 审查与配对 limma 模型
R/program_analysis.R         固定基因集映射和秩检验
analysis/run.R               分析命令行入口
analysis/reproduce_study.sh  三队列复现入口
analysis/fetch_geo.R         GEO 元数据命令行入口
config/study_sources.tsv     GEO / MSigDB 下载源和 SHA-256
environment/                 锁定的 R / Bioconductor 容器
data/metadata/GSE164441_pairs.tsv  公开的 10 对样本对应表
docs/                        输入说明、方法和队列审查
examples/demo/               已运行的模拟示例与输出
examples/GSE126848-qc/       真实队列 QC 与来源映射
examples/GSE126848-de/       已完成的真实 DESeq2 分析快照
examples/three-cohort/      三队列程序结果、样本统计与运行出处
tests/run_tests.R            输入边界、设计、评分和解析检查
.github/workflows/checks.yml CI 与示例报告产物
```

[研究方案](docs/three_cohort_plan.md)说明对比与基因集选择。[GSE126848 交互 QC](examples/GSE126848-qc/report.html)可下载后查看样本 PCA。代码采用 [MIT License](LICENSE)。

本项目输出为研究分析结果。横断面关联不直接证明纵向进展或因果关系；转录签名也不直接测量 YAP 核定位、蛋白活性或棕榈酰化。
