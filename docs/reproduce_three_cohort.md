# 复现三队列分析

在 Linux 或 macOS 的仓库根目录运行以下命令。Windows 可在启用了 Docker Desktop 集成的 WSL2 终端使用相同命令。分析容器固定为 Linux x86-64。需要 Docker、curl 和 sha256sum。

```bash
bash analysis/fetch_study.sh data/raw
docker build --platform linux/amd64 --tag masld-analysis environment
mkdir -p results
docker run --rm --network none --user "$(id -u):$(id -g)" \
  -v "$PWD:/study:ro" -v "$PWD/results:/output" -w /study masld-analysis \
  Rscript --vanilla analysis/run_study.R data/raw /output/three-cohort
```

已经安装相应 R 和 Bioconductor 包的环境，可以直接用 `bash analysis/reproduce_study.sh results/three-cohort` 下载、验证并运行。输出目录需要使用新路径，成功标记为 `_SUCCESS`。

[输入清单](../config/study_sources.tsv)锁定 11 个 GEO 和 MSigDB 文件的 URL 与 SHA-256。分析容器固定 Ubuntu 镜像摘要和 2026-11-17 UTC 的 Ubuntu 软件包快照，验证 R 4.3.3、DESeq2 1.42.0 和 limma 3.58.1。安装 TLS 信任证书使用 Ubuntu 基础镜像的软件源，分析依赖使用冻结快照。运行产物记录安装包清单、R 会话、代码 MD5 和输入来源。固定容器是三队列分析的复现入口。

## 模型与结果表

| 队列 | 方法 | 输出 |
|---|---|---|
| [GSE126848](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126848) | 57 例一次拟合 `~ sex + group`；从同一拟合提取 NASH / Obese、NASH / NAFL、NASH / Normal-weight | 每组分母和分子、差异基因效应量、Wald 95% CI、对比内 BH 值、程序秩检验 |
| [GSE225740](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225740) | 按患者 ID 核对 GEO SOFT、表型表和表达列；队列内 log2(CPM+1) 基因 z 分数均值；Spearman rho 与 2,000 次患者自助重采样的百分位区间 | NAS 与纤维化的有效样本、缺失或重复表达记录、rho、95% CI、联合 BH 值 |
| [GSE164441](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164441) | 10 对患者；Cuffdiff 20 列逐样本 FPKM；对 log2(FPKM+0.5) 拟合 `~ patient_id + group`，limma 经验贝叶斯趋势校正 | 单基因配对效应和 95% CI、程序秩检验、FPKM +1 的敏感性结果 |

发现队列只调整了可核实的性别字段。性别调整并不能消除采样点、体重和其他未提供协变量的影响。GSE225740 的公开表型 CSV 包含 NAS、纤维化和脂肪性肝炎编码，缺少可用于本模型的年龄、性别或批次信息，因此实施未调整等级相关。Bootstrap 区间以构建好的队列内程序分数为条件。对完全重复的表达列会记录并排除这组样本。HCC 首先保留唯一、可量化的 Ensembl 基因，要求 FPKM ≥1 至少出现在 5 个样本，再做配对拟合。公开 count 表包含小数估计值，原始测序读段未公开。

四个预定程序来自 [MSigDB 2025.1.Hs](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/) 的 H 和 C6 集合。脂质程序源自 Hallmark 集合，YAP 使用 [Cordenonsi 保守签名](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/CORDENONSI_YAP_CONSERVED_SIGNATURE.html)。项目保存集合 IDs、来源、版本和覆盖率，运行时从 Broad 下载成员表和固定 CHIP 映射。GSE126848 与 GSE164441 使用 Ensembl 映射，GSE225740 使用 NCBI Gene ID 映射。歧义 ID 和多对一映射不进入程序计算。[MSigDB 许可](https://gsea-msigdb.org/gsea/msigdb_license_terms.jsp)

发现队列用有符号 Wald 统计量、HCC 用配对模型的有符号 moderated t 统计量做 rank-based CAMERA-PR。主设置基因间相关系数 0.01，发现队列再用 0.05 检查敏感性。平均百分位秩减 0.5 用于表示程序排序的方向与幅度，数值不等同于 log2FC。发现队列 12 个程序检验、组织学队列 8 个检验、HCC 4 个检验分别形成各自的 BH 校正家族。[limma CAMERA-PR 方法说明](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf)

`cross_cohort_programs.tsv` 在同一张表中展示各队列的方向、效应、覆盖和 FDR。其 `estimate` 对等级相关是 rho，对排序检验是平均百分位秩偏移。病例对照、组织学关联与肿瘤配对使用不同的统计量。并列结果用于判断方向与证据是否一致，横断面数据无法建立个体随时间的进展顺序。
