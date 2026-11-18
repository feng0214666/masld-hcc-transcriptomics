# 三队列结果快照

由 [GitHub Actions run 44291619378](https://github.com/feng0214666/masld-hcc-transcriptomics/actions/runs/44291619378) 在提交 `ccd7598d403f967ca95c08601d11b097e15752dc` 上生成。代码 MD5 对应运行提交，数据来源 URL 和 SHA-256 见 [source_manifest.tsv](source_manifest.tsv)。

- [三队列程序总表](cross_cohort_programs.tsv)包含 24 项分析、方向、效应量、FDR、基因覆盖率和组织学区间。
- [发现队列差异摘要](discovery_contrasts.tsv)与[程序检验](discovery_programs.tsv)包含三个对比及基因间相关系数敏感性检查。
- [组织学关联](histology_programs.tsv)与[有效样本](histology_cohort.tsv)记录 NAS 和纤维化效应。
- [HCC 配对程序结果](hcc_programs.tsv)、[基因过滤规模](hcc_cohort.tsv)与[Cuffdiff 文件审查](hcc_source_format.txt)记录主分析与伪计数敏感性检查。
- [完整运行产物](https://github.com/feng0214666/masld-hcc-transcriptomics/actions/runs/44291619378/artifacts/14420874689)包含四张逐基因结果表、95% CI、三个样本映射、设计矩阵、软件版本与代码校验和。旧产物过期后可按[复现指南](../../docs/reproduce_three_cohort.md)重新生成。
