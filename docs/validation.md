# 验证记录

日期：2026-11-17。本地环境：Windows，R 4.6.1。远端集成：Ubuntu 24.04，R 4.3.3，DESeq2 1.42.0。

## 基础与扩展检查

- tests/run_tests.R：32 项基础检查通过，覆盖非法计数、ID、元数据重排、QC、PCA、评分、配对和 SOFT 解析。
- tests/test_extended.R：29 项扩展检查通过，覆盖协变量、完全混杂、配对共线、字符 ID、组织学关联、缺失值、重复患者、显式 GEO 映射及报告转义。
- tests/test_association_cli.R：GMT → 调整表型关联 → HTML 输出的完整 CLI 集成通过。
- tests/check_report.R：模拟标签、报告链接、输入读回、校验和、外部 TSV 输入及输出覆盖保护通过。
- 模拟和真实 QC 报告均在本地浏览器渲染检查。

## DESeq2 与真实队列

- tests/test_deseq2.R 在远端要求安装 DESeq2，禁止跳过。实际拟合完整配对、反向对比和年龄调整模型，均通过。
- 真实 GSE126848 官方计数和 family SOFT 已下载成功；导入 19,786 个基因、57 个样本，明确检查前导零和 GSM 对应。
- 真实 NASH / Normal-weight 的性别调整分析成功，完整结果保存至 examples/GSE126848-de。代码快照与运行链接见该目录 README。
- 先前下载超时已通过成功获取官方文件和真实导入消除；独立 fetch_geo.R 的网络行为仍取决于运行环境。在线 CI 通过 curl 下载官方文件并验证适配器。

## 尚待研究工作

GSE225740 组织学实证分析、GSE164441 表达量纲适配、经过来源核验的生物学签名、跨队列综合和敏感性分析尚未完成。

原 renv.lock 仅记录 renv。实际 DESeq2 版本以结果目录 sessionInfo.txt 为准；可选 Bioconductor 环境尚未完整锁定。
