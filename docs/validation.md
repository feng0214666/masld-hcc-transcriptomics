# 验证记录

日期：2026-11-17。环境：Windows，R 4.6.1。

## 已完成

- `Rscript --vanilla tests/run_tests.R`：32 项检查通过。覆盖输入非法值、重复及错误样本 ID、元数据重排、文库 QC、PCA 方差、评分覆盖率、配对设计、SOFT 解析、真实公开配对表和文件读取。
- `Rscript --vanilla analysis/run.R --demo --out examples/demo`：完整离线示例执行成功。
- `Rscript --vanilla tests/check_report.R`：报告所有文件链接存在，模拟数据标签存在，输入可读回，输入和代码 MD5 一致。

## 尚未验证

- 本地未安装可选 DESeq2，因此该集成测试明确跳过；不能将基础检查通过理解为真实差异表达已经完成。
- 在线 GEO family SOFT 下载在当前环境失败：沙箱内网络初始化失败，提权网络尝试超时。解析器已用固定样例验证，但在线获取仍需在可访问 NCBI 的环境中验证。
- CI 配置已添加；远端实际运行状态以 Actions 为准。
- 未执行真实队列的表达导入、差异表达、组织学关联或跨队列综合分析。

原始 `renv.lock` 继续保留，不声称它已覆盖可选 Bioconductor 依赖。
