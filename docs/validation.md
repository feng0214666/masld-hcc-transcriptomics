# 验证记录

2026-11-18 的 [GitHub Actions 三队列运行](https://github.com/feng0214666/masld-hcc-transcriptomics/actions/runs/44291619378)成功完成。该作业从固定 Ubuntu 镜像和软件包快照构建容器，校验全部 11 个源文件的 SHA-256，并在禁用网络的分析容器内完成三队列拟合。R 4.3.3、DESeq2 1.42.0 和 limma 3.58.1 已通过运行时版本检查。

- 32 项基础检查和 29 项扩展检查在本地通过。CLI 组织学流程、示例报告、输入校验和及输出覆盖保护也通过。
- 容器内实际运行 DESeq2 的配对、反向对比和年龄调整测试，随后拟合 GSE126848 四组模型。
- 真实数据检查覆盖发现队列 57 例、组织学队列 93/87 例、HCC 10 对的样本映射，以及 24 项程序结果的覆盖率、区间和 BH 校正家族。记录在 `tests/check_study.R`。
- 从 GitHub 下载完整产物后核对归档 SHA-256，并用同一检查脚本复核文件。仓库中的 [结果快照](../examples/three-cohort/) 与该次运行的对应文件逐字节核对。

[历史两组 GSE126848 验证](real_cohort_walkthrough.md)保留其原有运行出处。当前四组模型使用独立的结果目录。多模型、测序单位和版本选择见[复现指南](reproduce_three_cohort.md)。
