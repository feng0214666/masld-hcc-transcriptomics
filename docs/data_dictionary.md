# 输入与输出约定

所有命令从仓库根目录执行。正式分析以一个队列为单位运行。

## counts.tsv

Tab 分隔；第一列必须叫 `gene_id`，其余列名必须与 metadata 的 `sample_id` 精确对应。基因和样本 ID 均非空且唯一。

```text
gene_id	S01	S02	S03	S04	S05	S06
GENE_A	12	19	24	49	58	67
GENE_B	80	94	81	22	38	21
```

矩阵值必须为有限、非负、未转换的整数计数。**数值符合整数约束并不能证明其生物学来源合格**：仍需检查原始数据说明。TPM、FPKM、CPM、对数表达和 Cuffdiff 推定数量不可通过四舍五入伪装成原始计数。已有连续表达矩阵需要另外实现与其量纲相容的分析方法。

Ensembl ID 与 HGNC symbol 不能混用；GMT 必须采用相同 ID 空间。不要在未记录规则时自动删版本后缀、转换 ID 或合并重复基因。表达文件若含额外注释列，先审查后显式导出规范 TSV。

## samples.tsv

```text
sample_id	group	patient_id
S01	adjacent	P1
S02	adjacent	P2
S03	adjacent	P3
S04	tumor	P1
S05	tumor	P2
S06	tumor	P3
```

`sample_id` 和 `group` 必填，配对时 `patient_id` 必填。ID 按字符读取并保留前导零。允许元数据顺序不同，程序按 ID 重排；缺少样本或多出样本会报错。协变量通过 `--covariates age,sex,batch` 明确指定。数值列作为连续变量，字符列作为分类变量；用数值编码的分类协变量应先转换为明确的文字类别。协变量缺失不自动删除样本。患者 ID 以去标识化公开编号或本地替代编号表示。

GSE164441_pairs.tsv 的 sample_id 为 GSM accession，source_title 为 GEO 标题。该文件提供公开记录中的配对关系，不保证补充表达文件使用相同列名。应显式保存 `expression_column → GSM → patient_id` 映射，检查无遗漏、无重复后导入。邻近组织不等于健康肝脏。

## GMT

每行：集合名、来源描述、至少一个基因 ID，以 Tab 分隔。集合名不可重复。建议将版本、来源链接、下载日期、正负方向定义与最终文件校验和保存在旁附说明中。本评分函数使用无方向的正向基因集均值；若签名含上调和下调成分，需分别定义和验证评分方法。

输出 `program_coverage.tsv` 的 requested 为集合去重后基因数，matched 为通过过滤且可计算 z 分数的重叠基因数，coverage = matched/requested。未达阈值的集合 `scored=FALSE`，不产生样本分数；全无重叠时导出仅有表头的分数表。

## 运行产物

| 文件 | 含义 |
|---|---|
| report.html | 可离线打开的结果索引和 QC 表格 |
| sample_qc.tsv | 全部输入基因的文库大小、非零表达基因数 |
| metadata_aligned.tsv | 与表达矩阵列顺序一致的元数据 |
| gene_filter.tsv | 每个基因是否通过过滤 |
| pca_coordinates.tsv / pca_variance.tsv | 样本 PC1/PC2 坐标和全部主成分解释方差 |
| sample_correlations.tsv | 过滤后 log-CPM 的 Pearson 样本相关 |
| program_scores.tsv / program_coverage.tsv | 可选描述性程序评分及覆盖率 |
| differential_expression.tsv | 可选 DESeq2 结果；case/reference 方向见运行参数 |
| model_matrix.tsv / model_samples.tsv / model_formula.txt | DE 实际设计矩阵、对比样本及公式 |
| de_gene_filter.tsv / de_summary.tsv | 仅指定对比的基因过滤决策与显式阈值摘要 |
| trait_associations.tsv / trait_sample_inclusion.tsv | 可选数值表型关联；方法、样本数、状态、效应、CI、原始和校正 p 值及缺失值记录 |
| input_checksums.tsv / code_checksums.tsv | 输入文件及运行脚本的 MD5，用于复核版本 |
| run_info.tsv / sessionInfo.txt | 命令参数、方法标签和实际软件会话 |

输出目录若已存在则拒绝运行，防止混入上一次分析产物。失败时保留的目录应检查原因后使用新目录重试。MD5 用作内容一致性标识，不作为安全签名。

## 数值表型关联

`--trait NAS` 指向一个数值列。未指定协变量时计算 Spearman rho；指定协变量时拟合 `score ~ trait + covariates`，estimate 为每单位表型对应的分数变化，CI 是 OLS 的 95% 区间。缺失或非有限的表型值被记录后排除，协变量缺失会阻止运行。

至少需要 5 个有效独立样本；常量、样本不足和不可估计设计记录为相应 status，并保留 NA。重复患者会报错。BH 仅校正该次单个表型运行中有可用检验的程序。多个表型、多次运行或探索性选择的整体多重比较策略需另行定义。
