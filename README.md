# masld-hcc-transcriptomics
跨队列分析MASLD不同严重阶段及MASLD相关肝细胞癌中的脂质代谢和Hippo–YAP转录程序。Cross-cohort transcriptomic analysis of lipid metabolism and Hippo–YAP programs across MASLD severity and MASLD-associated hepatocellular carcinoma.
# MASLD–HCC Transcriptomics

Cross-cohort transcriptomic analysis of lipid metabolism and Hippo–YAP programs across metabolic dysfunction-associated steatotic liver disease (MASLD) severity and MASLD-associated hepatocellular carcinoma (HCC).

> **Project status:** Study design and data preparation  
> **Project type:** Reanalysis of publicly available human transcriptomic datasets  
> **Primary language:** R

## Overview

Metabolic dysfunction-associated steatotic liver disease (MASLD) encompasses a spectrum ranging from hepatic steatosis to metabolic dysfunction-associated steatohepatitis (MASH), fibrosis, and hepatocellular carcinoma (HCC).

Altered lipid metabolism is a central feature of MASLD, while Hippo–YAP signaling is involved in liver injury, regeneration, fibrosis, and tumor development. However, it remains unclear whether lipid-metabolic alterations and YAP-related transcriptional activity show consistent patterns across independent datasets representing different stages of MASLD and MASLD-associated HCC.

This project will reanalyse publicly available human liver transcriptomic datasets to characterize these molecular programs across disease severity and cancer status.

## Research question

**Do lipid-metabolic and Hippo–YAP transcriptional programs change consistently across MASLD severity and MASLD-associated hepatocellular carcinoma in independent human transcriptomic cohorts?**

中文：在不同的人肝脏转录组队列中，脂质代谢与 Hippo–YAP 转录程序是否随 MASLD 严重程度及 MASLD 相关肝细胞癌呈现一致变化？

## Working hypothesis

Increasing MASLD severity and MASLD-associated HCC may be accompanied by coordinated alterations in lipid-metabolic pathways and selected YAP/TEAD-related transcriptional programs.

This is a working hypothesis to be tested rather than an established conclusion.

## Objectives

1. Characterize transcriptomic changes across healthy liver, steatosis, and MASH.
2. Examine associations between pathway activity and histological disease severity.
3. Compare MASLD-associated HCC tumors with adjacent non-tumor tissues.
4. Evaluate whether lipid-metabolic and Hippo–YAP signals are reproducible across independent cohorts.
5. Generate a reproducible R workflow and publication-ready figures.

## Public datasets

| Dataset | Planned use | Study design |
|---|---|---|
| [GSE126848](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126848) | Discovery of expression and pathway changes across MASLD stages | Human liver RNA-seq including healthy, obese, steatosis, and MASH groups |
| [GSE225740](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225740) | Association with histological severity | Human liver biopsies with disease activity and fibrosis information |
| [GSE164441](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164441) | Evaluation in MASLD-associated HCC | Paired tumor and adjacent non-tumor liver tissues |

Dataset inclusion and final sample numbers may be revised after metadata review and quality control.

## Planned analysis

1. Download and curate public datasets and clinical metadata.
2. Conduct sample-level quality control and exploratory analysis.
3. Perform differential gene-expression analysis within each cohort.
4. Quantify lipid metabolism, fibrosis, inflammation, and Hippo–YAP pathway activity.
5. Compare tumors with matched adjacent non-tumor tissues.
6. Evaluate the consistency of gene- and pathway-level changes across cohorts.
7. Conduct sensitivity analyses and document all analytical decisions.

Count-based RNA-seq datasets will primarily be analysed using `DESeq2`. Pathway analyses may use `fgsea`, `GSVA`, `msigdbr`, and related R packages.

Each cohort will first be analysed independently. Expression matrices from different studies will not be directly pooled without a justified and validated integration strategy.

## Pathways of interest

Predefined biological programs include:

- Fatty-acid metabolism
- Adipogenesis
- Cholesterol homeostasis
- Bile-acid metabolism
- Oxidative phosphorylation
- Reactive oxygen species
- Inflammatory signaling
- Extracellular matrix and fibrosis
- TGF-β signaling
- Hippo–YAP/TAZ signaling
- YAP–TEAD target-gene activity
- Protein palmitoylation-related genes

Exploratory analyses will be clearly distinguished from these predefined analyses.

## Expected figures

1. **Study design and cohort overview**  
   Sample composition, clinical groups, and analytical workflow.

2. **Transcriptomic landscape of MASLD severity**  
   PCA, differential-expression summary, and pathway changes across healthy liver, steatosis, and MASH.

3. **Association with histological severity**  
   Relationships between pathway scores, disease activity, and fibrosis stage.

4. **MASLD-associated HCC analysis**  
   Paired comparison of tumor and adjacent tissues, focusing on lipid metabolism and Hippo–YAP programs.

5. **Cross-cohort synthesis**  
   Consistency of gene- and pathway-level changes across independent datasets.

These are planned outputs and may be modified according to data quality and metadata availability.

## Project boundaries

This project is based on secondary analysis of public, primarily cross-sectional transcriptomic datasets.

Therefore:

- The study can identify associations and reproducible transcriptional patterns.
- It cannot demonstrate longitudinal disease progression in the same individuals.
- It cannot establish causal relationships between lipid metabolism, YAP activity, and HCC development.
- Transcript abundance does not directly establish protein activity, YAP nuclear localization, or protein palmitoylation.
- Proposed mechanisms will require experimental or prospective clinical validation.
- Predictive models will not be presented as clinically applicable without external validation.

The project will use terms such as **associated with**, **consistent with**, and **suggests**, rather than making unsupported causal claims.

## Repository structure

<pre>
masld-hcc-transcriptomics/
├── README.md
├── masld-hcc-transcriptomics.Rproj
├── renv.lock
├── analysis_log.md
├── config/
│   └── datasets.yml
├── data/
│   ├── raw/
│   ├── metadata/
│   └── processed/
├── R/
│   ├── download_data.R
│   ├── curate_metadata.R
│   ├── quality_control.R
│   ├── differential_expression.R
│   ├── pathway_analysis.R
│   └── plotting.R
├── analysis/
│   ├── 01_data_preparation.Rmd
│   ├── 02_quality_control.Rmd
│   ├── 03_masld_analysis.Rmd
│   ├── 04_hcc_analysis.Rmd
│   └── 05_cross_cohort_analysis.Rmd
├── results/
│   ├── tables/
│   └── figures/
├── docs/
│   ├── study_protocol.md
│   └── data_dictionary.md
└── references/
</pre>

Large data files and intermediate analysis objects will be excluded from version control. Data-download and preprocessing scripts will be provided whenever permitted by the original repositories.

## Reproducibility

The analysis will be conducted in R. Package versions will be managed using [`renv`](https://rstudio.github.io/renv/).

After cloning the repository, the R environment can be restored with:

<pre>
install.packages("renv")
renv::restore()
</pre>

Exact reproduction instructions will be added after completion of the data-processing workflow.

## Planned deliverables

By the end of the initial one-month project period, this repository aims to contain:

- Curated dataset documentation
- Reproducible R analysis scripts
- A complete analysis log
- Quality-control reports
- Three to five principal scientific figures
- Processed summary tables
- A Chinese research report or manuscript outline
- A scientific poster draft

## Analysis log

Important analytical decisions, failed attempts, parameter changes, and daily progress will be documented in [`analysis_log.md`](analysis_log.md).

## Data and code availability

All datasets used in this project are publicly available through the NCBI Gene Expression Omnibus. This repository will contain analysis code and derived summary results but will not redistribute sensitive, controlled, or unnecessarily large source data.

## License

The code in this repository is released under the [MIT License](LICENSE).

## Disclaimer

This is an exploratory academic bioinformatics project. The results are not intended for clinical diagnosis, treatment selection, or medical decision-making.
