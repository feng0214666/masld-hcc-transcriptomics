# Analysis Log

This file records daily progress, analytical decisions, problems, and changes made during the project.

---

## 2026-07-18 — Day 1: Project setup

### Completed

- Created the public GitHub repository.
- Added the initial project README.
- Defined the research question and project boundaries.
- Selected three candidate GEO datasets: GSE126848, GSE225740, and GSE164441.
- Created the local RStudio project and directory structure.
- Initialized a reproducible R environment using `renv`.

### Research question

Do lipid-metabolic and Hippo–YAP transcriptional programs change consistently across MASLD severity and MASLD-associated hepatocellular carcinoma in independent human transcriptomic cohorts?

### Decisions

- Analyse each cohort independently before cross-cohort synthesis.
- Do not directly merge expression matrices from different studies.
- Distinguish predefined pathway analyses from exploratory analyses.
- Avoid causal language because the datasets are primarily cross-sectional.
- Store only public, deidentified data and derived results.

### Problems encountered

- Command-line access to GitHub was temporarily unavailable during local setup.

### Next step

Audit the sample metadata and downloadable expression files for GSE126848, GSE225740, and GSE164441.
