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

## 2026-11-17 — Runnable workflow and cohort audit

- Implemented a base-R command-line workflow for strict count/metadata validation, QC, PCA tables, user-GMT program scores and an offline HTML report.
- Added an optional DESeq2 contrast with explicit direction and complete-pair validation. DESeq2 is not installed in the local validation environment; its integration checks are skipped there.
- Added a GEO family SOFT metadata downloader and parser, input/code checksums and session provenance.
- Checked GEO Series records and curated GSE164441 public sample-title pair mappings (20 samples / 10 pairs).
- GSE164441 Cuffdiff summaries require a format/unit audit before any count-based analysis; raw sequencing data are not publicly supplied.
- Added offline tests, a reproducible synthetic example, CI and user documentation. Synthetic expression does not establish biological findings.
- Local GEO download attempts failed (network initialization in sandbox; timeout after network escalation). Parser tests use a small SOFT fixture; live downloading remains unverified.
- Kept the original renv lockfile; it only locks renv, not the optional DESeq2 stack.

Next: audit actual expression columns, annotation and covariates, select documented gene sets, and complete cohort-specific analyses.
