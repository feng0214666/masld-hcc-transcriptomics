#!/usr/bin/env Rscript
# Usage: Rscript --vanilla analysis/run.R --demo --out results/demo
source("R/workflow.R")
source("R/associations.R")
args <- commandArgs(trailingOnly = TRUE)
flags <- c("--demo", "--deseq2", "--paired")
options <- c("--out", "--counts", "--metadata", "--gmt", "--reference", "--case", "--covariates", "--trait")
opt <- list(); i <- 1L
while (i <= length(args)) {
  key <- args[i]
  if (key %in% flags) { opt[[key]] <- TRUE; i <- i + 1L
  } else if (key %in% options && i < length(args) && !startsWith(args[i + 1L], "--")) {
    opt[[key]] <- args[i + 1L]; i <- i + 2L
  } else stop("Unknown option or missing value: ", key)
}
demo <- isTRUE(opt[["--demo"]])
out <- opt[["--out"]]
if (is.null(out)) stop("Specify --out with a NEW output directory.")
if (file.exists(out)) stop("Output directory already exists; choose a new --out to preserve run provenance.")
covariates <- if (is.null(opt[["--covariates"]])) character() else strsplit(opt[["--covariates"]], ",", fixed = TRUE)[[1]]
if (!is.null(opt[["--covariates"]]) && !length(covariates)) stop("Empty covariate list.")
inputs <- character()
if (demo) {
  if (any(c("--counts", "--metadata", "--gmt") %in% names(opt))) stop("--demo cannot be mixed with external inputs.")
  d <- demo_inputs(); counts <- d$counts; metadata <- d$metadata
  sets <- list(DEMO_PROGRAM = rownames(counts)[1:12], DEMO_BACKGROUND = rownames(counts)[40:60])
} else {
  if (is.null(opt[["--counts"]]) || is.null(opt[["--metadata"]])) stop("Supply --counts and --metadata, or --demo.")
  counts <- read_counts(opt[["--counts"]])
  metadata <- read_metadata(opt[["--metadata"]])
  inputs <- c(opt[["--counts"]], opt[["--metadata"]])
  sets <- NULL
  if (!is.null(opt[["--gmt"]])) { sets <- read_gmt(opt[["--gmt"]]); inputs <- c(inputs, opt[["--gmt"]]) }
}
metadata <- validate_inputs(counts, metadata)
qc <- qc_counts(counts)
de <- NULL; design <- NULL
if (length(covariates) && !isTRUE(opt[["--deseq2"]]) && is.null(opt[["--trait"]])) stop("--covariates requires --deseq2 or --trait.")
if (isTRUE(opt[["--paired"]]) && !isTRUE(opt[["--deseq2"]])) stop("--paired requires --deseq2.")
if (isTRUE(opt[["--deseq2"]])) {
  if (is.null(opt[["--reference"]]) || is.null(opt[["--case"]])) stop("DESeq2 requires --reference and --case.")
  design <- make_design(metadata, opt[["--reference"]], opt[["--case"]], isTRUE(opt[["--paired"]]), covariates)
  de <- run_deseq2(counts, metadata, opt[["--reference"]], opt[["--case"]], isTRUE(opt[["--paired"]]), covariates)
}
programs <- if (!is.null(sets)) score_sets(qc$logcpm, sets) else NULL
association <- NULL
if (!is.null(opt[["--trait"]])) {
  if (is.null(programs) || !nrow(programs$scores)) stop("--trait requires scored GMT programs.")
  association <- associate_programs(programs$scores, metadata, opt[["--trait"]], covariates)
}
dir.create(out, recursive = TRUE)
if (demo) {
  write_tsv(data.frame(gene_id = rownames(counts), counts, check.names = FALSE), file.path(out, "demo_counts.tsv"))
  write_tsv(metadata, file.path(out, "demo_metadata.tsv"))
  inputs <- file.path(out, c("demo_counts.tsv", "demo_metadata.tsv"))
}
write_tsv(qc$samples, file.path(out, "sample_qc.tsv"))
write_tsv(metadata, file.path(out, "metadata_aligned.tsv"))
write_tsv(qc$pca, file.path(out, "pca_coordinates.tsv"))
write_tsv(qc$variance, file.path(out, "pca_variance.tsv"))
write_tsv(data.frame(gene_id = rownames(counts), retained = qc$keep), file.path(out, "gene_filter.tsv"))
cors <- cor(qc$logcpm)
write_tsv(data.frame(sample_id = rownames(cors), cors, check.names = FALSE), file.path(out, "sample_correlations.tsv"))
if (!is.null(programs)) {
  write_tsv(programs$coverage, file.path(out, "program_coverage.tsv"))
  write_tsv(programs$scores, file.path(out, "program_scores.tsv"))
}
if (!is.null(de)) {
  write_tsv(de, file.path(out, "differential_expression.tsv"))
  write_tsv(design$metadata, file.path(out, "model_samples.tsv"))
  mm <- model.matrix(design$formula, design$metadata)
  write_tsv(data.frame(sample_id = rownames(mm), mm, check.names = FALSE), file.path(out, "model_matrix.tsv"))
  writeLines(deparse(design$formula), file.path(out, "model_formula.txt"))
}
if (!is.null(association)) {
  write_tsv(association$results, file.path(out, "trait_associations.tsv"))
  write_tsv(association$exclusions, file.path(out, "trait_sample_inclusion.tsv"))
}
write_tsv(data.frame(input = inputs, md5 = unname(tools::md5sum(inputs))), file.path(out, "input_checksums.tsv"))
code <- c("analysis/run.R", "R/workflow.R", "R/associations.R")
write_tsv(data.frame(file = code, md5 = unname(tools::md5sum(code))), file.path(out, "code_checksums.tsv"))
write_tsv(data.frame(parameter = c("mode", "command", "filter", "score_method", "R"),
                    value = c(if (demo) "synthetic" else "user_counts", paste(args, collapse = " "),
                              "count>=10 in >=3 samples", "within-cohort mean gene z-score; >=3 genes, >=50% coverage", R.version.string)),
          file.path(out, "run_info.tsv"))
writeLines(sub("[[:blank:]]+$", "", capture.output(sessionInfo())), file.path(out, "sessionInfo.txt"))
write_report(out, qc, metadata, demo, !is.null(de), !is.null(programs))
message("Complete: ", file.path(out, "report.html"))
