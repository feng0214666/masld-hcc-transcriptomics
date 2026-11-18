# Validate artifacts from the real run independently of report generation.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- args[1]
table_at <- function(path) read.delim(file.path(out, path), check.names = FALSE)
discovery <- table_at("GSE126848/contrast_summary.tsv")
stopifnot(identical(discovery$contrast, c("nash_vs_obese", "nash_vs_nafl", "nash_vs_normal")),
          identical(discovery$n_reference, c(12L, 15L, 14L)), all(discovery$n_case == 16),
          nrow(table_at("GSE126848/model_samples.tsv")) == 57L)
hist <- table_at("GSE225740/trait_programs.tsv")
stopifnot(all(hist$n[hist$trait == "NAS"] == 93), all(hist$n[hist$trait == "fibrosis"] == 87),
          all(hist$ci_low <= hist$estimate & hist$estimate <= hist$ci_high))
pairs <- table_at("GSE164441/sample_mapping.tsv")
stopifnot(nrow(pairs) == 20, all(table(pairs$patient_id, pairs$group) == 1),
          !anyDuplicated(pairs$source_column))
all <- table_at("cross_cohort_programs.tsv")
stopifnot(nrow(all) == 24L, !anyDuplicated(all[c("dataset", "analysis", "program")]),
          all(all$coverage >= .5), all(is.finite(all$pvalue)), all(is.finite(all$qvalue)))
for (dataset in unique(all$dataset)) {
  rows <- all[all$dataset == dataset, ]
  stopifnot(isTRUE(all.equal(rows$qvalue, p.adjust(rows$pvalue, "BH"), tolerance = 1e-12)))
}
stopifnot(file.exists(file.path(out, "_SUCCESS")))
cat("Real cohort sample mappings, coverage, interval and FDR families passed\n")
