#!/usr/bin/env Rscript
files <- c("workflow", "geo", "import_geo", "associations", "program_analysis", "discovery", "histology", "hcc_paired")
for (file in files) source(file.path("R", paste0(file, ".R")))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript --vanilla analysis/run_study.R VERIFIED_RAW_DIR NEW_OUTPUT_DIR")
raw <- args[1]
out <- args[2]
if (file.exists(out)) stop("Output path already exists: ", out)
for (pkg in c("DESeq2", "limma")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing study dependency: ", pkg)
}
inputs <- read.delim("config/study_sources.tsv", colClasses = "character")
input_paths <- file.path(raw, inputs$file)
if (!all(file.exists(input_paths))) stop("Missing source files: ", paste(inputs$file[!file.exists(input_paths)], collapse = ", "))
programs <- load_study_programs(raw)
import <- tempfile("gse126848-")
import_gse126848(file.path(raw, "GSE126848_Gene_counts_raw.txt.gz"),
                 file.path(raw, "GSE126848_family.soft.gz"), import)
discovery <- list(counts = read_counts(file.path(import, "counts.tsv")),
                  metadata = read_metadata(file.path(import, "samples.tsv")),
                  mapping = read.delim(file.path(import, "sample_mapping.tsv"), colClasses = "character"))
histology <- import_gse225740(raw)
dir.create(out, recursive = TRUE)
result <- study_discovery(discovery, programs, file.path(out, "GSE126848"))
traits <- study_histology(histology, programs, file.path(out, "GSE225740"))
hcc <- study_hcc(raw, programs, file.path(out, "GSE164441"))

a <- result$pathways
b <- traits
c <- hcc
discovery_n <- setNames(result$contrasts$n_case + result$contrasts$n_reference, result$contrasts$contrast)
summary <- rbind(
  data.frame(dataset = "GSE126848", analysis = a$contrast, program = a$program,
    method = "CAMERA_PR_Wald_rank", n = unname(discovery_n[a$contrast]), unit = "samples_in_contrast",
    matched = a$matched, coverage = a$coverage, direction = a$direction, estimate = a$mean_rank_shift,
    ci_low = NA_real_, ci_high = NA_real_, pvalue = a$pvalue, qvalue = a$qvalue_all_discovery, status = a$status),
  data.frame(dataset = "GSE225740", analysis = b$trait, program = b$program, method = "Spearman_rho",
    n = b$n, unit = "patients", matched = b$matched, coverage = b$coverage,
    direction = ifelse(is.na(b$estimate), NA_character_, ifelse(b$estimate > 0, "Up", ifelse(b$estimate < 0, "Down", "Zero"))),
    estimate = b$estimate, ci_low = b$ci_low, ci_high = b$ci_high,
    pvalue = b$pvalue, qvalue = b$qvalue_all_histology, status = b$status),
  data.frame(dataset = "GSE164441", analysis = "tumor_vs_adjacent", program = c$program, method = "CAMERA_PR_moderated_t_rank",
    n = 10L, unit = "pairs", matched = c$matched, coverage = c$coverage, direction = c$direction,
    estimate = c$mean_rank_shift, ci_low = NA_real_, ci_high = NA_real_, pvalue = c$pvalue,
    qvalue = c$qvalue, status = c$status)
)
write_tsv(summary, file.path(out, "cross_cohort_programs.tsv"))
write_tsv(programs$definitions, file.path(out, "program_definitions.tsv"))
inputs$observed_md5 <- unname(tools::md5sum(input_paths))
write_tsv(inputs, file.path(out, "source_manifest.tsv"))
code <- c("analysis/run_study.R", file.path("R", paste0(files, ".R")),
          "config/programs.tsv", "config/discovery_contrasts.tsv", "config/study_sources.tsv",
          "data/metadata/GSE164441_pairs.tsv")
write_tsv(data.frame(file = code, md5 = unname(tools::md5sum(code))), file.path(out, "code_checksums.tsv"))
packages <- installed.packages()
write_tsv(data.frame(package = packages[, "Package"], version = packages[, "Version"], row.names = NULL),
          file.path(out, "installed_packages.tsv"))
writeLines(sub("[[:blank:]]+$", "", capture.output(sessionInfo())), file.path(out, "sessionInfo.txt"))
writeLines("All three cohorts finished", file.path(out, "_SUCCESS"))
print(result$contrasts, row.names = FALSE)
print(summary[, c("dataset", "analysis", "program", "n", "direction", "estimate", "qvalue")], row.names = FALSE)
