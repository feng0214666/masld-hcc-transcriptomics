# Exercise GMT parsing, covariate typing, trait switches and report outputs together.
source("R/workflow.R")
d <- demo_inputs(); m <- d$metadata; m$patient_id <- NULL
m$NAS <- rep(0:5, each = 2)
m$sex <- rep(c("F", "M"), 6)
root <- tempfile("masld-association-test-"); dir.create(root)
write_tsv(data.frame(gene_id = rownames(d$counts), d$counts, check.names = FALSE), file.path(root, "counts.tsv"))
write_tsv(m, file.path(root, "samples.tsv"))
writeLines(paste(c("SIMULATED_PROGRAM", "synthetic-test-only", rownames(d$counts)[1:12]), collapse = "\t"), file.path(root, "sets.gmt"))
rscript <- file.path(R.home("bin"), paste0("Rscript", if (.Platform$OS.type == "windows") ".exe" else ""))
out <- file.path(root, "output")
args <- c("--vanilla", "analysis/run.R", "--counts", shQuote(file.path(root, "counts.tsv")),
           "--metadata", shQuote(file.path(root, "samples.tsv")), "--gmt", shQuote(file.path(root, "sets.gmt")),
           "--trait", "NAS", "--covariates", "sex", "--dataset", "SYNTHETIC_TEST", "--out", shQuote(out))
result <- system2(rscript, args, stdout = TRUE, stderr = TRUE)
if (!is.null(attr(result, "status"))) stop(paste(result, collapse = "\n"))
a <- read.delim(file.path(out, "trait_associations.tsv"))
stopifnot(nrow(a) == 1, a$method == "adjusted_OLS", a$n == 12, a$status == "ok", is.finite(a$padj))
html <- paste(readLines(file.path(out, "report.html")), collapse = "\n")
stopifnot(grepl("Trait associations", html, fixed = TRUE), grepl("SYNTHETIC_TEST", html, fixed = TRUE))
cat("GMT + adjusted trait association CLI and report integration passed\n")
