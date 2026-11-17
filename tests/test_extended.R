source("R/workflow.R"); source("R/geo.R"); source("R/import_geo.R"); source("R/associations.R")
checks <- 0L
check <- function(ok) { stopifnot(ok); checks <<- checks + 1L }
fails <- function(expr) check(tryCatch({ force(expr); FALSE }, error = function(e) TRUE))
d <- demo_inputs(); m <- d$metadata
m$sex <- rep(c("F", "F", "M", "M"), 3)
m$age <- c(41, 44, 39, 46, 50, 53, 48, 56, 60, 57, 62, 65)
design <- make_design(m, "adjacent", "tumor", covariates = c("sex", "age"))
check(qr(model.matrix(design$formula, design$metadata))$rank == 4)
fails(make_design(m, "adjacent", "tumor", TRUE, "sex")) # sex is patient-constant and aliased
fails(make_design(m, "adjacent", "tumor", covariates = "unknown"))
fails(make_design(m, "adjacent", "tumor", covariates = c("age", "age")))
fails(make_design(m, "adjacent", "tumor", covariates = "age + sex"))
bad <- m; bad$age[1] <- NA; fails(make_design(bad, "adjacent", "tumor", covariates = "age"))
bad <- m; bad$age <- 50; fails(make_design(bad, "adjacent", "tumor", covariates = "age"))
bad <- m; bad$batch <- bad$group; fails(make_design(bad, "adjacent", "tumor", covariates = "batch"))
tmp <- tempfile(); writeLines(c("sample_id\tgroup\tpatient_id\tage", "001\tA\t007\t42", "002\tB\t008\t51"), tmp)
parsed <- read_metadata(tmp)
check(identical(parsed$sample_id, c("001", "002")) && identical(parsed$patient_id, c("007", "008")))
check(is.numeric(parsed$age)); unlink(tmp)
tmp <- tempfile(); writeLines(c("gene_id\tS1\tS2", "001\t2\t3"), tmp)
check(identical(rownames(read_counts(tmp)), "001")); unlink(tmp)

# Controlled data: positive signal, missing trait, repeated subjects and confounding.
m$patient_id <- NULL
m$NAS <- c(0, 1, 1, 2, 2, 3, 4, 4, 5, 6, 7, 7)
scores <- data.frame(program = rep(c("positive", "constant"), each = 12), sample_id = rep(m$sample_id, 2),
                     score = c(m$NAS + c(.1, -.1, .2, 0, -.2, .1, 0, -.1, .2, -.1, 0, .1), rep(1, 12)))
a <- associate_programs(scores, m, "NAS")
check(a$results$estimate[1] > .95 && a$results$pvalue[1] < .001)
check(a$results$status[2] == "constant_variable" && is.na(a$results$padj[2]))
check(abs(a$results$padj[1] - a$results$pvalue[1]) < 1e-10)
a <- associate_programs(scores, m, "NAS", "sex")
check(a$results$estimate[1] > .8 && a$results$ci_low[1] > 0 && a$results$method[1] == "adjusted_OLS")
bad <- m; bad$NAS[1:2] <- NA
a <- associate_programs(scores, bad, "NAS")
check(a$results$n[1] == 10 && sum(!a$exclusions$included) == 2)
bad <- m; bad$NAS[1:10] <- NA
check(associate_programs(scores, bad, "NAS")$results$status[1] == "insufficient_samples")
bad <- m; bad$age <- bad$NAS
check(associate_programs(scores, bad, "NAS", "age")$results$status[1] == "nonestimable_design")
fails(associate_programs(scores, d$metadata, "NAS"))
fails(associate_programs(scores, m, "group"))
fails(associate_programs(scores, m, "NAS", "NAS"))
bad <- m; bad$patient_id <- rep(1:6, each = 2); fails(associate_programs(scores, bad, "NAS"))
bad <- scores; bad$score[1] <- NA; fails(associate_programs(bad, m, "NAS"))

# The adapter maps by documented numeric IDs, not matrix position or title order.
soft <- data.frame(sample_id = paste0("GSM", 1:57), description = as.character(1:57),
                    organism = "Homo sapiens", title = paste0("NAFL_", 1:57),
                    characteristics = "gender: Male | tissue: Liver")
mapped <- map_gse126848(sprintf("%04d", 57:1), soft)
check(mapped$sample_id[1] == "GSM57" && mapped$source_description[1] == "57")
bad <- soft; bad$description[1] <- "2"; fails(map_gse126848(sprintf("%04d", 1:57), bad))
bad <- soft; bad$title[1] <- "unexpected"; fails(map_gse126848(sprintf("%04d", 1:57), bad))
fails(map_gse126848(as.character(1:56), soft))
cat(checks, "extended checks passed\n")
