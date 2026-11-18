read_study_gzip <- function(path, reader, ...) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  withCallingHandlers(reader(con, ...), warning = function(w) {
    stop(basename(path), ": ", conditionMessage(w), call. = FALSE)
  })
}

import_gse225740 <- function(raw) {
  counts_path <- file.path(raw, "GSE225740_NAFLD_RNA_counts.csv.gz")
  traits_path <- file.path(raw, "GSE225740_NAFLD_phenotypes.csv.gz")
  soft_path <- file.path(raw, "GSE225740_family.soft.gz")
  counts <- read_study_gzip(counts_path, read.csv, check.names = FALSE, strip.white = TRUE,
                            colClasses = c(EntrezGene = "character"))
  traits <- read_study_gzip(traits_path, read.csv, check.names = FALSE, strip.white = TRUE)
  soft <- parse_geo_soft(read_study_gzip(soft_path, readLines))
  required <- c("pt. ID", "nas", "fibrosis", "steatohep")
  if (names(counts)[1] != "EntrezGene" || !identical(trimws(names(traits)), required) ||
      !all(vapply(counts[-1], is.numeric, logical(1)))) stop("Unexpected GSE225740 files.")
  names(traits) <- required
  ids <- trimws(names(counts)[-1])
  for (observed in list(ids, as.character(traits[["pt. ID"]]), soft$description)) {
    if (length(observed) != 93L || anyDuplicated(observed) ||
        !setequal(observed, as.character(seq_len(93)))) stop("GSE225740 patient identifiers are not bijective.")
  }
  soft <- soft[match(ids, soft$description), ]
  traits <- traits[match(ids, traits[["pt. ID"]]), ]
  if (!all(soft$organism == "Homo sapiens" & soft$source == soft$source[1]) ||
      !all(grepl("^GSM[0-9]+$", soft$sample_id))) stop("Unexpected GSE225740 sample labels.")
  read_trait <- function(name, domain) {
    value <- traits[[name]]
    if (!is.numeric(value) || any(!is.na(value) & !value %in% domain))
      stop("Unexpected phenotype levels for ", name)
    fields <- strsplit(soft$characteristics, " | ", fixed = TRUE)
    geo <- vapply(fields, function(x) {
      prefix <- paste0(name, ": ")
      matched <- x[startsWith(x, prefix)]
      if (length(matched) != 1L) stop("Missing or ambiguous SOFT field: ", name)
      substring(matched, nchar(prefix) + 1L)
    }, character(1))
    if (!all(geo == ifelse(is.na(value), "NA", as.character(value))))
      stop("Disagreement between phenotype CSV and SOFT for ", name)
    value
  }
  metadata <- data.frame(sample_id = soft$sample_id, group = "liver_biopsy", patient_id = ids,
    NAS = read_trait("nas", 0:8), fibrosis = read_trait("fibrosis", 0:4),
    steatohepatitis_code = read_trait("steatohep", 0:2))
  x <- as.matrix(counts[-1])
  dimnames(x) <- list(counts[[1]], metadata$sample_id)
  list(counts = x, metadata = validate_inputs(x, metadata),
       mapping = data.frame(expression_column = ids, sample_id = metadata$sample_id, source_title = soft$title))
}

expression_duplicates <- function(x) {
  # Group on cheap signatures first, then compare full columns to rule out signature collisions.
  by_sum <- split(seq_len(ncol(x)), colSums(x))
  match_group <- rep(NA_character_, ncol(x))
  for (indices in by_sum[lengths(by_sum) > 1]) {
    for (j in seq_along(indices)) {
      i <- indices[j]
      if (!is.na(match_group[i])) next
      matched <- indices[vapply(indices, function(k) identical(x[, i], x[, k]), logical(1))]
      if (length(matched) > 1) match_group[matched] <- colnames(x)[matched[1]]
    }
  }
  data.frame(sample_id = colnames(x), duplicate_group = match_group, excluded = !is.na(match_group))
}

study_histology <- function(study, programs, out) {
  dir.create(out, recursive = TRUE)
  duplicates <- expression_duplicates(study$counts)
  retain <- !duplicates$excluded
  counts <- study$counts[, retain, drop = FALSE]
  metadata <- study$metadata[retain, , drop = FALSE]
  qc <- qc_counts(counts)
  mapping <- map_programs(rownames(qc$logcpm), programs, "entrez")
  scored <- score_study_programs(qc$logcpm, mapping)
  set.seed(225740)
  results <- lapply(c("NAS", "fibrosis"), function(trait) {
    association <- associate_programs(scored$scores, metadata, trait)
    result <- merge(scored$coverage, association$results, by = "program", all.x = TRUE, sort = FALSE)
    result$trait <- trait
    result$status[!result$scored] <- "insufficient_coverage"
    for (i in which(result$status == "ok")) {
      score <- scored$scores[scored$scores$program == result$program[i], ]
      y <- metadata[[trait]][match(score$sample_id, metadata$sample_id)]
      keep <- is.finite(y)
      x <- score$score[keep]
      y <- y[keep]
      boot <- replicate(2000, {
        ids <- sample.int(length(y), replace = TRUE)
        if (length(unique(y[ids])) < 2L || length(unique(x[ids])) < 2L) return(NA_real_)
        cor(x[ids], y[ids], method = "spearman")
      })
      result[i, c("ci_low", "ci_high")] <- as.list(quantile(boot, c(0.025, 0.975), na.rm = TRUE))
    }
    coverage <- data.frame(trait = trait, sample_id = study$metadata$sample_id,
      included = retain & is.finite(study$metadata[[trait]]),
      reason = ifelse(!retain, "duplicate_expression_profile",
        ifelse(is.finite(study$metadata[[trait]]), "included", "missing_trait")))
    write_tsv(coverage, file.path(out, paste0(trait, "_sample_inclusion.tsv")))
    result
  })
  results <- do.call(rbind, results)
  results$qvalue_all_histology <- p.adjust(results$pvalue, "BH", n = nrow(results))
  write_tsv(results, file.path(out, "trait_programs.tsv"))
  write_tsv(study$metadata, file.path(out, "sample_metadata.tsv"))
  write_tsv(study$mapping, file.path(out, "sample_mapping.tsv"))
  write_tsv(duplicates, file.path(out, "duplicate_profiles.tsv"))
  write_tsv(scored$scores, file.path(out, "program_scores.tsv"))
  write_tsv(scored$coverage, file.path(out, "program_coverage.tsv"))
  write_tsv(data.frame(input_genes = nrow(study$counts), retained_genes = sum(qc$keep),
    input_samples = nrow(study$metadata), duplicate_samples = sum(!retain),
    nas_n = sum(is.finite(metadata$NAS)), fibrosis_n = sum(is.finite(metadata$fibrosis))),
    file.path(out, "cohort_summary.tsv"))
  results
}
