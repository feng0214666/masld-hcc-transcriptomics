load_cuffdiff_table <- function(path, expected) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  found <- FALSE
  for (i in seq_len(30)) {
    line <- readLines(con, n = 1, warn = FALSE)
    if (!length(line)) break
    if (startsWith(line, "#test_id\t")) { found <- TRUE; break }
  }
  if (!found) stop("Missing Cuffdiff tabular header: ", basename(path))
  fields <- strsplit(substring(line, 2L), "\t", fixed = TRUE)[[1]]
  required <- c("test_id", "gene_id", expected)
  if (anyDuplicated(fields) || !all(required %in% fields)) stop("Unexpected Cuffdiff columns: ", basename(path))
  types <- ifelse(fields %in% expected, "numeric", ifelse(fields %in% c("test_id", "gene_id"), "character", "NULL"))
  withCallingHandlers(read.delim(con, header = FALSE, col.names = fields, colClasses = types,
    quote = "", comment.char = "", check.names = FALSE, row.names = NULL),
    warning = function(w) stop(basename(path), ": ", conditionMessage(w), call. = FALSE))
}

study_hcc <- function(raw, programs, out) {
  m <- read.delim("data/metadata/GSE164441_pairs.tsv", colClasses = "character")
  if (nrow(m) != 20L || anyDuplicated(m$sample_id) || !all(table(m$patient_id, m$group) == 1L))
    stop("Expected 10 complete GEO pairs.")
  condition <- c(adjacent = "N", tumor = "T")
  m$source_column <- paste0(m$patient_id, condition[m$group], "_FPKM")
  fpkm <- load_cuffdiff_table(file.path(raw, "GSE164441_RNAseq_10Tvs10NT_cuffdiff_FPKM.txt.gz"), m$source_column)
  x <- as.matrix(fpkm[, m$source_column])
  colnames(x) <- m$sample_id
  same_id <- !is.na(fpkm$test_id) & !is.na(fpkm$gene_id) & fpkm$test_id == fpkm$gene_id
  stable <- sub("\\.[0-9]+$", "", fpkm$gene_id)
  annotated <- same_id & grepl("^ENSG[0-9]+$", stable)
  ambiguous <- rep(FALSE, length(stable))
  ambiguous[annotated] <- duplicated(stable[annotated]) | duplicated(stable[annotated], fromLast = TRUE)
  quantified <- rowSums(!is.finite(x) | x < 0) == 0L
  expressive <- rowSums(x >= 1, na.rm = TRUE) >= 5L
  keep <- annotated & !ambiguous & quantified & expressive
  count_path <- file.path(raw, "GSE164441_RNAseq_10Tvs10NT_cuffdiff_count.txt.gz")
  count_con <- gzfile(count_path)
  count_head <- tryCatch(readLines(count_con, n = 10), finally = close(count_con))
  count_header <- grep("^#?test_id\t", count_head, value = TRUE)
  if (length(count_header) != 1) stop("Could not identify the Cuffdiff count header.")
  dir.create(out, recursive = TRUE)
  writeLines(c("Files published by GEO as Cuffdiff estimates.",
    paste("Observed sample-level FPKM columns:", paste(m$source_column, collapse = ", ")),
    paste("Cuffdiff count table header:", count_header),
    "Analysis input: per-sample FPKM. Cuffdiff count file is audited for layout only."),
    file.path(out, "source_format.txt"))
  write_tsv(data.frame(source_rows = nrow(x), ensembl_gene_rows = sum(annotated),
    ambiguous_id_rows = sum(ambiguous), kept = sum(keep), complete_pairs = nrow(m) / 2),
    file.path(out, "cohort_summary.tsv"))
  x <- x[keep, , drop = FALSE]
  rownames(x) <- stable[keep]
  sets <- map_programs(rownames(x), programs, "ensembl")
  m$group <- factor(m$group, levels = c("adjacent", "tumor"))
  m$patient_id <- factor(m$patient_id)
  mm <- model.matrix(~ patient_id + group, m)
  rownames(mm) <- m$sample_id
  co <- match("grouptumor", colnames(mm))
  estimate <- function(offset) {
    fit <- limma::eBayes(limma::lmFit(log2(x + offset), mm), trend = TRUE)
    limma::topTable(fit, coef = co, number = Inf, sort.by = "none", confint = TRUE)
  }
  result <- estimate(0.5)
  sensitivity <- estimate(1)
  res <- data.frame(gene_id = rownames(result), symbol = unname(sets$symbols[rownames(result)]), result)
  write_tsv(res[order(res$adj.P.Val), ], file.path(out, "paired_genes.tsv"))
  p <- preranked_programs(setNames(result$t, rownames(result)), sets)
  alt <- preranked_programs(setNames(sensitivity$t, rownames(sensitivity)), sets)
  p$sensitivity_offset_1_direction <- alt$direction
  p$sensitivity_offset_1_qvalue <- alt$qvalue
  write_tsv(p, file.path(out, "ranked_programs.tsv"))
  write_tsv(m, file.path(out, "sample_mapping.tsv"))
  write_tsv(data.frame(sample_id = rownames(mm), mm, check.names = FALSE), file.path(out, "model_matrix.tsv"))
  p
}
