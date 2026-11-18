study_discovery <- function(input, programs, out) {
  counts <- input$counts
  m <- validate_inputs(counts, input$metadata)
  m$group <- factor(m$group, levels = c("Normal-weight", "Obese", "NAFL", "NASH"))
  m$sex <- factor(m$sex)
  rownames(m) <- m$sample_id
  formula <- ~ sex + group
  mm <- model.matrix(formula, m)
  if (qr(mm)$rank != ncol(mm) || any(table(m$group) < 3)) stop("Invalid discovery study design.")
  rule <- rowSums(counts >= 10) >= 3
  x <- counts[rule, , drop = FALSE]
  storage.mode(x) <- "integer"
  dds <- DESeq2::DESeqDataSetFromMatrix(x, m, formula)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  contrasts <- read.delim("config/discovery_contrasts.tsv", colClasses = "character")
  sets <- map_programs(rownames(x), programs, "ensembl")
  dir.create(out, recursive = TRUE)
  write_tsv(m, file.path(out, "model_samples.tsv"))
  write_tsv(input$mapping, file.path(out, "sample_mapping.tsv"))
  write_tsv(data.frame(sample_id = rownames(mm), mm, check.names = FALSE), file.path(out, "model_matrix.tsv"))
  write_tsv(data.frame(gene_id = rownames(counts), retained = rule), file.path(out, "gene_filter.tsv"))
  writeLines(deparse(formula), file.path(out, "model_formula.txt"))
  summaries <- vector("list", nrow(contrasts))
  pathways <- summaries
  for (i in seq_len(nrow(contrasts))) {
    row <- contrasts[i, ]
    res <- DESeq2::results(dds, contrast = c("group", row$case, row$reference), alpha = 0.05)
    de <- add_wald_intervals(data.frame(gene_id = rownames(res), symbol = unname(sets$symbols[rownames(res)]), as.data.frame(res)))
    de <- de[order(de$padj, na.last = TRUE), ]
    write_tsv(de, file.path(out, paste0(row$contrast, "_genes.tsv")))
    sig <- !is.na(de$padj) & de$padj < 0.05 & is.finite(de$log2FoldChange)
    summaries[[i]] <- data.frame(row, n_case = sum(m$group == row$case), n_reference = sum(m$group == row$reference),
      fitted_genes = nrow(de), adjusted_pvalues = sum(!is.na(de$padj)),
      up_fdr05_abs_lfc1 = sum(sig & de$log2FoldChange >= 1),
      down_fdr05_abs_lfc1 = sum(sig & de$log2FoldChange <= -1))
    statistic <- setNames(de$stat, de$gene_id)
    a <- preranked_programs(statistic, sets, correlation = 0.01)
    b <- preranked_programs(statistic, sets, correlation = 0.05)
    a$pvalue_corr05 <- b$pvalue
    pathways[[i]] <- cbind(contrast = row$contrast, a)
    message(row$contrast, ": ", summaries[[i]]$up_fdr05_abs_lfc1, " up, ", summaries[[i]]$down_fdr05_abs_lfc1, " down")
  }
  pathways <- do.call(rbind, pathways)
  pathways$qvalue_all_discovery <- p.adjust(pathways$pvalue, "BH", n = nrow(pathways))
  pathways$qvalue_all_discovery_corr05 <- p.adjust(pathways$pvalue_corr05, "BH", n = nrow(pathways))
  write_tsv(do.call(rbind, summaries), file.path(out, "contrast_summary.tsv"))
  write_tsv(pathways, file.path(out, "ranked_programs.tsv"))
  list(pathways = pathways, contrasts = do.call(rbind, summaries))
}
