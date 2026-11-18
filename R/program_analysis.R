load_study_programs <- function(raw) {
  definitions <- read.delim("config/programs.tsv", colClasses = "character")
  hallmark <- read_gmt(file.path(raw, "h.all.v2025.1.Hs.symbols.gmt"))
  oncogenic <- read_gmt(file.path(raw, "c6.all.v2025.1.Hs.symbols.gmt"))
  by_collection <- list(H = hallmark, C6 = oncogenic)
  sets <- Map(function(id, collection) {
    genes <- by_collection[[collection]][[id]]
    if (is.null(genes)) stop("Required MSigDB set is missing: ", id)
    genes
  }, definitions$set_id, definitions$collection)
  names(sets) <- definitions$program
  definitions$release <- "2025.1.Hs"
  definitions$source_genes <- lengths(sets)
  load_chip <- function(file) {
    chip <- read.delim(file.path(raw, file),
                      check.names = FALSE, colClasses = "character", quote = "")
    required <- c("Probe Set ID", "Gene Symbol")
    if (!all(required %in% names(chip))) stop("Unexpected MSigDB CHIP columns.")
    pairs <- unique(chip[required])
    names(pairs) <- c("id", "symbol")
    pairs <- pairs[!is.na(pairs$symbol) & nzchar(pairs$symbol) &
                   !is.na(pairs$id) & nzchar(pairs$id), ]
    ambiguous <- duplicated(pairs$id) | duplicated(pairs$id, fromLast = TRUE)
    setNames(pairs$symbol[!ambiguous], pairs$id[!ambiguous])
  }
  list(sets = sets, definitions = definitions,
       ensembl = load_chip("Human_Ensembl_Gene_ID_MSigDB.v2025.1.Hs.chip"),
       entrez = load_chip("Human_NCBI_Gene_ID_MSigDB.v2025.1.Hs.chip"))
}

map_programs <- function(gene_ids, programs, id_type) {
  symbols <- switch(id_type,
    ensembl = unname(programs$ensembl[sub("\\.[0-9]+$", "", gene_ids)]),
    entrez = unname(programs$entrez[gene_ids]),
    symbol = gene_ids,
    stop("Unknown expression gene identifier type: ", id_type))
  # Exclude ambiguous many-to-one probe or Ensembl mappings symmetrically.
  ambiguous <- !is.na(symbols) & (duplicated(symbols) | duplicated(symbols, fromLast = TRUE))
  symbols[ambiguous] <- NA_character_
  sets <- lapply(programs$sets, function(genes) gene_ids[!is.na(symbols) & symbols %in% genes])
  list(sets = sets, symbols = setNames(symbols, gene_ids), expected = lengths(programs$sets))
}

preranked_programs <- function(statistic, mapping, correlation = 0.01) {
  statistic <- statistic[is.finite(statistic)]
  if (is.null(names(statistic)) || anyDuplicated(names(statistic)) || length(statistic) < 20)
    stop("Expected at least 20 uniquely named finite gene statistics.")
  sets <- lapply(mapping$sets, intersect, y = names(statistic))
  observed <- lengths(sets)
  eligible <- observed >= 10 & observed / mapping$expected >= 0.5 & observed < length(statistic)
  result <- data.frame(program = names(sets), source_genes = unname(mapping$expected), matched = unname(observed),
    coverage = unname(observed / mapping$expected), ranked_genes = length(statistic),
    direction = NA_character_, mean_rank_shift = NA_real_, pvalue = NA_real_, qvalue = NA_real_,
    status = ifelse(eligible, "tested", "insufficient_coverage"), row.names = NULL)
  if (!any(eligible)) return(result)
  test <- limma::cameraPR(statistic, sets[eligible], use.ranks = TRUE, inter.gene.cor = correlation, sort = FALSE)
  idx <- match(result$program[eligible], rownames(test))
  result$direction[eligible] <- as.character(test$Direction[idx])
  result$pvalue[eligible] <- test$PValue[idx]
  result$qvalue[eligible] <- p.adjust(result$pvalue[eligible], "BH", n = nrow(result))
  ranks <- setNames((rank(statistic) - 0.5) / length(statistic), names(statistic))
  result$mean_rank_shift[eligible] <- vapply(sets[eligible], function(g) mean(ranks[g]) - 0.5, numeric(1))
  result
}

score_study_programs <- function(log_expression, mapping) {
  eligible <- apply(log_expression, 1, sd) > 0
  available <- lapply(mapping$sets, intersect, y = rownames(log_expression)[eligible])
  coverage <- data.frame(program = names(available), source_genes = unname(mapping$expected),
                         matched = unname(lengths(available)), row.names = NULL)
  coverage$coverage <- coverage$matched / coverage$source_genes
  coverage$scored <- coverage$matched >= 10 & coverage$coverage >= 0.5
  scored <- available[coverage$scored]
  if (!length(scored)) stop("No designated program has enough observed genes to score.")
  scores <- score_sets(log_expression, scored)$scores
  list(scores = scores, coverage = coverage)
}

add_wald_intervals <- function(result, estimate = "log2FoldChange", se = "lfcSE") {
  margin <- qnorm(0.975) * result[[se]]
  result$ci_low <- result[[estimate]] - margin
  result$ci_high <- result[[estimate]] + margin
  result
}
