# Explicit adapter for GSE126848, rather than a heuristic for arbitrary studies.
map_gse126848 <- function(expression_columns, soft_metadata) {
  if (nrow(soft_metadata) != 57 || any(soft_metadata$organism != "Homo sapiens"))
    stop("GSE126848 adapter expects 57 human samples; inspect a changed source record.")
  canonical <- function(ids) {
    if (anyNA(ids) || any(!grepl("^[0-9]+$", ids))) stop("Expected numeric source sample identifiers.")
    result <- sub("^0+", "", ids)
    if (any(!nzchar(result)) || anyDuplicated(result)) stop("Ambiguous source IDs after leading-zero normalization.")
    result
  }
  observed <- canonical(expression_columns)
  declared <- canonical(soft_metadata$description)
  if (!setequal(observed, declared)) stop("Expression columns and SOFT descriptions do not match exactly.")
  m <- soft_metadata[match(observed, declared), , drop = FALSE]
  group <- sub("_[0-9]+$", "", m$title)
  if (any(!group %in% c("NAFL", "NASH", "Normal-weight", "Obese"))) stop("Unexpected source group labels.")
  sex <- trimws(sub("^.*gender: ([^|]*).*$", "\\1", m$characteristics))
  if (any(!sex %in% c("Male", "Female"))) stop("Missing/unexpected sex labels; inspect source metadata.")
  data.frame(expression_column = expression_columns, sample_id = m$sample_id, source_title = m$title,
             source_description = m$description, group = group, sex = sex,
             mapping_rule = "numeric_source_id_leading_zeros_removed_for_matching_only", stringsAsFactors = FALSE)
}

import_gse126848 <- function(counts_path, soft_path, out) {
  if (file.exists(out)) stop("Import directory already exists.")
  ccon <- gzfile(counts_path, "rt"); on.exit(close(ccon), add = TRUE)
  scon <- gzfile(soft_path, "rt"); on.exit(close(scon), add = TRUE)
  strict <- function(w) stop("Source read warning: ", conditionMessage(w))
  d <- withCallingHandlers(read.delim(ccon, check.names = FALSE, colClasses = c(key = "character")), warning = strict)
  if (names(d)[1] != "key" || !all(vapply(d[-1], is.numeric, logical(1)))) stop("Unexpected GSE126848 expression layout.")
  metadata_raw <- parse_geo_soft(withCallingHandlers(readLines(scon), warning = strict))
  mapping <- map_gse126848(names(d)[-1], metadata_raw)
  counts <- as.matrix(d[-1]); rownames(counts) <- d[[1]]; colnames(counts) <- mapping$sample_id
  metadata <- mapping[, c("sample_id", "group", "sex", "source_title")]
  metadata <- validate_inputs(counts, metadata)
  dir.create(out, recursive = TRUE)
  write_tsv(data.frame(gene_id = rownames(counts), counts, check.names = FALSE), file.path(out, "counts.tsv"))
  write_tsv(metadata, file.path(out, "samples.tsv"))
  write_tsv(mapping, file.path(out, "sample_mapping.tsv"))
  write_tsv(data.frame(source_file = c(counts_path, soft_path), md5 = unname(tools::md5sum(c(counts_path, soft_path)))),
            file.path(out, "source_checksums.tsv"))
  message("Imported GSE126848: ", nrow(counts), " genes, ", ncol(counts), " samples.")
}
