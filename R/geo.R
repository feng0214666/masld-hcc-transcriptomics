# GEO family SOFT metadata only; retain raw records and never infer case labels silently.
parse_geo_soft <- function(lines) {
  starts <- grep("^\\^SAMPLE = GSM[0-9]+$", lines)
  if (!length(starts)) stop("No GSM sample records found in SOFT; inspect the downloaded response.")
  ends <- c(starts[-1] - 1L, length(lines))
  do.call(rbind, lapply(seq_along(starts), function(i) {
    block <- lines[starts[i]:ends[i]]
    field <- function(key) {
      values <- sub(paste0("^!Sample_", key, " = "), "", grep(paste0("^!Sample_", key, " = "), block, value = TRUE))
      paste(values, collapse = " | ")
    }
    data.frame(sample_id = sub("^\\^SAMPLE = ", "", block[1]), title = field("title"),
               organism = field("organism_ch1"), source = field("source_name_ch1"), description = field("description"),
               characteristics = field("characteristics_ch1"), stringsAsFactors = FALSE)
  }))
}

fetch_geo_metadata <- function(accession, out) {
  if (!grepl("^GSE[0-9]+$", accession)) stop("Expected a GSE accession.")
  if (dir.exists(out)) stop("Metadata output already exists; choose a new directory.")
  digits <- substring(accession, 4)
  prefix <- paste0("GSE", if (nchar(digits) > 3) substr(digits, 1, nchar(digits) - 3) else "", "nnn")
  url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/soft/%s_family.soft.gz", prefix, accession, accession)
  tmp <- tempfile(fileext = ".soft.gz")
  on.exit(unlink(tmp), add = TRUE)
  download.file(url, tmp, mode = "wb", quiet = FALSE)
  con <- gzfile(tmp, "rt"); on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  metadata <- parse_geo_soft(lines)
  if (anyDuplicated(metadata$sample_id)) stop("Duplicate GSM identifiers in SOFT.")
  dir.create(out, recursive = TRUE)
  file.copy(tmp, file.path(out, paste0(accession, "_family.soft.gz")))
  write_tsv(metadata, file.path(out, "metadata_raw.tsv"))
  write_tsv(data.frame(accession = accession, samples_observed = nrow(metadata),
                       url = url, retrieved_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                       md5 = unname(tools::md5sum(tmp))), file.path(out, "provenance.tsv"))
  message(accession, ": ", nrow(metadata), " samples exported; curate labels before analysis.")
}
