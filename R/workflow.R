# Small, dependency-free building blocks; run scripts from the repository root.
read_counts <- function(path) {
  d <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
  if (ncol(d) < 3L || names(d)[1] != "gene_id")
    stop("Counts TSV must start with gene_id and contain at least two samples.")
  ids <- d[[1]]
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    stop("gene_id must be nonempty and unique; resolve annotation explicitly.")
  if (anyDuplicated(names(d)) || any(!nzchar(names(d)))) stop("Duplicate/empty count columns.")
  for (nm in names(d)[-1]) d[[nm]] <- type.convert(d[[nm]], as.is = TRUE)
  if (!all(vapply(d[-1], is.numeric, logical(1)))) stop("Counts must be numeric.")
  x <- as.matrix(d[-1]); rownames(x) <- ids
  x
}

validate_inputs <- function(counts, metadata) {
  if (!is.matrix(counts) || !is.numeric(counts) || nrow(counts) < 2 || ncol(counts) < 3)
    stop("Expected a numeric matrix with at least two genes and three samples.")
  if (any(!is.finite(counts)) || any(counts < 0)) stop("Counts must be finite and nonnegative.")
  if (any(counts != round(counts)) || any(counts > .Machine$integer.max))
    stop("Expected integer counts; do not round TPM, FPKM or estimated abundances.")
  for (ids in list(rownames(counts), colnames(counts))) {
    if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
      stop("Count row/column IDs must be unique and nonempty.")
  }
  if (!all(c("sample_id", "group") %in% names(metadata))) stop("Metadata needs sample_id and group.")
  if (anyNA(metadata$sample_id) || anyNA(metadata$group) ||
      any(!nzchar(metadata$sample_id)) || any(!nzchar(metadata$group)) ||
      anyDuplicated(metadata$sample_id)) stop("Invalid metadata IDs/groups.")
  if (!setequal(colnames(counts), metadata$sample_id)) stop("Counts and metadata sample IDs must match exactly.")
  if (any(colSums(counts) == 0)) stop("Zero-library samples must be investigated before analysis.")
  metadata[match(colnames(counts), metadata$sample_id), , drop = FALSE]
}

qc_counts <- function(counts, min_count = 10L, min_samples = 3L) {
  keep <- rowSums(counts >= min_count) >= min_samples
  if (sum(keep) < 2L) stop("Fewer than two genes pass the expression filter.")
  # Library sizes use ALL supplied genes; this is exploratory log-CPM, not DESeq2 normalization.
  logcpm <- log2(sweep(counts[keep, , drop = FALSE], 2, colSums(counts), "/") * 1e6 + 1)
  variable <- apply(logcpm, 1, sd) > 0
  if (sum(variable) < 2L) stop("Fewer than two variable genes for PCA.")
  pca <- prcomp(t(logcpm[variable, , drop = FALSE]), center = TRUE, scale. = FALSE)
  list(logcpm = logcpm, keep = keep,
       samples = data.frame(sample_id = colnames(counts), library_size = colSums(counts),
                            detected_genes = colSums(counts > 0), row.names = NULL),
       pca = data.frame(sample_id = rownames(pca$x), pca$x[, 1:2], row.names = NULL),
       variance = data.frame(component = paste0("PC", seq_along(pca$sdev)),
                             variance_fraction = pca$sdev^2 / sum(pca$sdev^2)))
}

read_gmt <- function(path) {
  rows <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  if (!length(rows) || any(lengths(rows) < 3)) stop("GMT needs name, source and gene IDs.")
  ids <- vapply(rows, `[`, character(1), 1)
  if (any(!nzchar(ids)) || anyDuplicated(ids)) stop("GMT set names must be unique.")
  sets <- lapply(rows, function(x) unique(x[-c(1, 2)]))
  names(sets) <- ids
  if (any(vapply(sets, function(x) any(!nzchar(x)), logical(1)))) stop("Empty GMT gene IDs.")
  sets
}

score_sets <- function(logcpm, sets, min_genes = 3L, min_coverage = 0.5) {
  sds <- apply(logcpm, 1, sd)
  z <- sweep(sweep(logcpm[sds > 0, , drop = FALSE], 1,
                   rowMeans(logcpm)[sds > 0], "-"), 1, sds[sds > 0], "/")
  coverage <- do.call(rbind, lapply(names(sets), function(nm) {
    genes <- intersect(sets[[nm]], rownames(z))
    data.frame(program = nm, requested = length(sets[[nm]]), matched = length(genes),
               coverage = length(genes) / length(sets[[nm]]))
  }))
  coverage$scored <- coverage$matched >= min_genes & coverage$coverage >= min_coverage
  scores <- do.call(rbind, lapply(coverage$program[coverage$scored], function(nm) {
    genes <- intersect(sets[[nm]], rownames(z))
    data.frame(program = nm, sample_id = colnames(z), score = colMeans(z[genes, , drop = FALSE]))
  }))
  if (is.null(scores)) scores <- data.frame(program = character(), sample_id = character(), score = numeric())
  list(coverage = coverage, scores = scores)
}

read_metadata <- function(path) {
  m <- read.delim(path, check.names = FALSE, colClasses = "character")
  if (anyDuplicated(names(m))) stop("Duplicate metadata columns.")
  for (nm in setdiff(names(m), c("sample_id", "group", "patient_id")))
    m[[nm]] <- type.convert(m[[nm]], as.is = TRUE)
  m
}

validate_covariates <- function(metadata, covariates, reserved = c("sample_id", "group", "patient_id")) {
  if (anyDuplicated(covariates) || any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", covariates)) ||
      any(covariates %in% reserved)) stop("Invalid, duplicate or reserved covariate name.")
  if (!all(covariates %in% names(metadata))) stop("Missing covariate columns: ", paste(setdiff(covariates, names(metadata)), collapse = ", "))
  for (nm in covariates) {
    v <- metadata[[nm]]
    if (anyNA(v) || any(!nzchar(as.character(v))) || (is.numeric(v) && any(!is.finite(v))))
      stop("Missing/nonfinite covariate: ", nm, "; explicitly curate samples before fitting.")
    if (length(unique(v)) < 2L) stop("Constant covariate: ", nm)
  }
  invisible(TRUE)
}

make_design <- function(metadata, reference, case, paired = FALSE, covariates = character()) {
  if (identical(reference, case)) stop("Reference and case must differ.")
  m <- metadata[metadata$group %in% c(reference, case), , drop = FALSE]
  if (any(table(factor(m$group, levels = c(reference, case))) < 3))
    stop("At least three samples in each contrast group are required.")
  m$group <- factor(m$group, levels = c(reference, case))
  validate_covariates(m, covariates)
  design <- reformulate(c(covariates, "group"))
  if (paired) {
    if (!"patient_id" %in% names(m) || anyNA(m$patient_id) || any(!nzchar(m$patient_id)))
      stop("Paired analysis requires complete patient_id.")
    m$patient_id <- factor(m$patient_id)
    if (any(table(m$patient_id, m$group) != 1L)) stop("Every patient needs exactly one sample per group.")
    design <- reformulate(c("patient_id", covariates, "group"))
  }
  mm <- model.matrix(design, m)
  if (qr(mm)$rank != ncol(mm) || nrow(mm) <= ncol(mm)) stop("Design is rank deficient or has no residual degrees of freedom.")
  rownames(m) <- m$sample_id
  list(metadata = m, formula = design)
}

run_deseq2 <- function(counts, metadata, reference, case, paired = FALSE, covariates = character()) {
  metadata <- validate_inputs(counts, metadata)
  design <- make_design(metadata, reference, case, paired, covariates)
  if (!requireNamespace("DESeq2", quietly = TRUE))
    stop("Install DESeq2 with BiocManager::install('DESeq2') before using --deseq2.")
  x <- counts[, design$metadata$sample_id, drop = FALSE]
  keep <- rowSums(x >= 10) >= 3
  if (sum(keep) < 2) stop("Too few genes for DESeq2.")
  x <- x[keep, , drop = FALSE]; storage.mode(x) <- "integer"
  dds <- DESeq2::DESeqDataSetFromMatrix(x, design$metadata, design$formula)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  result <- as.data.frame(DESeq2::results(dds, contrast = c("group", case, reference), alpha = 0.05))
  result <- data.frame(gene_id = rownames(result), result, row.names = NULL)
  result[order(result$padj, na.last = TRUE), ]
}

html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub('"', "&quot;", x, fixed = TRUE)
}

html_table <- function(d) {
  paste0("<table><thead><tr>", paste0("<th>", html_escape(names(d)), "</th>", collapse = ""),
         "</tr></thead><tbody>", paste(apply(d, 1, function(row)
           paste0("<tr>", paste0("<td>", html_escape(row), "</td>", collapse = ""), "</tr>")), collapse = ""),
         "</tbody></table>")
}

write_report <- function(out, qc, metadata, synthetic, has_de, has_sets, dataset = NULL) {
  label <- if (synthetic) "SYNTHETIC DEMO / 模拟数据" else "USER-SUPPLIED COUNTS / 用户提供计数"
  composition <- as.data.frame(table(metadata$group)); names(composition) <- c("group", "samples")
  body <- c('<!doctype html><html lang="en"><meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width,initial-scale=1">',
    '<title>MASLD–HCC | Analysis report</title><style>',
    'body{font:16px/1.65 system-ui,sans-serif;color:#17333d;background:#f3f7f8;max-width:1120px;margin:40px auto;padding:24px}',
    'header{background:#143b47;color:white;border-radius:18px;padding:32px}h1{margin:4px 0}small{letter-spacing:.1em}',
    '.tag{display:inline-block;padding:6px 14px;background:#d8f2e7;color:#155341;border-radius:30px}',
    'section{background:white;padding:24px;margin:20px 0;border-radius:14px;overflow:auto}',
    'table{border-collapse:collapse;width:100%;font-size:14px}td,th{padding:10px;text-align:left;border-bottom:1px solid #dae6e8}',
    'a{color:#087e83}code{background:#edf4f5;padding:3px}footer{color:#546b73}</style>',
    '<header><small>MASLD–HCC TRANSCRIPTOMICS</small><h1>Cohort analysis report</h1>',
    paste0('<span class="tag">', label, '</span>'),
    if (!is.null(dataset)) paste0('<p>Dataset: ', html_escape(dataset), '</p>'),
    '<p>Input validation → library QC → exploratory PCA → optional program scores / DESeq2</p></header>',
    '<section><h2>Run overview</h2>',
    paste0('<p>', nrow(metadata), ' samples · ', length(qc$keep), ' input genes · ', sum(qc$keep), ' retained genes.</p>'),
    if (synthetic) '<p>All expression values and gene programs in this run are simulated. They establish software behavior and contain no evidence about MASLD or HCC.</p>' else
      '<p>Input provenance and suitability must be reviewed before biological interpretation. QC does not automatically exclude samples.</p>',
    '<p>Filter: count ≥10 in ≥3 samples. PCA uses centered log2(CPM+1), with library sizes computed before filtering.</p></section>',
    '<section><h2>Sample composition</h2>', html_table(composition), '</section>',
    '<section><h2>Exploratory PCA</h2>', pca_widget(qc, metadata), '</section>',
    '<section id="sample-qc"><h2>Library quality</h2><details><summary>Expand sample statistics and search</summary><label>Find a sample <input data-filter="sample-qc" placeholder="Sample ID or value" style="padding:8px;margin:8px"></label>', html_table(qc$samples), '</details></section>',
    '<section><h2>PCA variance</h2>', html_table(head(qc$variance, 6)),
    '<p>Coordinates: <a href="pca_coordinates.tsv">pca_coordinates.tsv</a>. PCA is descriptive and is not a group significance test.</p></section>',
    '<section><h2>Downloads</h2><ul>',
    '<li><a href="sample_qc.tsv">Library QC</a> / <a href="metadata_aligned.tsv">Aligned metadata</a></li>',
    '<li><a href="gene_filter.tsv">Gene filter decisions</a> / <a href="sample_correlations.tsv">Sample correlations</a></li>',
    if (has_sets) '<li><a href="program_scores.tsv">Within-cohort mean gene z-scores</a> / <a href="program_coverage.tsv">Gene-set coverage</a></li>',
    if (has_de) '<li><a href="differential_expression.tsv">DESeq2 results (case / reference)</a></li>',
    if (has_de) '<li><a href="model_matrix.tsv">Fitted design matrix</a> / <a href="model_samples.tsv">Contrast samples</a> / <a href="model_formula.txt">Model formula</a></li>',
    if (has_de) '<li><a href="de_summary.tsv">DE summary (FDR&lt;0.05, |LFC|≥1)</a> / <a href="de_gene_filter.tsv">Contrast-specific gene filter</a></li>',
    if (file.exists(file.path(out, "trait_associations.tsv"))) '<li><a href="trait_associations.tsv">Trait associations and BH-adjusted p-values</a> / <a href="trait_sample_inclusion.tsv">Trait missingness decisions</a></li>',
    '<li><a href="run_info.tsv">Parameters</a> / <a href="input_checksums.tsv">Input MD5</a> / <a href="code_checksums.tsv">Code MD5</a> / <a href="sessionInfo.txt">R session</a></li>',
    '</ul></section>',
    if (has_sets) c('<section><h2>Program coverage</h2>', html_table(read.delim(file.path(out, "program_coverage.tsv"))), '</section>'),
    if (has_de) c('<section><h2>Differential expression · first 30 rows</h2>', html_table(head(read.delim(file.path(out, "differential_expression.tsv")), 30)), '</section>'),
    if (file.exists(file.path(out, "trait_associations.tsv"))) c('<section><h2>Trait associations</h2>', html_table(read.delim(file.path(out, "trait_associations.tsv"))),
      '<p>Spearman estimates are rank correlations. Adjusted OLS estimates are score changes per trait unit; confidence intervals apply only to OLS. BH correction covers the programs tested for this trait in this run.</p></section>'),
    '<footer>MASLD–HCC Transcriptomics · Cohorts are analysed independently. Program scores are descriptive, not GSVA/ssGSEA or proof of YAP protein activity.</footer>',
    '<script>document.querySelectorAll("[data-filter]").forEach(input=>input.addEventListener("input",()=>{const q=input.value.toLowerCase();document.getElementById(input.dataset.filter).querySelectorAll("tbody tr").forEach(row=>{row.hidden=!row.textContent.toLowerCase().includes(q)});}));</script></html>')
  writeLines(body, file.path(out, "report.html"), useBytes = TRUE)
}

write_tsv <- function(x, path) {
  # Binary connection keeps LF bytes stable across Windows/Linux and Git checkouts.
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", eol = "\n")
}

demo_inputs <- function() {
  set.seed(20261117)
  n <- 12L
  metadata <- data.frame(sample_id = sprintf("DEMO_%02d", 1:n),
                         group = rep(c("adjacent", "tumor"), 6),
                         patient_id = rep(sprintf("P%02d", 1:6), each = 2))
  mu <- matrix(rep(seq(40, 400, length.out = 120), n), nrow = 120)
  mu[1:12, metadata$group == "tumor"] <- mu[1:12, metadata$group == "tumor"] * 3
  counts <- matrix(rnbinom(length(mu), mu = mu, size = 15), nrow = 120,
                   dimnames = list(sprintf("DEMO_GENE_%03d", 1:120), metadata$sample_id))
  list(counts = counts, metadata = metadata)
}
