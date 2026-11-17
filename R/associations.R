# One independent observation per person; no implicit missing-value exclusion.
associate_programs <- function(scores, metadata, trait, covariates = character()) {
  if (!trait %in% names(metadata) || !is.numeric(metadata[[trait]])) stop("Trait must be a numeric metadata column.")
  if (trait %in% covariates) stop("Trait cannot also be an adjustment covariate.")
  if (anyDuplicated(metadata$sample_id)) stop("Duplicate metadata sample IDs.")
  if ("patient_id" %in% names(metadata) && anyDuplicated(metadata$patient_id[!is.na(metadata$patient_id)]))
    stop("Trait association requires independent subjects; repeated patients need a separate model.")
  validate_covariates(metadata, covariates, reserved = c("sample_id", "group", "patient_id", "score", "trait_value"))
  if (!all(scores$sample_id %in% metadata$sample_id)) stop("Score samples are absent from metadata.")
  if (anyDuplicated(scores[c("program", "sample_id")])) stop("Duplicate program/sample scores.")
  if (!is.numeric(scores$score) || any(!is.finite(scores$score))) stop("Program scores must be finite numeric values.")
  exclusion <- data.frame(sample_id = metadata$sample_id,
                          included = is.finite(metadata[[trait]]),
                          reason = ifelse(is.finite(metadata[[trait]]), "included", "missing_or_nonfinite_trait"))
  results <- lapply(unique(scores$program), function(program) {
    s <- scores[scores$program == program, ]
    m <- metadata[match(s$sample_id, metadata$sample_id), , drop = FALSE]
    d <- data.frame(score = s$score, trait_value = m[[trait]], m[, covariates, drop = FALSE], check.names = FALSE)
    d <- d[is.finite(d$trait_value), , drop = FALSE]
    row <- data.frame(program = program, trait = trait, n = nrow(d),
                       method = if (length(covariates)) "adjusted_OLS" else "Spearman",
                       estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                       pvalue = NA_real_, status = "insufficient_samples")
    if (nrow(d) < 5) return(row)
    if (length(unique(d$trait_value)) < 2 || sd(d$score) == 0) { row$status <- "constant_variable"; return(row) }
    if (!length(covariates)) {
      test <- cor.test(d$score, d$trait_value, method = "spearman", exact = FALSE)
      row$estimate <- unname(test$estimate); row$pvalue <- test$p.value
    } else {
      formula <- reformulate(c("trait_value", covariates), response = "score")
      mm <- tryCatch(model.matrix(formula, d), error = function(e) NULL)
      if (is.null(mm) || qr(mm)$rank < ncol(mm) || nrow(mm) <= ncol(mm) + 1) { row$status <- "nonestimable_design"; return(row) }
      fit <- lm(formula, data = d)
      coefficient <- summary(fit)$coefficients["trait_value", ]
      interval <- confint(fit, "trait_value", level = 0.95)
      row$estimate <- coefficient[1]; row$ci_low <- interval[1]; row$ci_high <- interval[2]; row$pvalue <- coefficient[4]
    }
    row$status <- "ok"; row
  })
  if (!length(results)) return(list(results = data.frame(), exclusions = exclusion))
  result <- do.call(rbind, results)
  result$padj <- p.adjust(result$pvalue, method = "BH")
  list(results = result, exclusions = exclusion)
}
