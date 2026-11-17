# Dependency-free report widgets; these are exploratory dashboard views.
pca_widget <- function(qc, metadata) {
  d <- qc$pca
  group <- metadata$group[match(d$sample_id, metadata$sample_id)]
  groups <- sort(unique(group))
  colors <- rep(c("#147d92", "#db7446", "#6373ba", "#248866", "#bb5284", "#8c751b"), length.out = length(groups))
  names(colors) <- groups
  limits <- function(v) { r <- range(v); if (diff(r) == 0) r <- r + c(-1, 1); r + c(-1, 1) * diff(r) * .08 }
  xr <- limits(d$PC1); yr <- limits(d$PC2)
  px <- function(v) 70 + (v - xr[1]) / diff(xr) * 490
  py <- function(v) 320 - (v - yr[1]) / diff(yr) * 270
  ticks <- paste(vapply(seq(0, 1, length.out = 5), function(t) {
    x <- 70 + t * 490; y <- 320 - t * 270
    sprintf('<path d="M %.2f 50 V 320 M 70 %.2f H 560" stroke="#dce7e9"/><text x="%.2f" y="343" text-anchor="middle">%.1f</text><text x="58" y="%.2f" text-anchor="end">%.1f</text>',
            x, y, x, xr[1] + t * diff(xr), y + 4, yr[1] + t * diff(yr))
  }, character(1)), collapse = "")
  points <- paste(vapply(seq_len(nrow(d)), function(i) {
    sprintf('<circle tabindex="0" cx="%.2f" cy="%.2f" r="5.5" fill="%s" stroke="white" stroke-width="1.4"><title>%s · %s · PC1 %.2f · PC2 %.2f</title></circle>',
            px(d$PC1[i]), py(d$PC2[i]), colors[[group[i]]], html_escape(d$sample_id[i]), html_escape(group[i]), d$PC1[i], d$PC2[i])
  }, character(1)), collapse = "")
  legend <- paste(vapply(groups, function(g) paste0('<span style="margin-right:20px"><b style="color:', colors[[g]], '">●</b> ', html_escape(g), '</span>'), character(1)), collapse = "")
  paste0('<div class="pca-widget"><svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 640 390" role="img" aria-label="Exploratory PCA by source group" style="width:100%;max-width:850px;font:12px system-ui;fill:#48616b">',
          ticks, '<path d="M 70 50 V 320 H 560" stroke="#607d86" fill="none"/>', points,
          sprintf('<text x="315" y="376" text-anchor="middle">PC1 (%.1f%%)</text><text transform="translate(18 185) rotate(-90)" text-anchor="middle">PC2 (%.1f%%)</text>',
                  qc$variance$variance_fraction[1] * 100, qc$variance$variance_fraction[2] * 100),
          '</svg><p>', legend, '</p><p>Hover or focus a sample for its ID. Colors follow supplied group labels; apparent separation alone is not a statistical test.</p></div>')
}
