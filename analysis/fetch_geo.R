#!/usr/bin/env Rscript
source("R/workflow.R")
source("R/geo.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript --vanilla analysis/fetch_geo.R GSE126848 data/metadata/GSE126848_run1")
fetch_geo_metadata(args[1], args[2])
