#!/usr/bin/env Rscript
source("R/workflow.R"); source("R/geo.R"); source("R/import_geo.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: analysis/import_gse126848.R counts.txt.gz family.soft.gz NEW_OUTPUT_DIR")
import_gse126848(args[1], args[2], args[3])
