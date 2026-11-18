#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 1 ]]; then
  echo "Usage: bash analysis/reproduce_study.sh NEW_OUTPUT_DIR" >&2
  exit 2
fi
if [[ -e "$1" || -L "$1" ]]; then
  echo "Output path already exists: $1" >&2
  exit 2
fi
raw=${MASLD_RAW_DIR:-data/raw}
bash analysis/fetch_study.sh "$raw"
Rscript --vanilla analysis/run_study.R "$raw" "$1"
