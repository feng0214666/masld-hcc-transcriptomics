#!/usr/bin/env bash
set -euo pipefail
raw=${1:-data/raw}
mkdir -p "$raw"
tail -n +2 config/study_sources.tsv | while IFS=$'\t' read -r file url expected || [[ -n "$file" ]]; do
  destination="$raw/$file"
  if [[ -s "$destination" ]] && printf '%s  %s\n' "$expected" "$destination" | sha256sum --check --status; then
    continue
  fi
  pending=$(mktemp "$raw/.download.XXXXXX")
  trap 'rm -f -- "$pending"' EXIT
  curl --fail --location --silent --show-error --retry 3 --retry-max-time 180 --connect-timeout 10 --max-time 120 --output "$pending" "$url"
  printf '%s  %s\n' "$expected" "$pending" | sha256sum --check --status || {
    echo "SHA-256 mismatch: $file" >&2
    exit 1
  }
  mv -- "$pending" "$destination"
  trap - EXIT
done
