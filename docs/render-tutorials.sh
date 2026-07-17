#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

for qmd in integration_4body_b2ddKpi.qmd lb2lc3pi-model.qmd cascade-vs-dpd.qmd decayamplitude-compatibility.qmd; do
    echo "Rendering ${qmd}..."
    quarto render "${qmd}" --to gfm
done

echo "Rendered tutorials:"
ls -1 *.qmd | while read -r qmd; do
    md="${qmd%.qmd}.md"
  if [[ -f "${md}" ]]; then
    echo "  ${md}"
  fi
done

echo "Check generated .md files for missing expected output or stray return values."
