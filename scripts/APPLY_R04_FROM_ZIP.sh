#!/usr/bin/env bash
set -euo pipefail

# Apply the local CMDA R04 package to this repository.
# Usage from the repository root:
#   bash scripts/APPLY_R04_FROM_ZIP.sh /path/to/CMDA_FIBER_R04_AUTOCITACAO_20261008.zip
#
# This script does not publish a release and does not create a Zenodo DOI.

ZIP_PATH="${1:-}"
if [ -z "$ZIP_PATH" ]; then
  echo "Usage: bash scripts/APPLY_R04_FROM_ZIP.sh /path/to/CMDA_FIBER_R04_AUTOCITACAO_20261008.zip" >&2
  exit 2
fi
if [ ! -f "$ZIP_PATH" ]; then
  echo "ERROR: zip not found: $ZIP_PATH" >&2
  exit 2
fi

WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT
unzip -q "$ZIP_PATH" -d "$WORKDIR"
R04="$WORKDIR/CMDA_FIBER_R04_20261008"
if [ ! -d "$R04" ]; then
  echo "ERROR: expected folder not found inside zip: CMDA_FIBER_R04_20261008" >&2
  exit 2
fi

# Preserve a backup of historical DOCX/HAL artifacts if present.
mkdir -p _historical_v132_not_r04
for p in paper/main.docx paper/main_docx_source.tex paper/DOCX_NOTES.md hal_submission/METADATA.md hal_submission/main.pdf hal_submission/source_bundle.zip; do
  if [ -e "$p" ]; then
    mkdir -p "_historical_v132_not_r04/$(dirname "$p")"
    cp -a "$p" "_historical_v132_not_r04/$p"
  fi
done

# Copy R04 manuscript and revised figures.
cp "$R04/working/paper/main.tex" paper/main.tex
cp "$R04/working/paper/main.pdf" paper/main.pdf
cp "$R04/working/paper/figures/phase_map.png" paper/figures/phase_map.png
cp "$R04/working/paper/figures/dynamical_relevance.png" paper/figures/dynamical_relevance.png
cp "$R04/working/paper/figures/lobe_structure.png" paper/figures/lobe_structure.png

# Copy code snapshot shipped with R04. This preserves the numerical engines used in the audit.
cp "$R04/working/code/"*.py code/

# Do not keep stale DOCX/HAL deliverables as if they were R04 deliverables.
rm -f paper/main.docx paper/main_docx_source.tex paper/DOCX_NOTES.md
rm -f hal_submission/METADATA.md hal_submission/main.pdf hal_submission/source_bundle.zip

# Copy audit logs into a private/revision folder for traceability.
mkdir -p docs/revision_r04_audit
cp "$R04/00_START_HERE_R04.md" docs/revision_r04_audit/
cp "$R04/audit/"*.log docs/revision_r04_audit/ 2>/dev/null || true
cp "$R04/deliverables/AUDITORIA_R04_AUTOCITACAO_PT.pdf" docs/revision_r04_audit/ 2>/dev/null || true

# Rebuild once locally if pdflatex is available.
if command -v pdflatex >/dev/null 2>&1; then
  (cd paper && pdflatex -interaction=nonstopmode main.tex >/tmp/cmda_r04_pdflatex1.log && pdflatex -interaction=nonstopmode main.tex >/tmp/cmda_r04_pdflatex2.log)
fi

# Basic checks.
grep -R "companion manuscript\|GrisaCertificate\|129,238 subcells\|Zenodo.XXXXXXX\|seu-usuario" paper README.md CITATION.cff 2>/dev/null && {
  echo "ERROR: stale placeholder/dependency text found above." >&2
  exit 1
} || true

sha256sum paper/main.tex paper/main.pdf paper/figures/*.png > R04_FINAL_LOCAL_SHA256SUMS.txt

echo "R04 files applied locally. Review 'git diff', run tests, then commit."
echo "Suggested next commands:"
echo "  git status"
echo "  git diff --stat"
echo "  cd code && python fiber_elimination_audit.py && python rotation_term_audit.py && cd .."
echo "  git add paper code README.md CITATION.cff FILL_BEFORE_DEPOSIT.md docs scripts R04_FINAL_LOCAL_SHA256SUMS.txt _historical_v132_not_r04"
echo "  git commit -m 'Apply CMDA R04 author-review manuscript package'"
echo "  git push origin cmda-r04-author-review-20261008"
