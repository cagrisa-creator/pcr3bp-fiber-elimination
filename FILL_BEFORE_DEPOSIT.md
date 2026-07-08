# Fill before deposit (1 minute)

The three identifiers below are the depositor's own and are intentionally
placeholders in the package. From the repository root, run (with your values):

    NEW_GH="https://github.com/USER/REPO"
    NEW_DOI="https://doi.org/10.5281/zenodo.NNNNNNN"
    grep -rl "seu-usuario" --include="*.tex" --include="*.md" . | xargs sed -i "s|https://github.com/seu-usuario/seu-repositorio|$NEW_GH|g"
    grep -rl "zenodo.XXXXXXX" --include="*.tex" --include="*.md" . | xargs sed -i "s|https://doi.org/10.5281/zenodo.XXXXXXX|$NEW_DOI|g"
    # HAL id (hal-XXXXXXXX) appears only in docs/HAL.md instructions.

Then recompile the paper (3x pdflatex), regenerate main.docx (pandoc),
refresh hal_submission/ (main.pdf + source_bundle.zip) and results/SHA256SUMS.
