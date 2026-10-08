# Fill before Zenodo deposit

The GitHub repository URL is already fixed:

    https://github.com/cagrisa-creator/pcr3bp-fiber-elimination

The remaining identifier is the Zenodo DOI. Do **not** publish the Zenodo record yet.

Recommended workflow:

1. Create a Zenodo draft manually.
2. Reserve a DOI.
3. Copy the reserved DOI into:
   - `paper/main.tex`
   - `CITATION.cff`
   - Zenodo metadata / related identifiers, as appropriate
4. Rebuild `paper/main.pdf` from the updated TeX source.
5. Commit the DOI-synchronized files.
6. Only then tag/release and publish the Zenodo record.

Useful replacement command after the DOI has been reserved:

    NEW_DOI="https://doi.org/10.5281/zenodo.NNNNNNN"
    grep -rl "Zenodo DOI" --include="*.tex" --include="*.md" .

At the current R04 review stage, the manuscript intentionally says that the verified Zenodo DOI and archival release are still pending. Replace that revision-stage note only after the DOI is actually reserved.
