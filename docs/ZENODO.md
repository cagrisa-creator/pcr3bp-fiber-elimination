> **v1.3.0:** deposit this corrected version as a NEW Zenodo version of the
> same record (Zenodo versioning), so the v1 DOI keeps resolving with a
> visible newer version; mention ERRATUM.md in the version notes.

# Archiving this repository on Zenodo (for the DOI in the paper)

The paper's Data and Code Availability statement has two placeholders:

```
...publicly available on GitHub at https://github.com/seu-usuario/seu-repositorio.
...archived and is available on Zenodo under DOI https://doi.org/10.5281/zenodo.XXXXXXX.
```

This file walks through filling in both. `CITATION.cff` in the repository root
is already prepared with the paper's metadata (title, author, ORCID field
ready to fill, license, keywords) — Zenodo reads it automatically, so you
should not need to retype anything on the Zenodo side.

## 1. Push this repository to GitHub

```bash
cd pcr3bp-fiber-elimination
git init
git add .
git commit -m "Initial release: exact fiber elimination for PCR3BP convexity gates"
git branch -M main
git remote add origin https://github.com/YOUR-USERNAME/YOUR-REPOSITORY.git
git push -u origin main
```

Replace `YOUR-USERNAME/YOUR-REPOSITORY` with your actual GitHub path.

## 2. Update the GitHub URL in the paper

Edit `paper/main.tex`, search for the "Data and Code Availability" section, and
replace `https://github.com/seu-usuario/seu-repositorio` with the real URL from
step 1. Do the same in `paper/main_docx_source.tex` and regenerate
`paper/main.docx` (see `paper/DOCX_NOTES.md`). Recompile the PDF.

Commit and push this change before continuing, so the archived Zenodo snapshot
(next step) already contains the correct GitHub URL.

## 3. Connect Zenodo to your GitHub account

1. Go to [zenodo.org](https://zenodo.org) and log in (you can sign in directly
   with your GitHub account, or link GitHub from your Zenodo profile page).
2. Go to your GitHub-linked repositories page in Zenodo
   (Profile → GitHub, or `https://zenodo.org/account/settings/github/`).
3. Find `YOUR-USERNAME/YOUR-REPOSITORY` in the list and toggle it **On**.
   This sets up a webhook on the GitHub repository; Zenodo will now archive
   every new GitHub release automatically.

## 4. Make a GitHub release

A Zenodo archive is created from a GitHub **release** (a tagged snapshot), not
from ordinary commits.

1. On GitHub, go to the repository → **Releases** → **Create a new release**.
2. Tag: `v1.0.0` (or your preferred version).
3. Title and description: e.g. "Initial release: exact fiber elimination for
   the PCR3BP convexity gates".
4. Click **Publish release**.

Within a minute or two, Zenodo downloads the repository as a zip, creates a
new record, and mints a DOI, using `CITATION.cff` (and `LICENSE`) for the
metadata.

## 5. Retrieve the DOI and finish the paper

1. On your Zenodo GitHub integration page, click through to the new record
   (or check `https://zenodo.org/me/uploads`).
2. Copy the DOI shown, e.g. `10.5281/zenodo.1234567`.
3. Edit `paper/main.tex` (and `main_docx_source.tex` → regenerate `main.docx`)
   and replace `10.5281/zenodo.XXXXXXX` with the real DOI. Recompile the PDF.
4. Optional but recommended: add the Zenodo DOI badge to `README.md`
   (Zenodo shows the exact Markdown snippet on the record page).

## 6. New versions later

If you revise the code after publication, tag a new GitHub release
(e.g. `v1.1.0`); Zenodo automatically creates a new version of the same
record, with its own DOI, while a "concept DOI" continues to resolve to the
latest version. You do not need to repeat the toggle step.

## Notes

- Pre-reserving a DOI *before* the first GitHub release is not supported via
  the GitHub integration; if you need the DOI in hand before making the
  release public, use Zenodo's manual upload flow instead (New Upload →
  Reserve DOI), then later switch to the GitHub-integration workflow for
  subsequent versions.
- Zenodo's default license metadata is overridden by the `LICENSE` file found
  in the repository root (already MIT for the code in this repository, as
  described in `LICENSE`).
