> **v1.3.0:** on HAL, upload the corrected PDF as a new VERSION (V2) of the
> existing deposit (Modifier > Déposer une nouvelle version), keeping the
> same HAL id; put a one-line erratum summary in the comment field.

# Submitting to HAL

[HAL](https://hal.science) (Hyper Articles en Ligne, run by the CCSD/CNRS) is
France's national open archive. Foreign researchers and institutions are
welcome to deposit. This file walks through submitting this paper, and points
to the self-contained bundle in [`hal_submission/`](../hal_submission/).

## 1. Create or log in to your HAL account

Go to [hal.science](https://hal.science) and sign in (an ORCID-linked account
is convenient — HAL will offer to import your identity and past publications).

## 2. Start a new deposit

From your HAL dashboard: **Déposer** (Deposit) → **Déposer un document**. Since
this paper has not yet been peer-reviewed at the time of writing, choose a
document type of **Pré-print** or **Document de travail** (working paper); if
you later publish it in a journal, you can add a link between the HAL record
and the published version, or deposit a second, "article" record for the
accepted manuscript at that point (subject to the journal's self-archiving
policy — check via [Sherpa/Romeo](https://v2.sherpa.ac.uk/romeo/) before
depositing a version that has already been formally published elsewhere).

## 3. Upload the file

Upload `hal_submission/main.pdf` (or `paper/main.pdf` — same file) as the main
document. HAL recommends PDF for the primary deposit. You may optionally also
attach the LaTeX source as a supplementary file — `hal_submission/` includes a
ready-made `source_bundle.zip` (main.tex + figures/) for this purpose.

## 4. Choose a license — mandatory since February 2026

As of February 2026, HAL **requires** a license to be selected for every
deposited file (previously optional for text documents). The options are the
Creative Commons licenses, the Etalab open license, or Copyright/all-rights-
reserved. **CC-BY is the license recommended by the CCSD** for open-science
deposits — it permits reuse and adaptation, requiring only attribution. Unless
you have a specific reason to restrict reuse, select CC-BY when prompted
during the deposit form. This choice is essentially irrevocable once the
document is public, so decide deliberately; it applies to the version
deposited in HAL and does not constrain what license you later agree to with
a journal for the published version.

## 5. Fill in the metadata

The fields below are ready to copy-paste into HAL's deposit form; they are
also collected in [`hal_submission/METADATA.md`](../hal_submission/METADATA.md).

- **Title:** Exact fiber elimination for the Levi–Civita convexity gates of
  the planar restricted three-body problem, with a global convexity
  reconnaissance
- **Author:** César Augusto Grisa — affiliation: Instituto de Computação,
  Universidade Federal Fluminense (UFF), Niterói, RJ, Brazil
  (add your idHAL/ORCID if you have one, so the deposit attaches to your
  author record)
- **Abstract:** see `hal_submission/METADATA.md` (copied from the paper)
- **Keywords:** restricted three-body problem; Levi–Civita regularization;
  strict convexity; global surface of section; Birkhoff conjecture;
  computer-assisted proof
- **Scientific domain:** Mathematics / Dynamical Systems (`math.DS`); you may
  also tag Mathematical Physics if HAL's domain list offers it
- **Language:** English

## 6. Optional: transfer to arXiv

HAL offers a checkbox to automatically forward the deposit to arXiv. If you
want this paper on arXiv as well, tick it and provide an English abstract
(already the case here, since the paper is written in English) — arXiv
moderation still applies independently.

## 7. Submit for moderation

HAL deposits are moderated (by the CCSD team or your institution's HAL
correspondent) before going live — this is usually fast but is not
instantaneous. Once live, HAL assigns a permanent `hal-XXXXXXXX` identifier
and URL.

## 8. Cross-link with GitHub/Zenodo

Once you have the HAL URL, consider adding it to this repository's `README.md`
alongside the Zenodo DOI (see [`docs/ZENODO.md`](ZENODO.md)), so all three
records (GitHub, Zenodo, HAL) point to each other.

## What's in `hal_submission/`

```
hal_submission/
├── main.pdf           # the file to upload as the primary document
├── source_bundle.zip  # main.tex + figures/, for the optional source upload
└── METADATA.md         # every field above, ready to copy-paste
```
