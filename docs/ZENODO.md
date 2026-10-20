# Same-record version update

Current published version: 1.4.0, DOI 10.5281/zenodo.23244862.
Concept DOI: 10.5281/zenodo.23244861. Candidate: 1.4.1.

Use the existing record's **New version** action after signing in to its
owner account. Do not create an unrelated New upload. Reserve/record the
new version-specific DOI shown by Zenodo; never infer it from record numbers.
Confirm the existing family and version lineage in the draft before upload.

Before publishing, replace candidate wording and update CITATION.cff,
paper/main.tex Data and Code Availability, the compiled paper/main.pdf and
the manifests/archives with the actual reserved DOI. Keep previous DOIs as
historical references only. Recompile and verify the identifier-only diff.
Confirm the actual GitHub commit/tag after updating the existing repository.

The repository has an older unmerged draft PR. Review/update that branch and
merge the revision into the existing default branch only with author approval.
Check whether GitHub-to-Zenodo archiving is enabled before releasing, to avoid
a second record family. Public metadata does not establish that setting.
Do not publish a GitHub release or Zenodo record without final authorization.
