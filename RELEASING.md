# Releasing bugphyzz

A release has two halves: a data release here (a tagged export, deposited on Zenodo) and a package release in waldronlab/bugphyzz (the new deposit as the `importBugphyzz()` default, pushed to Bioconductor). Do them in this order, before the Bioconductor deadline for packages passing build and check ([release schedule](https://bioconductor.org/developers/release-schedule/)). Track the work in a checklist issue in waldronlab/bugphyzz.

## 1. Validation, only if curated inputs changed

Re-run the cross-validation in waldronlab/taxPProValidation when sources, sheets or the vocabulary changed enough to move the scores; a re-run of the export alone does not need it. It runs on the supermicro for several hours: use a private R library with bugphyzz devel, run `bugphyzz::physiologies("all", forceDownload = TRUE)` first, and run detached. Copy `validation_attributes.tsv`, `validation_states.tsv` and `validation_summary.tsv` into `bugphyzz/inst/extdata/` by a bugphyzz PR, with the commit recorded in `inst/scripts/README.md`.

## 2. Export

Merge every data and rule change into `devel`, then run the export:

```sh
gh workflow run export-bugphyzz.yml -R waldronlab/bugphyzzExports --ref devel
```

It takes about two hours, installs bugphyzz from its `devel` branch, and commits the files to `devel` as "Update exports". Read `log_file.log`: the number of imputed annotations kept per attribute should match the validation tables.

## 3. Check the data

With bugphyzz devel installed, run its network tests against the new export:

```r
testthat::test_file("tests/testthat/test-importBugphyzz.R")  # from the bugphyzz checkout
```

`checkCuration()` must pass for every attribute. A failure means a value outside `attributes.tsv`; fix the vocabulary in bugphyzz, merge, and go back to step 2.

## 4. Release commit

Open a PR on `devel` that changes only:

- `NEWS.md`: date and name the section, `## YYYY-MM-DD, vX.Y.Z (Zenodo DOI added after the release)`;
- `DESCRIPTION`: `Version: X.Y.Z`;
- `.zenodo.json`: creators and affiliations.

Merge it after the export from step 2 has landed, so the release commit sits on top of the files to be deposited. The weekly export runs on Sundays at 22:00 UTC; tag straight after the merge, or disable the workflow until the tag is pushed (`gh workflow disable export-bugphyzz.yml -R waldronlab/bugphyzzExports`, then `enable`).

## 5. Tag and publish

```sh
git tag -a vX.Y.Z <merge-commit> -m "bugphyzz data vX.Y.Z"
git push origin vX.Y.Z
# release notes: the vX.Y.Z section of NEWS.md, without its heading
awk '/^## .*vX\.Y\.Z/{f=1; next} /^## /{f=0} f' NEWS.md > notes.md
gh release create vX.Y.Z -R waldronlab/bugphyzzExports --title "bugphyzz vX.Y.Z" --notes-file notes.md
```

Publishing the release makes the Zenodo deposit under the concept DOI 10.5281/zenodo.10980653 within a minute or two.

## 6. Check the deposit

Download the zip from the new record. The three `bugphyzz_*.csv` files, `bugphyzz_asr_dropped.csv.gz` and the GMT files must be the data, megabytes each, not Git LFS pointers of about 130 bytes. Then:

```r
bp <- bugphyzz::importBugphyzz(version = "10.5281/zenodo.<record>", forceDownload = TRUE)
```

must return the same number of rows as `version = "devel"`. Then, in a separate commit after the tag, replace "Zenodo DOI added after the release" in the NEWS heading with the record DOI.

## 7. Package

In waldronlab/bugphyzz, one PR:

- the record DOI as the default `version` of `importBugphyzz()`, in the roxygen text and `man/importBugphyzz.Rd`;
- the pinned-DOI test in `tests/testthat/test-importBugphyzz.R`;
- a version bump and a NEWS entry.

NEWS keeps one entry per Bioconductor release: before the release branch is created, collapse the development entries into one `# bugphyzz X.Y.0` section (`Y` even) with short entries under `SIGNIFICANT USER-VISIBLE CHANGES`, `NEW FEATURES` and `BUG FIXES AND MINOR IMPROVEMENTS`.

After the merge, push to Bioconductor. Bioconductor's `devel` must be an ancestor of GitHub's:

```sh
git fetch upstream devel
git merge-base --is-ancestor upstream/devel devel && git push upstream devel:devel
```

If it is not, merge `upstream/devel` into `devel` on GitHub first. Check the next day's build report.

## 8. Manuscript

Re-render the bugphyzzAnalyses articles with `bugphyzz_version` set to the record DOI, and the slides in waldronlab/bugphyzz-manuscript with `refresh: true`. Every number quoted in the manuscript comes from that deposit.
