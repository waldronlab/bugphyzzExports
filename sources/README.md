# Programmatic sources

Curation data that comes from a script rather than a curator, plus `sheets/`, a snapshot of the hand-curation Google Sheets (ADR 0007). Each programmatic source has its own directory with the script that produced it, its output as one TSV per physiology in the curation-sheet column layout, and a README with provenance. `bugphyzz::physiologies()` reads `manifest.tsv` and the files it lists from the `devel` branch of this repository, and binds them to the sheet rows (read from `sheets/` once waldronlab/bugphyzz#306 is merged, until then live from Google).

Rules:
- Never edit a `*.tsv` output by hand. Put corrections in the source's `overrides.tsv` so a re-run keeps them.
- Rows in `sheets/` copied from a cited source are not edited either, and neither are the frozen Google Sheets they came from. A wrong row gets a line in `sheet_overrides.tsv` (drop only, with a reason naming the issue); the export drops it and logs it, and `physiologies()` still shows it (ADR 0005).
- Pin upstream versions (commit, record, model) in the script and the README.
- Look taxids up with taxoniq and a pinned `ncbi-taxon-db` release (https://github.com/lwaldron/taxoniq/releases), as `ncbi_assemblies/python.sh` does, rather than a raw taxdump. taxoniq has no synonym lookup; where old names must be matched, use `names.dmp` from the NCBI `taxdump_archive` file of the same date as the release, as `microbe_directory/` does.
- Add a row to `manifest.tsv` for each new file; `attribute_type` must match `bugphyzz/inst/extdata/spreadsheet_links.tsv` for that physiology. `sheets/` is exempt: `physiologies()` reads those files as the sheet data itself, so a manifest row would read every row twice.
- Attribute values must be in `bugphyzz/inst/extdata/attributes.tsv`; the export drops anything else.

| Source | What | Script |
|---|---|---|
| `bergeys_2026/` | Bergey's Manual genus-chapter abstracts, extracted with an LLM | `extract.py` then `build_traits.py` |
| `madin_2020/` | Madin et al. 2020 condensed species traits | `import.R` |
| `bacdive/` | Every BacDive strain record from the public API | `run.sh` (`fetch.py` then `build.py`) |
| `ncbi_assemblies/` | Genome size, G+C content and coding genes of one RefSeq assembly per species | `run.sh` (`build.py`) |
| `microbe_directory/` | The Microbe Directory 2018 table, eleven traits | `run.sh` (`build.py`) |
| `wikipedia_2026/` | Wikipedia introductory paragraphs carried by taxoniq, extracted with an LLM (preliminary, waldronlab/bugphyzz#296) | `collect.py`, `extract.py`, then `build_traits.py` |
| `sheets/` | Snapshot of the hand-curation Google Sheets, one CSV per physiology (not in the manifest) | `snapshot.sh`, one-time |
