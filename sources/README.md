# Programmatic sources

Curation data that comes from a script rather than a curator. Each source has its own directory with the script that produced it, its output as one TSV per physiology in the curation-sheet column layout, and a README with provenance. `bugphyzz::physiologies()` reads `manifest.tsv` and the files it lists from the `devel` branch of this repository, and binds them to the Google Sheet rows.

Rules:
- Never edit a `*.tsv` output by hand. Put corrections in the source's `overrides.tsv` so a re-run keeps them.
- Pin upstream versions (commit, record, model) in the script and the README.
- Add a row to `manifest.tsv` for each new file; `attribute_type` must match `bugphyzz/inst/extdata/spreadsheet_links.tsv` for that physiology.
- Attribute values must be in `bugphyzz/inst/extdata/attributes.tsv`; the export drops anything else.
- Taxonomy lookups (rank, name, merged ids) use taxoniq with a pinned `ncbi-taxon-db` release from https://github.com/lwaldron/taxoniq/releases; record the release tag in the source's README. (The Bergey's and BacDive builders predate this rule and parse a raw taxdump; switch them when next touched.)

| Source | What | Script |
|---|---|---|
| `bergeys_2026/` | Bergey's Manual genus-chapter abstracts, extracted with an LLM | `extract.py` then `build_traits.py` |
| `madin_2020/` | Madin et al. 2020 condensed species traits | `import.R` |
