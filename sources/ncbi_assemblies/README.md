# NCBI RefSeq assemblies: genome size, G+C content and coding genes

One RefSeq assembly per bacterial or archaeal species, with the genome size, G+C percent and protein-coding gene count NCBI reports for it. Supersedes the Kegg genome size and coding gene rows (4,303 and 4,298 species, with plasmid-sized outliers, waldronlab/bugphyzz#262) and extends Madin (14,149 species) to every species with a RefSeq genome.

Citation: Goldfarb T, Kodali VK, Pujar S, et al. NCBI RefSeq: reference sequence standards through 25 years of curation and annotation. Nucleic Acids Res. 2025;53(D1):D243-D257. doi:10.1093/nar/gkae1038. NCBI data carry no copyright restrictions (https://www.ncbi.nlm.nih.gov/home/about/policies/).

## Provenance

- Upstream: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt (columns described in https://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt). Downloaded 5 November 2026; 569,714 assemblies, 551,414 of them bacterial or archaeal.
- Selection (`build.py`): bacterial and archaeal assemblies with `excluded_from_refseq` empty (drops 50,007 flagged assemblies, mostly derived from metagenomes), full genome representation and the latest version; 501,407 candidates. One per `species_taxid`: the RefSeq reference genome first, then an assembly from type material (or a designated neotype or reftype), then the most complete assembly level (Complete Genome, Chromosome, Scaffold, Contig), then the newest release date. No filter on genome size.
- Taxids: each `species_taxid` is checked with taxoniq against ncbi-taxon-db 2026.11.1 (https://github.com/lwaldron/taxoniq/releases/tag/ncbi-taxon-db-2026.11.1, built from the NCBI new_taxdump of 1 November 2026); its current scientific name is the `Taxon_name`. 319 species taxids newer than that release were dropped, and one no longer at species rank.
- Values are copied from the summary: `genome_size` (total length, base pairs), `gc_percent` (mol%, one decimal), `protein_coding_gene_count` (genes). `Note` holds the assembly accession, its level and its reference or type-material flags.
- Nothing is edited by hand.

## Output

| File | Species | Notes |
|---|---|---|
| `genome_size.tsv` | 77,673 | `Genome size`, base pair |
| `gc_content.tsv` | 77,673 | `gc content`, mol%; one value, not a range |
| `coding_genes.tsv` | 77,673 | `genes`, protein-coding genes |

Of the 77,673 assemblies chosen, 23,617 are reference genomes and 19,851 are from type material; 15,601 are complete genomes, 1,110 chromosome level, 26,392 scaffolds and 34,570 contigs. 53,551 species are NCBI placeholder species (`Pseudomonas sp. XYZ`); they are kept because bugphyzz keys on the taxid.

## Refreshing

```sh
sh run.sh <scratch_dir>      # downloads the summary (about 250 MB, not committed) and builds the tables
sh python.sh test_build.py   # fixture test, no network
```

`python.sh` runs Python through `uv run` with the lwaldron fork of taxoniq (the taxonomy-only `ncbi-taxon-db` wheels need it) and the pinned wheel. To move to a newer taxonomy, change the wheel in `python.sh` and `DB_VERSION` in `build.py` together, then update the dates and counts above from the build output.
