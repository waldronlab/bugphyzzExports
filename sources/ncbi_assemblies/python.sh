#!/bin/sh
# Python with taxoniq and the pinned ncbi-taxon-db release, for build.py and test_build.py.
# Usage: sh python.sh <script.py> [args...]
# The taxonomy-only ncbi-taxon-db wheels need the lwaldron fork of taxoniq (accession indexes optional);
# both are pinned here and DB_VERSION in build.py must match the wheel.
exec uv run --no-project \
    --with "taxoniq @ git+https://github.com/lwaldron/taxoniq@cdaf0ebefc2ec76b0bbbf4e47d6bd60e6d94d4ca" \
    --with "https://github.com/lwaldron/taxoniq/releases/download/ncbi-taxon-db-2026.11.1/ncbi_taxon_db-2026.11.1-py3-none-any.whl" \
    python -I "$@"
