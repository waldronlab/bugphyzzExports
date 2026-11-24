#!/usr/bin/env bash
# One-time migration tool: snapshot the bugphyzz curation Google Sheets into
# this directory as CSV (waldronlab/bugphyzz#305, #113).
#
# Do not re-run once bugphyzz reads these files instead of the sheets:
# a re-run would overwrite the frozen snapshot with whatever the sheets
# serve at that moment.
#
# Downloads the CSV each sheet serves through its publish-to-the-web link,
# as listed in bugphyzz's inst/extdata/spreadsheet_links.tsv and
# spreadsheet_customlinks.tsv, to <physiology>.csv, with spaces as
# underscores. The CSVs come with CRLF line endings and no final newline;
# this script adds the final newline and
# .gitattributes has git store the files with LF.
set -euo pipefail
cd "$(dirname "$0")"

# Pinned to the bugphyzz commit the README cites, so the documented snapshot
# can be reproduced; devel's link column will point back at this directory
# after waldronlab/bugphyzz#306.
BUGPHYZZ_RAW=${BUGPHYZZ_RAW:-https://raw.githubusercontent.com/waldronlab/bugphyzz/e33fb4a/inst/extdata}

fetch() { # fetch <physiology> <url>
    local file="${1// /_}.csv"
    echo "$file"
    curl -fsSL "$2" -o "$file"
    [ -n "$(tail -c 1 "$file")" ] && echo >> "$file"
    true
}

# spreadsheet_links.tsv: physiology in column 1, CSV link in column 5
curl -fsSL "$BUGPHYZZ_RAW/spreadsheet_links.tsv" |
    awk -F'\t' 'NR > 1 && $5 != "" { print $1 "\t" $5 }' |
    while IFS=$'\t' read -r physiology link; do
        fetch "$physiology" "$link"
    done

# spreadsheet_customlinks.tsv: physiology in column 1, link in column 3
curl -fsSL "$BUGPHYZZ_RAW/spreadsheet_customlinks.tsv" |
    awk -F'\t' 'NR > 1 && $3 != "" { print $1 "\t" $3 }' |
    while IFS=$'\t' read -r physiology link; do
        fetch "$physiology" "$link"
    done
