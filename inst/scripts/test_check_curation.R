#!/usr/bin/env Rscript
## Test for check_curation.R: run with `Rscript inst/scripts/test_check_curation.R`
## from the repository root. Needs bugphyzz installed (its extdata tables are
## the rules) and git; no network.
source("inst/scripts/check_curation.R")

dir <- file.path(tempfile("check_curation"), "sheets")
dir.create(dir, recursive = TRUE)
sheet <- function(name, lines, d = dir) {
    path <- file.path(d, name)
    writeLines(lines, path)
    path
}
head10 <- paste0(
    "NCBI_ID,Genome_ID,Accession_ID,Taxon_name,Attribute,Attribute_value,",
    "Attribute_source,Frequency,Note,Strain"
)
## problems as "line kind", to compare with what is expected
found <- function(p) paste(p$line, p$kind)
has <- function(p, line, pattern) any(p$line %in% line & grepl(pattern, p$message))

## a clean discrete file: a record over two lines, squished and upper-case
## values, a synonym on a TRUE row, a synonym mapped to nothing
clean <- sheet("aerophilicity.csv", c(
    head10,
    "562,,,Escherichia coli,facultatively anaerobic,TRUE,BacMap,usually,\"two",
    "lines\",K-12",
    "unknown,,,Some bacterium, aerobic ,FALSE,ProTraits,Always,,",
    "1280,,,Staphylococcus aureus,facultative aerobe,TRUE,BacMap,unknown,,"
))
p <- check_file(clean)
stopifnot(nrow(p) == 0)
p <- check_file(sheet("shape.csv", c(
    head10, "1,,,A b,rod-shaped,TRUE,BacMap,always,,",
    "1,,,A b,star-shaped,TRUE,BacMap,always,,"
)))
stopifnot(nrow(p) == 0)

## one problem of each kind, each on the line its record starts
dir2 <- tempfile("bad")
dir.create(dir2)
bad <- sheet("aerophilicity.csv", c(
    head10,                                                             # 1
    "562,,,Escherichia coli,aerobic,TRUE,BacMap,usually,\"two",         # 2
    "lines\",",                                                         # 3
    "abc,,,Escherichia coli,aerobic,TRUE,BacMap,usually,,",             # 4
    "562,,,Escherichia coli,flying,TRUE,BacMap,usually,,",              # 5
    "562,,,Escherichia coli,aerobic,yes,BacMap,usually,,",              # 6
    "562,,,Escherichia coli,aerobic,TRUE,BacMap,often,,",               # 7
    "562,,,Escherichia coli,aerobic,TRUE,Nobody_1900,usually,,",        # 8
    "562,,,Escherichia coli,aerobic,TRUE,BacMap",                       # 9
    "562,,,Escherichia coli,facultative aerobe,FALSE,BacMap,usually,,", # 10
    "562,,,,aerobic,TRUE,BacMap,usually,\"a",                           # 11
    "",                                                                 # 12
    "b\",",                                                             # 13
    "562,,,Escherichia coli,bacillus,TRUE,BacMap,usually,,",            # 14
    "562,,,Escherichia coli,,,BacMap,usually,,x,y"                      # 15
), d = dir2)
p <- check_file(bad)
stopifnot(identical(found(p), c(
    "4 value", "5 vocabulary", "6 value", "7 value", "8 value", "9 parse",
    "10 vocabulary", "11 value", "14 vocabulary", "15 parse"
)))
stopifnot(
    has(p, 4, "^NCBI_ID \"abc\" does not match"),
    has(p, 5, "^Attribute \"flying\" is not in bugphyzz's attributes.tsv$"),
    has(p, 6, "Attribute_value \"yes\" is not TRUE or FALSE"),
    has(p, 7, "^Frequency \"often\""),
    has(p, 8, "Attribute_source \"Nobody_1900\""),
    has(p, 9, "^7 fields, expected 10$"),
    ## synonyms are for TRUE rows only
    has(p, 10, "facultative aerobe"),
    has(p, 11, "^Taxon_name is empty$"), p$end[p$line == 11] == 13,
    ## a known attribute of another physiology
    has(p, 14, "\"bacillus\" is not in .* under \"aerophilicity\"$"),
    has(p, 15, "^11 fields, expected 10$")
)

## columns: missing, then out of order; both on line 1
p <- check_file(sheet("gram_stain.csv", c(
    "NCBI_ID,Genome_ID,Taxon_name,Attribute,Attribute_value,Attribute_source",
    "1,,A b,gram stain positive,TRUE,BacMap"
)))
stopifnot(
    identical(found(p), c("1 columns", "1 columns")),
    has(p, 1, "required column Accession_ID is missing"),
    has(p, 1, "required column Frequency is missing")
)
p <- check_file(sheet("gram_stain.csv", c(
    "Taxon_name,NCBI_ID,Genome_ID,Accession_ID,Attribute,Attribute_value,Attribute_source,Frequency",
    "A b,1,,,gram stain positive,TRUE,BacMap,always"
)))
stopifnot(identical(found(p), "1 columns"), has(p, 1, "out of order"))

## numeric and range values
head8 <- sub(",Note,Strain", "", head10)
p <- check_file(sheet("growth_temperature.csv", c(
    head8, paste0("1,,,A b,growth temperature,", c("37", "-2.5", "warm", "1-2", ""),
                  ",BacMap,unknown")
)))
stopifnot(identical(found(p), c("4 value", "5 value", "6 value")),
          has(p, 6, "Attribute_value is empty"))
p <- check_file(sheet("width.csv", c(
    head8, paste0("1,,,A b,width,", c("0.5-1", "<5", ">0.2", "0.5", "wide", "< 4"),
                  ",Bergey's Manual,always"),
    "1,,,A b,girth,1-2,Bergey's Manual,always"
)))
stopifnot(identical(found(p), c("6 value", "7 value", "8 vocabulary")))

## the wide fatty acid sheet: columns and template values only
p <- check_file(sheet("fatty_acid_composition.csv", c(
    "NCBI_ID,Genome_ID,Accession_ID,Taxon_name,Attribute,Attribute_source,Frequency,C10:0",
    "1,,,A b,fatty acid composition,Microbial Fatty Acid Compositions,unknown,3.2",
    "x,,,A b,fatty acid composition,Microbial Fatty Acid Compositions,unknown,3.2"
)))
stopifnot(identical(found(p), "3 value"))

## problems with no line
p <- check_file(sheet("no_such_trait.csv", c(head10, "1,,,A b,x,TRUE,BacMap,always,,")))
stopifnot(nrow(p) == 1, is.na(p$line), p$kind == "file")
p <- check_file(sheet("shape.csv", c(head10, "1,,,A b,bacillus,TRUE,BacMap,always,\"open,")))
stopifnot(nrow(p) == 1, is.na(p$line), p$kind == "parse")
p <- check_file(file.path(dir, "habitat.csv"))
stopifnot(nrow(p) == 1, grepl("not found", p$message))

## changed lines from a diff
stopifnot(identical(
    diff_added_lines(c(
        "diff --git a/x b/x", "@@ -3 +3 @@", "-a", "+b", "@@ -10,0 +11,2 @@ ctx",
        "+c", "+d", "@@ -20,2 +22,0 @@", "-e", "-f"
    )),
    c(3L, 11L, 12L)
))

## --base in a real repository: a file with an old bad row gets a new bad
## row, a new good one, and a new last line on the old two-line record
repo <- tempfile("repo")
dir.create(file.path(repo, "sources", "sheets"), recursive = TRUE)
g <- function(...)
    invisible(git(repo, c("-c", "user.name=t", "-c", "user.email=t@t", ...)))
f <- file.path(repo, "sources", "sheets", "aerophilicity.csv")
old <- c(
    head10,
    "562,,,Escherichia coli,flying,TRUE,BacMap,usually,,",              # 2
    ",,,Escherichia coli,aerobic,TRUE,BacMap,usually,\"two",            # 3
    "lines\","                                                          # 4
)
writeLines(old, f)
g("init", "-q", "-b", "main")
g("add", ".")
g("commit", "-q", "-m", "one")
g("checkout", "-q", "-b", "pr")
new <- c(
    old[1:3], "more lines\",",
    "562,,,Escherichia coli,aerobic,maybe,BacMap,usually,,",            # 5
    "562,,,Escherichia coli,anaerobic,FALSE,BacMap,usually,,"           # 6
)
writeLines(new, f)
writeLines(c(sub(",Strain", "", head10), "1,,,A b,coccus,TRUE,BacMap,always,"),
           file.path(repo, "sources", "sheets", "shape.csv"))
g("add", ".")
g("commit", "-q", "-m", "two")
stopifnot(identical(found(check_file(f)), c("2 vocabulary", "3 value", "5 value")))
stopifnot(identical(changed_lines(f, "main"), 4:6))
p <- check_curation(f, base = "main")
stopifnot(identical(found(p), c("3 value", "5 value")))
## uncommitted edits count too, and an unchanged header is not reported
writeLines(c(sub("Frequency", "Freq", new[1]), new[-1]), f)
stopifnot(identical(found(check_curation(f, base = "main")),
                    c("1 columns", "3 value", "5 value")))
writeLines(new, f)
owd <- setwd(repo)
stopifnot(identical(
    sort(changed_files("main")),
    file.path("sources/sheets", c("aerophilicity.csv", "shape.csv"))
))
setwd(owd)

## the command line: exit status, output, GitHub annotations
run <- function(args, env = character(0), wd = NULL) {
    script <- normalizePath("inst/scripts/check_curation.R")
    if (!is.null(wd)) { owd <- setwd(wd); on.exit(setwd(owd)) }
    out <- suppressWarnings(system2(
        file.path(R.home("bin"), "Rscript"), c(shQuote(script), shQuote(args)),
        stdout = TRUE, stderr = TRUE, env = env
    ))
    list(status = attr(out, "status") %||% 0L, out = out)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
r <- run(clean, "GITHUB_ACTIONS=false")
stopifnot(r$status == 0, grepl("^0 problem\\(s\\)", r$out))
r <- run(bad, "GITHUB_ACTIONS=false")
stopifnot(
    r$status == 1, !any(grepl("^::error", r$out)),
    paste0(bad, ":5: Attribute \"flying\" is not in bugphyzz's attributes.tsv") %in% r$out,
    grepl("^10 problem\\(s\\) in 1 of 1 file\\(s\\) checked \\(file 0, parse 2, columns 0, value 5, vocabulary 3\\)$",
          r$out[length(r$out)])
)
r <- run(bad, "GITHUB_ACTIONS=true")
stopifnot(r$status == 1, sum(grepl("^::error file=.*,line=[0-9]+::", r$out)) == 10)
## --base with no files finds the changed ones
r <- run(c("--base", "main"), "GITHUB_ACTIONS=false", wd = repo)
stopifnot(
    r$status == 1,
    "sources/sheets/aerophilicity.csv:5: Attribute_value \"maybe\" is not TRUE or FALSE" %in% r$out,
    !any(grepl("flying", r$out)),
    grepl("^2 problem\\(s\\) in 1 of 2 file\\(s\\)", r$out[length(r$out)])
)
r <- run(character(0))
stopifnot(r$status != 0, any(grepl("usage", r$out)))

cat("check_curation tests passed\n")
