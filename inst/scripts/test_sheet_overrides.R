#!/usr/bin/env Rscript
## Test for sheet_overrides.R: run with `Rscript inst/scripts/test_sheet_overrides.R`
## from the repository root. Also checks that sources/sheet_overrides.tsv parses.
source("inst/scripts/sheet_overrides.R")
quiet <- function(x) invisible(NULL)
l <- list(
    `gram stain` = data.frame(
        NCBI_ID = c("1", "1", "1", "2"), Taxon_name = "x",
        Attribute = c("gram stain positive", "gram stain positive",
                      "gram stain positive", "gram stain negative"),
        Attribute_value = c(TRUE, FALSE, TRUE, TRUE),
        Attribute_source = c("ProTraits", "ProTraits", "BacMap", "ProTraits")
    ),
    `genome size` = data.frame(
        NCBI_ID = c("3", "4"), Taxon_name = "y", Attribute = "Genome size",
        Attribute_value_min = c(53560, 3e6), Attribute_value_max = c(53560, 3e6),
        Attribute_source = "Kegg"
    )
)
ov <- function(...) {
    d <- data.frame(
        physiology = "gram stain", NCBI_ID = "1", Taxon_name = "x",
        Attribute = "gram stain positive", Attribute_value = "",
        Attribute_source = "ProTraits", reason = "test"
    )
    a <- list(...)
    d[names(a)] <- a
    d
}
fails <- function(expr, pattern) {
    msg <- tryCatch({ expr; "" }, error = conditionMessage)
    stopifnot(grepl(pattern, msg))
}

## value given: only the TRUE ProTraits row goes; FALSE and BacMap stay
out <- apply_sheet_overrides(l, ov(Attribute_value = "TRUE"), log = quiet)
g <- out[["gram stain"]]
stopifnot(nrow(g) == 3, !any(g$Attribute_source == "ProTraits" &
    g$NCBI_ID == "1" & g$Attribute_value))
## value blank: every ProTraits row for that taxon and attribute goes
out <- apply_sheet_overrides(l, ov(), log = quiet)
stopifnot(nrow(out[["gram stain"]]) == 2)
## numeric physiology, no Attribute_value column
out <- apply_sheet_overrides(l, ov(physiology = "genome size", NCBI_ID = "3",
    Taxon_name = "y", Attribute = "Genome size", Attribute_source = "Kegg"),
    log = quiet)
stopifnot(identical(out[["genome size"]]$NCBI_ID, "4"))
## error paths
fails(apply_sheet_overrides(l, ov(NCBI_ID = "9"), log = quiet), "matched nothing")
fails(apply_sheet_overrides(l, ov(Attribute_source = "Kegg"), log = quiet), "matched nothing")
fails(apply_sheet_overrides(l, ov(physiology = "shape"), log = quiet), "not imported")
fails(apply_sheet_overrides(l, ov(physiology = "genome size", NCBI_ID = "3",
    Attribute = "Genome size", Attribute_source = "Kegg",
    Attribute_value = "53560"), log = quiet), "has none")
## the committed file parses
invisible(read_sheet_overrides("sources/sheet_overrides.tsv"))
cat("sheet_overrides tests passed\n")
