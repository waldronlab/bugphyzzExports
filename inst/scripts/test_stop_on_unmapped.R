## Run from the repo root: Rscript inst/scripts/test_stop_on_unmapped.R
## Needs bugphyzz installed (reads its attributes.tsv).
source("inst/scripts/stop_on_unmapped.R")
ok <- function(values, group = "antimicrobial resistance", type = "binary") {
    list(x = data.frame(Attribute = values, Attribute_type = type,
                        stringsAsFactors = FALSE)) |> setNames(group)
}
## values in the vocabulary pass, and the list is returned invisibly.
## These two have been in attributes.tsv all along, so this test does not
## depend on the names added for #99; bugphyzz's own tests assert those.
res <- stop_on_unmapped(ok(c("resistance to beta-lactam", "sensitive to amikacin")))
stopifnot(is.list(res), nrow(res[[1]]) == 2)
## a value outside the vocabulary stops, and the message names it
e <- tryCatch(stop_on_unmapped(ok("sensitive to unobtainium")), error = identity)
stopifnot(inherits(e, "error"), grepl("unobtainium", conditionMessage(e)))
## NA is ignored by the check: the export guards against it separately
stop_on_unmapped(ok(c("resistance to beta-lactam", NA)))
## numeric attributes are skipped
stop_on_unmapped(ok("whatever", group = "gc content", type = "numeric"))
cat("ok\n")
