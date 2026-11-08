## Run from the repo root: Rscript inst/scripts/test_filter_asr.R
suppressPackageStartupMessages(library(dplyr))
source("inst/scripts/filter_asr.R")
validation <- data.frame(
    physiology = c("aerophilicity", "aerophilicity", "motility", "length"),
    attribute = c("aerobic", "anaerobic", "motility", "length"),
    value = c(0.66, 0.85, 0.52, 0.29)
)
validation_whole <- dplyr::filter(validation, physiology == attribute)
aero <- data.frame(
    Attribute = "aerophilicity", Attribute_type = "multistate-intersection",
    Attribute_value = c("aerobic", "aerobic", "anaerobic", "anaerobic", "facultatively anaerobic", "aerobic"),
    Evidence = c("exp", "asr", "asr", "asr", "asr", NA),
    Frequency = c("always", "usually", "usually", "never", "usually", "always"),
    stringsAsFactors = FALSE
)
res <- filterAsr(aero, validation, validation_whole, 0.7, 0.5)
## curated rows and the NA-evidence row stay; anaerobic/usually passes;
## aerobic (0.66), anaerobic/never and the unvalidated facultative row drop
stopifnot(
    nrow(res$kept) == 3, nrow(res$dropped) == 3,
    nrow(res$kept) + nrow(res$dropped) == nrow(aero),
    res$n_before == 4, res$n_after == 1,
    !"value" %in% names(res$kept), "Validation" %in% names(res$dropped),
    identical(res$dropped$Attribute_value, c("aerobic", "anaerobic", "facultatively anaerobic")),
    identical(res$dropped$Validation, c(0.66, 0.85, NA)),   # never fails on frequency, not score
    !any(is.na(res$kept$Attribute_value))
)
## binary: validated as a whole, through validation_whole
mot <- data.frame(Attribute = "motility", Attribute_type = "binary",
    Attribute_value = c("TRUE", "FALSE"), Evidence = c("asr", "exp"),
    Frequency = "always", stringsAsFactors = FALSE)
res <- filterAsr(mot, validation, validation_whole, 0.7, 0.5)
stopifnot(nrow(res$kept) == 1, res$dropped$Validation == 0.52)
res <- filterAsr(mot, validation, validation_whole, 0.5, 0.5)
stopifnot(nrow(res$kept) == 2, nrow(res$dropped) == 0)
## numeric: R2 threshold
len <- data.frame(Attribute = "length", Attribute_type = "range",
    Attribute_value = "1-2", Evidence = "asr", Frequency = "always", stringsAsFactors = FALSE)
stopifnot(nrow(filterAsr(len, validation, validation_whole, 0.7, 0.5)$kept) == 0,
          nrow(filterAsr(len, validation, validation_whole, 0.7, 0.2)$kept) == 1)
## empty input
e <- filterAsr(aero[0, ], validation, validation_whole, 0.7, 0.5)
stopifnot(nrow(e$kept) == 0, nrow(e$dropped) == 0, e$n_before == 0)
## round trip through the gz side file as export_bugphyzz.R writes it
tmp <- tempfile(fileext = ".csv.gz")
d <- filterAsr(aero, validation, validation_whole, 0.7, 0.5)$dropped
con <- gzfile(tmp, "w"); cat("# header\n", file = con)
suppressWarnings(write.table(d, con, quote = TRUE, sep = ",", row.names = FALSE, append = TRUE, col.names = TRUE)); close(con)
back <- read.csv(gzfile(tmp), skip = 1, colClasses = "character")
stopifnot(nrow(back) == nrow(d), identical(names(back), names(d)))
cat("ok\n")
