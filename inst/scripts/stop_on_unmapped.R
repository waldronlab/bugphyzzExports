## Vocabulary check for discrete attribute values, sourced by
## export_bugphyzz.R and tested by test_stop_on_unmapped.R.
##
## stop_on_unmapped(phys) stops when a discrete value is not in bugphyzz's vocabulary
## (attributes.tsv). physiologies() has already mapped source terms with
## bugphyzz's attribute_synonyms.tsv; anything left would be dropped by
## taxPPro::filterData() with only a line in the log (#65).
stop_on_unmapped <- function(phys) {
    ## read exactly as taxPPro::filterDataDiscrete() reads it, so the two agree
    fname <- system.file("extdata", "attributes.tsv", package = "bugphyzz")
    if (!nzchar(fname)) stop("bugphyzz's attributes.tsv not found", call. = FALSE)
    vocab <- utils::read.table(fname, sep = "\t", header = TRUE)
    discrete <- c("binary", "multistate-intersection", "multistate-union")
    bad <- character(0)
    for (g in names(phys)) {
        df <- phys[[g]]
        if (!unique(df$Attribute_type)[1] %in% discrete) next
        rgx <- paste0("\\b", g, "\\b")
        valid <- vocab$attribute[grepl(rgx, vocab$attribute_group)]
        x <- setdiff(unique(df$Attribute), c(valid, NA))
        if (length(x))
            bad <- c(bad, paste0(g, ": ", paste0(x, collapse = " | ")))
    }
    if (length(bad))
        stop(
            "Attribute values not in bugphyzz's vocabulary; map them in ",
            "bugphyzz's inst/extdata/attribute_synonyms.tsv:\n",
            paste0(bad, collapse = "\n"), call. = FALSE
        )
    invisible(phys)
}
