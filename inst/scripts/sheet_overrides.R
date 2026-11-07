## Record-level corrections to rows the Google Sheets copy from other sources
## (ADR 0005). The sheets stay as they are and physiologies() still shows the
## rows; the export drops them. Each row of sources/sheet_overrides.tsv drops
## the rows of one physiology that match NCBI_ID, Attribute and
## Attribute_source, and Attribute_value too when it is filled in. A row that
## matches nothing stops the export, so stale corrections show up.

read_sheet_overrides <- function(path = "sources/sheet_overrides.tsv") {
    ov <- utils::read.delim(
        path, colClasses = "character", na.strings = character(0)
    )
    need <- c(
        "physiology", "NCBI_ID", "Taxon_name", "Attribute", "Attribute_value",
        "Attribute_source", "reason"
    )
    miss <- setdiff(need, names(ov))
    if (length(miss))
        stop("sheet_overrides.tsv lacks columns: ", toString(miss), call. = FALSE)
    if (any(!nzchar(ov$reason)))
        stop("every sheet override needs a reason", call. = FALSE)
    ov
}

## `l`: named list of data frames from physiologies(). Returns `l` without
## the matching rows; `log` is called with one line per override.
apply_sheet_overrides <- function(l, ov, log = message) {
    for (i in seq_len(nrow(ov))) {
        o <- ov[i, ]
        df <- l[[o$physiology]]
        if (is.null(df))
            stop("sheet override for a physiology not imported: ",
                 o$physiology, call. = FALSE)
        hit <- as.character(df$NCBI_ID) == o$NCBI_ID &
            df$Attribute == o$Attribute &
            df$Attribute_source == o$Attribute_source
        if (nzchar(o$Attribute_value)) {
            if (!"Attribute_value" %in% names(df))
                stop("sheet override gives Attribute_value but ",
                     o$physiology, " has none", call. = FALSE)
            hit <- hit & as.character(df$Attribute_value) == o$Attribute_value
        }
        hit <- !is.na(hit) & hit
        if (!any(hit))
            stop("sheet override matched nothing: ",
                 paste(o[1, 1:6], collapse = " | "), call. = FALSE)
        l[[o$physiology]] <- df[!hit, , drop = FALSE]
        log(paste0(
            o$physiology, ": dropped ", sum(hit), " ", o$Attribute_source,
            " row(s) for ", o$Taxon_name, " (", o$NCBI_ID, "): ", o$reason
        ))
    }
    l
}
