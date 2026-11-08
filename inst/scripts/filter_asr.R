## Validation filter for imputed (ASR) annotations, sourced by
## export_bugphyzz.R and tested by test_filter_asr.R.
##
## filterAsr(df, validation, validation_whole, mcc_threshold, r2_threshold)
## takes one attribute's rows and returns list(kept, dropped, n_before,
## n_after). `kept` is df without the ASR rows that fail: no validation
## entry, a score below the threshold for the attribute type, or Frequency
## "rarely"/"never". `dropped` is those rows with the score in a
## `Validation` column, for bugphyzz_asr_dropped.csv.gz. Curated rows
## are never filtered. `validation` has columns physiology, attribute,
## value (lower case); `validation_whole` is its physiology == attribute
## subset (binary and numeric attributes are validated as a whole).

filterAsr <- function(df, validation, validation_whole, mcc_threshold, r2_threshold) {
    attr_type <- unique(df$Attribute_type)
    attr_type <- attr_type[!is.na(attr_type)]
    df <- dplyr::mutate(df, .attr_key = tolower(Attribute))
    if (any(attr_type %in% c("numeric", "range", "binary"))) {
        thr <- if (any(attr_type == "binary")) mcc_threshold else r2_threshold
        df <- dplyr::left_join(
            df, dplyr::select(validation_whole, .attr_key = attribute, value),
            by = ".attr_key"
        )
    } else {
        thr <- mcc_threshold
        df <- df |>
            dplyr::mutate(.value_key = tolower(Attribute_value)) |>
            dplyr::left_join(
                dplyr::select(
                    validation, .attr_key = physiology,
                    .value_key = attribute, value
                ),
                by = c(".attr_key", ".value_key")
            ) |>
            dplyr::select(-.value_key)
    }
    df <- dplyr::select(df, -.attr_key)
    is_asr <- df$Evidence %in% "asr"
    passes <- !is.na(df$value) & df$value >= thr &
        !df$Frequency %in% c("rarely", "never")
    keep <- !is_asr | (passes %in% TRUE)
    list(
        kept = dplyr::select(df[keep, , drop = FALSE], -value),
        dropped = dplyr::rename(df[!keep, , drop = FALSE], Validation = value),
        n_before = sum(is_asr),
        n_after = sum(is_asr & keep)
    )
}
