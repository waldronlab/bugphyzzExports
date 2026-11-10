## Validation filter for imputed (ASR) annotations, sourced by
## export_bugphyzz.R and tested by test_filter_asr.R. Implements ADR 0006.
##
## filterAsr(df, validation, kappa_threshold, r2_threshold, score_threshold)
## takes one attribute's rows and returns list(kept, dropped, n_before,
## n_after). An attribute's imputed states are kept or dropped together: an
## ASR row is kept when the attribute's cross-validation score (lower 95%
## bound of kappa, or of R2 for numeric attributes) is at least the
## threshold for its type and the row's Score is at least score_threshold.
## `dropped` holds the ASR rows that fail either rule, with the attribute's
## score in a `Validation` column, for bugphyzz_asr_dropped.csv.gz. Curated
## rows are never filtered. `validation` has one row per attribute:
## `.attr_key` (physiology, lower case) and `value` (the score bound).

filterAsr <- function(df, validation, kappa_threshold, r2_threshold,
                      score_threshold) {
    attr_type <- unique(df$Attribute_type)
    attr_type <- attr_type[!is.na(attr_type)]
    is_numeric <- any(attr_type %in% c("numeric", "range"))
    thr <- if (is_numeric) r2_threshold else kappa_threshold
    df <- df |>
        dplyr::mutate(.attr_key = tolower(Attribute)) |>
        dplyr::left_join(validation, by = ".attr_key") |>
        dplyr::select(-.attr_key)
    is_asr <- df$Evidence %in% "asr"
    passes <- !is.na(df$value) & df$value >= thr &
        !is.na(df$Score) & df$Score >= score_threshold
    keep <- !is_asr | (passes %in% TRUE)
    list(
        kept = dplyr::select(df[keep, , drop = FALSE], -value),
        dropped = dplyr::rename(df[!keep, , drop = FALSE], Validation = value),
        n_before = sum(is_asr),
        n_after = sum(is_asr & keep)
    )
}
