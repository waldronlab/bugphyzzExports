filterDirectParents <- function(df) {
    vct <- dplyr::case_when(
        df$Rank == 'strain' & df$Parent_rank == 'species' ~ TRUE,
        df$Rank == 'species' & df$Parent_rank == 'genus' ~ TRUE,
        df$Rank == 'genus' & df$Parent_rank == 'family' ~ TRUE,
        df$Rank == 'family' & df$Parent_rank == 'order' ~ TRUE,
        df$Rank == 'order' & df$Parent_rank == 'class' ~ TRUE,
        df$Rank == 'class' & df$Parent_rank == 'phylum' ~ TRUE,
        TRUE ~ FALSE
    )
    return(df[vct,])
}

calcParentScores <- function(df) {
    attr_type <- unique(df$Attribute_type)
    attr_group <- unique(df$Attribute_group)
    if (attr_type == 'logical') {
        output <- df |>
            dplyr::group_by(
                .data$Parent_name, .data$Parent_NCBI_ID, .data$Parent_rank,
                .data$Attribute_type, .data$Attribute_group
            ) |>
            dplyr::mutate(
                n = dplyr::n(),
                total = sum(.data$Score),
                prop = .data$Score / .data$total
            ) |>
            dplyr::count(.data$Attribute, wt = .data$prop, name = 'Score') |>
            dplyr::mutate(Score = round(.data$Score, 1)) |>
            dplyr::ungroup() |>
            dplyr::mutate(
                Evidence = 'asr',
                Attribute_value = TRUE,
                Attribute_source = NA,
                Note = NA,
                Frequency = taxPPro::scores2Freq(Score),
                Attribute_type = attr_type,
                Attribute_group = attr_group
            )
        colnames(output)[which(colnames(output) == 'Parent_name')] <- 'Taxon_name'
        colnames(output)[which(colnames(output) == 'Parent_NCBI_ID')] <- 'NCBI_ID'
        colnames(output)[which(colnames(output) == 'Parent_rank')] <- 'Rank'
        # output$Score <- NULL
        output <- unique(output)
        return(output)
    }

    if (attr_type == 'range') {
        output <- df |>
            dplyr::group_by(
                .data$Parent_name, .data$Parent_NCBI_ID, .data$Parent_rank,
                .data$Attribute_type, .data$Attribute_group
            ) |>
            dplyr::mutate(
                min = min(.data$Attribute_value_min),
                max = max(.data$Attribute_value_max)
            ) |>
            dplyr::ungroup() |>
            dplyr::mutate(
                Evidence = 'asr',
                Attribute_value = TRUE,
                Attribute_type = attr_type,
                Attribute_group = attr_group
                # Frequency = taxPPro::scores2Freq(Score)
            )
        select_cols <- c(
            'Parent_name', 'Parent_NCBI_ID', 'Parent_rank',
            'Attribute', 'min', 'max', 'Evidence',
            'Frequency', 'Score', 'Attribute_type', 'Attribute_group'
        )

        output <- output[,select_cols]
        colnames(output)[which(colnames(output) == 'Parent_name')] <- 'Taxon_name'
        colnames(output)[which(colnames(output) == 'Parent_NCBI_ID')] <- 'NCBI_ID'
        colnames(output)[which(colnames(output) == 'Parent_rank')] <- 'Rank'
        colnames(output)[which(colnames(output) == 'min')] <- 'Attribute_value_min'
        colnames(output)[which(colnames(output) == 'max')] <- 'Attribute_value_max'
        # output$Score <- NULL
        output <- unique(output)
        return(output)
    }

}


getDataReady <- function(x) {
    x$NCBI_ID[which(is.na(x$NCBI_ID))] <- 'unknown'
    x$Parent_NCBI_ID[which(is.na(x$Parent_NCBI_ID))] <- 'unknown'
    x <- x[x$Parent_NCBI_ID != 'unknown',]
    x <- x[!is.na(x$Rank), ]
    x <- x[!is.na(x$Evidence), ]
    x <- x[!is.na(x$Frequency), ]
    x <- x[!is.na(x$Confidence_in_curation), ]
    x <- unique(x)
    x <- freq2Scores(x)
    x_yesid <- x[which(x$NCBI_ID != 'unknown'),]
    x_noid <- x[which(x$NCBI_ID == 'unknown'),]
    x_noid_asr <- calcParentScores(x_noid)
    x_new <- dplyr::bind_rows(x_yesid, x_noid_asr)
    attr_type <- unique(x$Attribute_type)
    if (attr_type == 'logical') {
        cols <- c(
            'NCBI_ID', 'Rank',
            'Attribute', 'Attribute_value', 'Attribute_source',
            'Evidence', 'Frequency',
            'Attribute_type', 'Attribute_group',
            'Confidence_in_curation', 'Score'
        )
    } else if (attr_type == 'range') {
        cols <- c(
            'NCBI_ID', 'Rank',
            'Attribute', 'Attribute_value_min', 'Attribute_value_max',
            'Attribute_source',
            'Evidence', 'Frequency',
            'Attribute_type', 'Attribute_group',
            'Confidence_in_curation', 'Score'
        )
    }
    x_new <- unique(x_new[,cols])
    x_ready <- prepareData2(x_new)
    resolvedAgreements <- resolveAgreements(x_ready)
    resolvedConflicts <- resolveConflicts(resolvedAgreements)
    return(resolvedConflicts)
}








