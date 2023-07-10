
#' Remove Accession_ID and Genome_ID columns
#'
#' \code{removeAcccessionAndGenomeID} removes Accession_ID and Genome_ID from
#' a bugphyzz dataset. The reason is that right now these columns can be
#' incomplete or inconsistent in some datasets, or just missing in some others.
#' I think a solution would be to implement a relational database in which we
#' have a data object (data.frame?) with all of the taxids.
#'
#' @param df A data.frame imported from bugphyzz.
#'
#' @return A data.frame.
#' @export
#'
removeAccessionAndGenomeID <- function(df) {
    if ('Accession_ID' %in% colnames(df)) {
        df <- dplyr::select(df, -.data[['Accession_ID']])
    }
    if ('Genome_ID' %in% colnames(df)) {
        df <- dplyr::select(df, -.data[['Genome_ID']])
    }
    dplyr::distinct(df)
}

#' Homogenize attribute names of aerophilicity
#'
#' \code{homogenizeAerophilicityAttribbuteNames} makes all levels in the
#' aerophilicity dataset at the same level in the GO tree.
#'
#' @param df A data.frame.
#'
#' @return A data.frame.
#' @export
#'
homogenizeAerophilicityAttributeNames <- function(df) {
    df |> dplyr::mutate(
        Attribute = dplyr::case_when(
            Attribute == 'obligately anaerobic' ~ 'anaerobic',
            Attribute == 'microaerophilic' ~ 'aerobic',
            Attribute == 'obligately aerobic' ~ 'aerobic',
            TRUE ~ Attribute
        )
    )
}
