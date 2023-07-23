
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

#' Check rank
#'
#' \code{checkRank} checks that the rank matches the NCBI ID (taxids).
#' If not, NA is returned.
#'
#' @param taxid A character string with a taxid.
#'
#' @return A character string. NA if error.
#' @export
#'
checkRank <- function(taxid) {
    tryCatch(
        error = function(e) NA,
        {
            taxizedb::taxid2rank(taxid, db = 'ncbi')
        }
    )
}

#' Check taxon name
#'
#' \code{checkTaxonName} checks that the taxon_name corresponds to the NCBI ID (taxids).
#' If not, NA is returned.
#'
#' @param taxid A character string with a taxid.
#'
#' @return A character string. NA if error.
#' @export
#'
checkTaxonName <- function(taxids) {
    tryCatch(
        error = function(e) NA,
        {
            taxizedb::taxid2name(taxid, db = 'ncbi')
        }
    )
}
