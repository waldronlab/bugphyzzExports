library(bugphyzz)
library(bugsigdbr)
library(readr)

# Lower bounds are excluded but upper bounds are included

.THRESHOLDS <- list(
    "acetate producing" = list(
        "small" = list(lower = NA, upper = 0.7),
        "moderate" = list(lower = 0.71, upper = 7.2),
        "large" = list(lower = 7.21, upper = 47.9)
    ), 
    "butyrate producing" = list(
        "low" = list(lower = NA, upper = 4.9),
        "medium" = list(lower = 4.91, upper = 12.7),
        "high" = list(lower = 12.71, upper = 21.3)
    ),
    "coding genes" = list(
        "very small" = list(lower = NA, upper = 473),
        "small" = list(lower = 474, upper = 600),
        "average" = list(lower = 601, upper = 6000), 
        "very large"= list(lower = 6001, upper = 999999999999)
    ),
    "cogem pathogenicity rating" = list(
        "does not belong to a species of which representatives are known to be pathogenic for humans, animals or plants" = list(lower = NA, upper = 1),
        "cause a disease in humans or animals" = list(lower = 1, upper = 2),
        "serious disease in humans or animals" = list(lower = 3, upper = 3)
    ),
    "fatty acid" = list(
        "minimally present" = list(lower = NA, upper = 0.99),
        "lower concentration" = list(lower = 0.991, upper = 23.2),
        "concentrated" = list(lower = 23.21, upper = 55.5),
        "very concentrated" = list(lower = 55.51, upper = 76.6)
    ),
    "genome size" = list(
        "small" = list(lower = NA, upper = 490885),
        "average" = list(lower = 490886, upper = 998123),
        "large" = list(lower = 998124, upper = 6997434),
        "very large" = list(lower = 6997435, upper = 16040666)
    ),
    "growth temperature" = list(
        "psychrophile" = list(lower = NA, upper = 24.9),
        "mesophile" = list(lower = 25, upper = 45),
        "thermophile" = list(lower = 46, upper = 60),
        "hyperthermophile" = list(lower = 61, upper = 121)
    ),
    "hydrogen gas producing" = list(
        "low" = list(lower = NA, upper = 1.7),
        "medium" = list(lower = 1.6, upper = 4.3),
        "high" = list(lower = 4.4, upper = 14.8)
    ),
    "lactate producing" = list(
        "low" = list(lower = NA, upper = 2.1),
        "medium" = list(lower = 2.1, upper = 7.1),
        "high" = list(lower = 7.1, upper = 26.5)
    ),
    "length" = list(
        "small" = list(lower = NA, upper = 3.8),
        "average" = list(lower = 3.9, upper = 22),
        "large" = list(lower = 23, upper = 60),
        "very large" = list(lower = 61, upper = 250)
    ),
    "mutation rate per site per generation" = list(
        "slow" = list(lower = NA, upper = 2.92),
        "medium" = list(lower = 2.93, upper = 16),
        "fast" = list(lower = 17, upper = 97.8)
    ),
    "mutation rate per site per year" = list(
        "slow" = list(lower = NA, upper = 7.5),
        "medium" = list(lower = 7.6, upper = 20),
        "medium fast" = list(lower = 21, upper = 54.2),
        "fast" = list(lower = 54.3, upper = 410)
    ),
    "optimal ph" = list(
        "acidic" = list(lower = NA, upper = 5),
        "neutral" = list(lower = 6, upper = 7),
        "alkaline" = list(lower = 8, upper = 9.75),
        "very alkaline" = list(lower = 9.76, upper = 11.25)
    ),
    "width" = list(
        "small" = list(lower = NA, upper = 0.9),
        "average" = list(lower = 0.91, upper = 3.5),
        "large" = list(lower = 3.51, upper = 12),
        "very large" = list(lower = 13, upper = 100)
    )
)

.TAX.LEVELS <- c("species", "genus")

.TAX.ID.TYPES <- c("Taxon_name", "NCBI_ID")

#' Returns if a physiology has manually curated thresholds
#'
#' @param physiology the name of a physiology
#' 
#' @return boolean
#'
#' @examples
#' .hasSpecialThreshold()
.hasSpecialThresholds <- function(physiology) {
    physiology %in% names(.THRESHOLDS)
}

#' Returns a physiology's manually curated thresholds
#'
#' @param physiology the name of a physiology
#' 
#' @return list of lists representing threshold ranges
#'
#' @examples
#' .getSpecialThreshold()
.getSpecialThresholds <- function(physiology) {
    .THRESHOLDS[[physiology]]
}

#' Get list of physiology signatures by attributes and attribute values
#' 
#' @param physiology data frame from bugphyzz
#' @param tax.id.type "Taxon_name" or "NCBI_ID"
#' @param tax.level "species" or "genus"
#' @return list of signature lists
#' 
#' @importFrom dplyr filter
#' @importFrom bugphyzz getSignatures
#' 
#' @example
#' ps <- physiologies()
#' p <- ps$`animal pathogen`
#' .getSignaturesByThreshold(p, "species", "Taxon_name")
.getSignaturesByThreshold <- function(physiology,
                                      tax.id.type = .TAX.ID.TYPES,
                                      tax.level = .TAX.LEVELS,
                                      lower.bound = NA,
                                      upper.bound = NA) {
    if (is.na(lower.bound) && is.na(upper.bound))
        makeSignatures(physiology, tax.id.type, tax.level)
    else if (is.na(lower.bound))
        makeSignatures(physiology, tax.id.type, tax.level, max = upper.bound)
    else if (is.na(upper.bound))
        makeSignatures(physiology, tax.id.type, tax.level, min = lower.bound)
    else {
        makeSignatures(physiology, tax.id.type, tax.level, min = lower.bound,
                       max = upper.bound)
    }
}

#' Make bugphyzz signatures
#'
#' @param ps list of data frames of physiologies from bugphyzz
#' @param tax.id.type a value in .TAX.ID.TYPES
#' @param tax.level a value in .TAX.LEVELS
#' @return a list of data frames representing signatures
#'
#' @importFrom bugphyzz fattyAcidComposition physiologies physiologiesList
#'             makeSignatures
#' @importFrom dplyr filter select
#'
#' @export
#'
#' @examples
#' .makeSignaturesByTaxIdAndLevel(fattyAcidComposition(), physiologies(), 
#'                               "species", "NCBI_ID")
.makeSignaturesByTaxIdAndLevel <- function(ps = physiologies(),
                                          tax.id.type = .TAX.ID.TYPES,
                                          tax.level = .TAX.LEVELS) {
    signatures <- list()
    for (p in names(ps)) {
        if (.hasSpecialThresholds(p)) {
            thresholds <- .getSpecialThresholds(p)
            for (threshold in names(thresholds)) {
                lower_bound <- thresholds[[threshold]]$lower
                upper_bound <- thresholds[[threshold]]$upper  
                signatures <- append(signatures, 
                                     .getSignaturesByThreshold(p,
                                                               tax.id.type,
                                                               tax.level,
                                                               lower_bound,
                                                               upper_bound))
            }
        } else {
            signatures <- append(signatures, makeSignatures(p,
                                                            tax.id.type,
                                                            tax.level))
        }
    }
    signatures <- append(signatures, signatures)
}

#' Header for bugphyzz files
#' 
#' @param identifier any character vector. Defaults to today.
#' 
#' @examples
#' .getHeader("3.17")
.getHeader <- function(identifier = format(Sys.time(), "%Y-%m-%d")) {
    paste0("# bugphyzz ", identifier,
           ", License: Creative Commons Attribution 4.0 International",
           ", URL: https://waldronlab.io/bugphyzz")
}

#' Write a header for a file given its path
#'
#' @param file_path path to the file
#' @param header a character vector representing the header
#' 
#' @importFrom readr read_lines write_lines
#'
#' @examples
#' writeHeader(file.path(tempdir(), "test.txt"))
.writeHeader <- function(file_path, header = header) {
    lines <- read_lines(file_path)
    write_lines(c(header, lines), file_path)
}

#' Write file given a signature with header
#' 
#' @param signatures to write to file
#' @param file_path to write a file
#' @param header a character vector representing the header
#'
#' @importFrom bugsigdbr writeGMT
#' @importFrom readr write_csv
#' 
#' @examples
#' writeFileWithHeader(signatures, "mysigs.gmt", .getHeader())
.writeFileWithHeader <- function(signatures, file_path, header = header) {
    writeGMT(signatures, file_path)
    .writeHeader(file_path, header)
}

#' Write a full dump csv
#' 
#' Include only similar columns
#'
#' @param file_path to write a file
#' @param ps physiologies data frame from bugphyzz
#' @param header a character vector representing the header
#'
#' @importFrom bugphyzz fattyAcidComposition physiologies
#' @importFrom readr write_csv
#'
#' @examples
#' makeFullDump(fattyAcidComposition(), physiologies())
.makeFullDump <- function(file_path, ps = physiologies(), header = header) {
    colname_intersection <- Reduce(intersect, lapply(ps, colnames))
    all_ps_ci <- lapply(ps, function(p) p[colname_intersection])
    for (i in 1:length(all_ps_ci)) {
        write_csv(all_ps_ci[[i]], file_path, col_names = i == 1, append = i != 1)
    }    
    .writeHeader(file_path, header)
}

#' Make all signatures
#' 
#' Write four files for the combinations of .TAX.LEVELS and .TAX.ID.TYPES
#'
#' @param header a character vector representing the header. Use default header with date.
#' 
#' @importFrom bugphyzz fattyAcidComposition physiologies
#' 
#' @export
#' 
#' @examples
#' makeAllSignatures()
makeAllSignatures <- function(header = .getHeader()) {
    signatures <- list()
    fac <- fattyAcidComposition()
    ps <- physiologies()
    all_ps <- append(ps, list(`fatty acid composition` = fac))
    for (tax.level in .TAX.LEVELS) {
        for (tax.id.type in .TAX.ID.TYPES) {
            file_path <- file.path(tolower(paste0("bugphyzz-", tax.id.type, "-",
                                                  tax.level, ".gmt")))
            signatures <- append(signatures,
                                 .makeSignaturesByTaxIdAndLevel(ps,
                                                                tax.id.type,
                                                                tax.level))
            .writeFileWithHeader(signatures, file_path, header)
        }
    }
    .makeFullDump("full_dump.csv", all_ps, header = header)
}

makeAllSignatures()
