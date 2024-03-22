
# Setup -------------------------------------------------------------------

library(purrr)

path2Files <- "../../" # Only for linux
fnames <- list.files(path = path2Files, pattern = ".csv", full.names = TRUE)

checkFnamesCSV <- function(x) {
    expected_fnames <- c(
        "bugphyzz_binary.csv",
        "bugphyzz_multistate.csv",
        "bugphyzz_numeric.csv"
    )
    if (!length(x))
        return(FALSE)
    fnames_ <- sub("^.*/bug", "bug", x)
    all(fnames_ == expected_fnames)
}

canReadCSV <- function(x) {
    if (!length(x))
        stop("Empty vector", call. = FALSE)
    purrr::map(fnames, ~ utils::read.csv(.x, skip = 1) )
}

expected_columns_multistate <- c(
    NCBI_ID = "integer", Taxon_name = "character",
    Rank = "character",
    Attribute = "character", Attribute_value = "character",
    Evidence = "character",
    Frequency = "character", Score = "double", Attribute_source = "character",
    Confidence_in_curation = "character", Attribute_type = "character"
)

expected_columns_binary <- c(
    NCBI_ID = "integer", Taxon_name = "character",
    Rank = "character",
    Attribute = "character", Attribute_value = "logical",
    Evidence = "character",
    Frequency = "character", Score = "double", Attribute_source = "character",
    Confidence_in_curation = "character", Attribute_type = "character"
)

expected_columns_numeric <- c(
    NCBI_ID = "integer", Taxon_name = "character",
    Rank = "character",
    Attribute = "character", Attribute_value = "double",
    Evidence = "character",
    Frequency = "character", Score = "double", Attribute_source = "character",
    Confidence_in_curation = "character", Attribute_type = "character",
    nsti = "double"
)

# Tests -------------------------------------------------------------------

test_that("All CSV files are present", {
    expect_true(checkFnamesCSV(fnames))
})

test_that("All CSV files can be read", {
    expect_true(length(fnames) > 0)
    expect_no_error({
        canReadCSV(fnames)
    })
})

dats <- map(fnames, ~ read.csv(.x, skip = 1))
names(dats) <- fnames

test_that("All columns are of right class and order", {
    ms <- dats[[grep("multistate", names(dats))]]
    expect_true(all(names(expected_columns_multistate) == colnames(ms)))
    expect_true(all(expected_columns_multistate == map_chr(ms, typeof)))

    bi <- dats[[grep("binary", names(dats))]]
    expect_true(all(names(expected_columns_binary) == colnames(bi)))
    expect_true(all(expected_columns_binary == map_chr(bi, typeof)))

    num <- dats[[grep("numeric", names(dats))]]
    expect_true(all(names(expected_columns_numeric) == colnames(num)))
    expect_true(all(expected_columns_numeric == map_chr(num, typeof)))
})
