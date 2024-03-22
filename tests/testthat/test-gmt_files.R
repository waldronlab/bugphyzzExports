
# Setup -------------------------------------------------------------------

library(clusterProfiler)
library(purrr)
path_to_files <- "../../" # Only for linux
checkFnamesGMT <- function(x) {
    if (!length(x))
        return(FALSE)
    expected_file_names <- c(
        "bugphyzz-genus-NCBI_ID.gmt", "bugphyzz-genus-Taxon_name.gmt",
        "bugphyzz-mixed-NCBI_ID.gmt",  "bugphyzz-mixed-Taxon_name.gmt",
        "bugphyzz-species-NCBI_ID.gmt", "bugphyzz-species-Taxon_name.gmt",
        "bugphyzz-strain-NCBI_ID.gmt", "bugphyzz-strain-Taxon_name.gmt"
    )
    fnames_ <- sub("^.*/bug", "bug", x)
    all(expected_file_names == fnames_)
}
canReadGMT <- function(x) {
    if (!length(x))
        stop("Empty vector", call. = FALSE)
    purrr::map(x, clusterProfiler::read.gmt)
}

# Tests -------------------------------------------------------------------

fnames <- list.files(path = path_to_files, pattern = "gmt", full.names = TRUE)

test_that("All gmt files are present", {
    expect_true(checkFnamesGMT(fnames))
})

test_that("All gmt files can be read", {
    expect_no_error({
        x <- canReadGMT(fnames)
    })
})

dats <- canReadGMT(fnames)
names(dats) <- fnames

test_that("NCBI_ID-GMT signatures can be converted to integers", {
    ncbi_dats <- dats[grep("NCBI_ID", names(dats))]
    are_int <- map_lgl(ncbi_dats, ~ {
        integ <- as.integer(.x$gene)
        any(is.na(integ))
    })
    expect_false(any(are_int))
})

test_that("Taxon_name-GMT signatures contain letters", {
    name_dats <- dats[grep("Taxon_name", names(dats))]
    are_ltrs <- map_lgl(name_dats, ~ {
       ltrs <- strsplit(.x$gene, "")[[1]] # just check the first element
        any(ltrs %in% letters) | any(ltrs %in% LETTERS)
    })
    expect_true(all(are_ltrs))
})
