library(qusage)
library(purrr)
expected_file_names <- c(
        "bugphyzz-genus-NCBI_ID.gmt", "bugphyzz-genus-Taxon_name.gmt",
        "bugphyzz-species-NCBI_ID.gmt", "bugphyzz-species-Taxon_name.gmt",
        "bugphyzz-strain-NCBI_ID.gmt", "bugphyzz-strain-Taxon_name.gmt",
        "bugphyzz-mixed-NCBI_ID.gmt", "bugphyzz-mixed-Taxon_name.gmt"
    )
fnames <- list.files(path = "../", pattern = "gmt")
sigs <- map(fnames, read.gmt)
names(sigs) <- fnames
## First line is a comment that should be removed
sigs <- map(sigs, function(x) discard(x, ~ !length(.x)))

test_that("All gmt files are present", {
    res <- all(fnames %in% expected_file_names)
    expect_true(res)
})

test_that("All gmt files can be read", {
    expect_no_error({
        gmt_files <- map(fnames, read.gmt)
    })
    expect_no_warning({
        gmt_files <- map(fnames, read.gmt)
    })
})

test_that("NCBI_ID-GMT signatures can be converted to integers", {
    ncbi_sigs <- sigs[grep("NCBI_ID", names(sigs))]
    ncbi_sigs <- list_flatten(ncbi_sigs)
    are_int <- map_lgl(ncbi_sigs, ~ {
        integ <- as.integer(.x)
        any(is.na(integ))
    })
    expect_false(any(are_int))
})

test_that("Taxon_name-GMT signatures contain letters", {
    name_sigs <- sigs[grep("Taxon_name", names(sigs))]
    name_sigs <- list_flatten(name_sigs)
    are_ltrs <- map_lgl(name_sigs, ~ {
       ltrs <- strsplit(.x, "")[[1]]
        any(ltrs %in% letters) | any(ltrs %in% LETTERS)
    })
    expect_true(all(are_ltrs))
})
