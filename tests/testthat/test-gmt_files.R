library(qusage)
expected_file_names <- c(
        "bugphyzz-genus-NCBI_ID.gmt", "bugphyzz-genus-Taxon_name.gmt",
        "bugphyzz-species-NCBI_ID.gmt", "bugphyzz-species-Taxon_name.gmt",
        "bugphyzz-strain-NCBI_ID.gmt", "bugphyzz-strain-Taxon_name.gmt",
        "bugphyzz-mixed-NCBI_ID.gmt", "bugphyzz-mixed-Taxon_name.gmt"
    )
fnames <- list.files(path = "../", pattern = "gmt")
sigs <- lapply(fnames, read.gmt)

test_that("All gmt files are present", {
    res <- all(fnames %in% expected_file_names)
    expect_true(res)
})

test_that("All gmt files can be read", {
    expect_no_error({
        gmt_files <- lapply(fnames, read.gmt)
    })
    expect_no_warning({
        gmt_files <- lapply(fnames, read.gmt)
    })
})
