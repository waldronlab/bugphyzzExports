
library(bugphyzz)
library(purrr)

bp <- importBugphyzz(version = 'devel', force_download = TRUE)
ranks <- c('genus', 'strain', 'species', 'mixed')
tax_id_types <- c('Taxon_name', 'NCBI_ID')
header <- paste0("# bugphyzz ", Sys.Date(),
                 ", License: Creative Commons Attribution 4.0 International",
                 ", URL: https://waldronlab.io/bugphyzz")
# helper function to add a header line to an already written GMT file
addHeader <- function(header, out.file)
{
    fconn <- file(out.file, "r+")
    lines <- readLines(fconn)
    header <- sub("\n$", "", header)
    writeLines(c(header, lines), con = fconn)
    close(fconn)
}
l <- length(ranks) * length(tax_id_types)
sigs <- vector('list', l)
counter <- 1
for (i in seq_along(ranks)) {
    for (j in seq_along(tax_id_types)) {
        names(sigs)[counter] <- paste0(ranks[i], '--', tax_id_types[j])

        gmt_file <- paste0(
            'bugphyzz-', ranks[i], '-', tax_id_types[j], '-', '.gmt'
        )
        sig <- getBugphyzzSignatures(
            df = bp, tax.level = ranks[i], tax.id.type = tax_id_types[j]
        )
        bugsigdbr::writeGMT(sigs = sig, gmt.file = gmt_file)
        addHeader(header, gmt_file)
        counter <- counter + 1
    }
}
