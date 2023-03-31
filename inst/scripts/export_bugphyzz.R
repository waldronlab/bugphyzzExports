
# library(taxPPro)
# library(bugphyzz)
# library(bugsigdbr)
# library(readr)
# library(purrr)

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
        "very large" = list(lower = 6001, upper = 999999999999)
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
#' @importFrom bugphyzz makeSignatures
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
    signatures
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
#' @param header a character vector representing the header. Use default header
#' with date.
#'
#' @importFrom bugphyzz fattyAcidComposition physiologies
#'
#' @export
#'
#' @examples
#' makeAllSignatures()
makeAllSignatures <- function(header = .getHeader()) {
    ps <- physiologies()
    # fac <- fattyAcidComposition()
    # ps[["fatty acid composition"]] <- fac
    for (tax.level in .TAX.LEVELS) {
        for (tax.id.type in .TAX.ID.TYPES) {
            file_path <- file.path(tolower(paste0("bugphyzz-", tax.id.type, "-",
                                                  tax.level, ".gmt")))

            .writeFileWithHeader(
                .makeSignaturesByTaxIdAndLevel(ps, tax.id.type, tax.level),
                file_path,
                header)
        }
    }
    .makeFullDump("full_dump.csv", ps, header = header)
}

# makeAllSignatures()

library(bugphyzz)
library(taxPPro)
library(purrr)
library(rlang)
library(dplyr)
library(data.tree)

phys_names <- c('aerophilicity', 'growth temperature')
# phys_names <- 'all'
phys <- physiologies(phys_names, remove_false = TRUE, full_source = FALSE)

## For now, removing these datasets
phys[['metabolite utilization']] <- NULL ## Contains a mix of binary attributes
phys[['metabolite production']] <- NULL ## Contains a mix of binary attributes
phys[['isolation site']] <- NULL ## Too many values
phys[['habitat']] <- NULL ## Some FALSE values are important

## The Accession_ID or Genome_ID columns are missing from some datasets.
## Or these columns are incomplete or inconsistent across datasets.
phys <- phys |>
    map( ~ {
        if ('Accession_ID' %in% colnames(.x)) {
            .x <- select(.x, -Accession_ID)
        }
        if ('Genome_ID' %in% colnames(.x)) {
            .x <- select(.x, -Genome_ID)
        }
        distinct(.x)
    })

## Ensure that only valid attribute values are included.
fname <- system.file('extdata/attributes.tsv', package = 'bugphyzz')
attributes <- read.table(fname, header = TRUE, sep = '\t')
phys <- map(phys, ~ filter(.x, Attribute %in% unique(attributes$attribute)))
phys <- keep(phys, ~ nrow(.x) > 0)

## Prepare data in an uniform format before running propagation with
## functions from the taxPPro package (currently at sdgamboa/taxPPro)
data_ready <- vector('list', length(phys))
for (i in seq_along(data_ready)) {
    message('Preparing ', names(phys)[i], '.')
    names(data_ready)[i] <- names(phys)[i]
    data_ready[[i]] <- tryCatch(
        error = function(e) e,
        {
            prepareDatForPropagation(phys[[i]])
        }
    )
}
vct_lgl <- map_lgl(data_ready, is_error)
if (any(vct_lgl))  {
    message('Removing data with errors.')
    data_ready <- discard(data_ready, is_error)
}

## Run propagation with functions from the taxPPro package.
## Output is a data.tree R6 object.
data('tree_list')
tree <- as.Node(tree_list)
propagated <- vector('list', length(data_ready))
for (i in seq_along(propagated)) {
    message('Propagating ', names(data_ready)[i], ' - ', Sys.time())
    names(propagated)[i] <- names(data_ready)[i]
    propagated[[i]] <- tryCatch(
        error = function(e) e,
        {
            propagate(data_tree = tree, df = data_ready[[i]])
        }
    )
}

## Convert data.tree R6 to a data.frame
## NCBI taxonomy is appended. Could be useful for creating metaphlan-like
## names for the taxids
ncbi_taxonomy <- get_ncbi_taxonomy()
dfs <- map(propagated, ~ toDataFrame(.x, ncbi_tax = ncbi_taxonomy))

## Append attribute group and attribute type information.
## This information was lost during the propagation step.
for (i in seq_along(dfs)) {
    if (names(dfs)[i] %in% names(phys)) {
        dfs[[i]]$Attribute_group <- unique(phys[[names(dfs[i])]]$Attribute_group)
        dfs[[i]]$Attribute_type <- unique(phys[[names(dfs[i])]]$Attribute_type)
    }
}

## Mix Attribute group and Attribute (just attributes with logical type)
data_ready <- data_ready |>
    map(~ {
        attr_type <- unique(.x$Attribute_type)
        if (attr_type == 'logical')  {
            .x$Attribute <- paste0(.x$Attribute_group, ':', .x$Attribute)
        }
        .x
    })

output <- vector('list', length(dfs))
for (i in seq_along(output)) {
    data <- dfs[[i]]
    names(output)[i] <- names(dfs)[i]
    attr_type <- unique(data$Attribute_type)
    if (attr_type == 'logical') {
        common_names <- grep('__', names(data), value = TRUE, invert = TRUE)
        unique_names <- grep('__', names(data), value = TRUE)
        attr_names <- unique(sub('__.*$', '', unique_names))
        data <- map(attr_names, ~ {
            names <- c(grep(.x, unique_names, value = TRUE), common_names)
            k <- data[,names]
            names(k) <- sub('.*__', '', names(k))
            k[['Attribute']] <- .x
            k[['Attribute_value']] <- TRUE
            cols <- c(
                'NCBI_ID', 'Taxon_name', 'Attribute', 'Attribute_value',
                'Evidence', 'Score', 'Rank'
            )
            k <- k[,cols]
            k <- k[which(!is.na(k$Evidence)),]
            k <- k[which(k$Evidence %in% c('asr', 'inh')),]
            k$Frequency <- scores2Freq(k$Score)
            k <- k[which(k$Frequency != 'unknown'),]
            k
        })
        names(data) <- attr_names
        data <- bind_rows(data)
        attr_grp <- unique(data_ready[[names(dfs)[i]]][['Attribute_group']])
        data[['Attribute_group']] <- attr_grp
        data[['Attribute_type']] <- 'logical'
        data <- bind_rows(data_ready[[names(dfs)[i]]], data)
        output[[i]] <- data
    } else if (attr_type == 'range') {
        cols <- c(
            'NCBI_ID', 'Taxon_name', 'Attribute',
            'Attribute_value_min', 'Attribute_value_max',
            'Evidence', 'Score', 'Rank'
        )
        attr_name <- unique(sub('__.*$', '', grep('__', names(data), value = TRUE)))
        colnames(data) <- sub('.*__', '', colnames(data))
        data[['Attribute']] <- sub('_', ' ', attr_name)
        data <- data[,cols]
        data <- data[which(!is.na(data$Evidence)),]
        data <- data[which(data$Evidence %in% c('asr', 'inh')),]
        data$Frequency <- scores2Freq(data$Score)
        data <- data[which(data$Frequency != 'unknown'),]
        attr_grp <- unique(data_ready[[names(dfs)[i]]][['Attribute_group']])
        data[['Attribute_group']] <- attr_grp
        data[['Attribute_type']] <- 'range'
        data <- bind_rows(data_ready[[names(dfs)[i]]], data)
        output[[i]] <- data
    }
}

## Add code here for automatically selecting thresholds, filtering, and
## changing numeric attributes to categorical/logical
test_output <- vector('list', length(output))
names(test_output) <- names(output)
for (i in seq_along(test_output)) {
    if (.hasSpecialThresholds(names(output)[i])) {
        thresholds <- .getSpecialThresholds(names(output[i]))
        if (!is.null(thresholds)) {
            test_output[[i]] <- thresholds
        }
    }
}

th <- discard(test_output, is.null)


map(th, ~ {
    names(.x)
})


subsetThreshold <- function(df, thr) {
    lower <- thr$lower
    upper <- thr$upper
    if (is.na(lower)) {
        df <- df[which(df$Attribute_value_max < upper),]
    } else if (is.na(upper)) {
        df <- df[which(df$Attribute_value_min >= lower),]
    } else {
        df <- df[which(df$Attribute_value_min >= lower * df$Attribute_value_max < upper),]
    }
}




gt <- output$`growth temperature`
gt <- gt |>
    mutate(
        Attribute = case_when(
            Attribute_value_max < 25 ~ 'psychrophile',
            Attribute_value_min >= 25 & Attribute_value_max < 45 ~ 'mesophile',
            Attribute_value_min >= 45 & Attribute_value_max < 60 ~ 'thermophile',
            Attribute_value_min >= 60 ~ 'hyperthermophile',
            TRUE ~ NA
        ),
        Attribute = paste0(Attribute_group, ':', Attribute)
    )




output$`growth temperature` <- gt
full_dump <- reduce(output, bind_rows)
full_dump <- full_dump |>
    mutate(
        Rank = case_when(
            grepl('t__', NCBI_ID) ~ 'strain',
            grepl('s__', NCBI_ID) ~ 'species',
            grepl('g__', NCBI_ID) ~ 'genus',
            grepl('f__', NCBI_ID) ~ 'family',
            grepl('o__', NCBI_ID) ~ 'order',
            grepl('c__', NCBI_ID) ~ 'class',
            grepl('p__', NCBI_ID) ~ 'phylum',
            grepl('d__', NCBI_ID) ~ 'domain',
            TRUE ~ Rank),
        Attribute = sub(' ', '_', Attribute)
    )

## Create and export dump file(s)
## A file with numeric values
## TODO

## A file with categorical labels for numeric values
full_dump$NCBI_ID <- sub('^[dpcofgst]__', '', full_dump$NCBI_ID)
fname <- paste0("full_dump_bugphyzz.csv.bz2")
unlink(fname)
con <- bzfile(fname, "w")
write.csv(full_dump, file = con, quote = TRUE, row.names = FALSE)
close(con)

## Code for creating signatures

## TODO
