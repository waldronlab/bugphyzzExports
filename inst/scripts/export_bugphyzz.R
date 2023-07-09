## Script to create bugphyzz exports dump files and signature files

library(bugphyzz)
library(taxPPro)
library(purrr)
library(rlang)
library(dplyr)
library(data.tree)
library(bugphyzzExports)

# library(bugsigdbr)
# library(readr)

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
        "hyperthermophile" = list(lower = 61, upper = NA)
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
        "very alkaline" = list(lower = 9.76, upper = NA)
    ),
    "width" = list(
        "small" = list(lower = NA, upper = 0.9),
        "average" = list(lower = 0.91, upper = 3.5),
        "large" = list(lower = 3.51, upper = 12),
        "very large" = list(lower = 13, upper = NA)
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


# Code starts -------------------------------------------------------------

message('>>>>>>> Importing data ', Sys.time(), ' <<<<<<')

phys_names <- c(
    ## Categorical
    'aerophilicity',
    'gram stain',
    # 'shape',
    # 'acetate producing',
    # 'animal pathogen',
    # 'arrangement',
    # 'biofilm forming',
    # 'biosafety level',
    # 'butyrate producing',
    # 'disease association',
    # 'halophily',

    ## Numeric
    'growth temperature'
    # 'optimal ph',
    # 'width'
)
# phys_names <- 'all'
phys <- physiologies(phys_names, remove_false = TRUE, full_source = FALSE)

## The Accession_ID or Genome_ID columns are missing from some datasets.
## These columns can be incomplete or inconsistent in some datasets.
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

## Make sure that only valid attribute values are included.
fname <- system.file('extdata/attributes.tsv', package = 'bugphyzz')
attributes <- read.table(fname, header = TRUE, sep = '\t')
phys <- map(phys, ~ filter(.x, Attribute %in% unique(attributes$attribute)))
phys <- keep(phys, ~ nrow(.x) > 0)

## This code is to make sure all annotations in the aerophilicity dataset
## are at the same level in the GO tree
phys$aerophilicity <- phys$aerophilicity |>
    mutate(
        Attribute = case_when(
            Attribute == 'obligately anaerobic' ~ 'anaerobic',
            Attribute == 'microaerophilic' ~ 'aerobic',
            Attribute == 'obligately aerobic' ~ 'aerobic',
            TRUE ~ Attribute
        )
    )

message('>>>>>>> Preparing data ', Sys.time(), ' <<<<<<')
## Prepare data in an uniform format before running propagation.
## Functions from the taxPPro package (currently at sdgamboa/taxPPro)
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
data_ready <- discard(data_ready, is_error)

message('>>>>>>> Propagating data ', Sys.time(), ' <<<<<<')
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
propagated <- discard(propagated, is_error)

## Convert data.tree R6 to a data.frame
## NCBI taxonomy is added. Could be useful for creating metaphlan-like
## names for the taxids
ncbi_taxonomy <- get_ncbi_taxonomy()
dfs <- map(propagated, ~ toDataFrame(.x, ncbi_tax = ncbi_taxonomy))

## Add attribute group and attribute type information to the data.frames.
## This information was lost during the propagation step, but is needed to treat
## differently logical/categorical attributes and numeric attributes.
for (i in seq_along(dfs)) {
    if (names(dfs)[i] %in% names(phys)) {
        dfs[[i]]$Attribute_group <- unique(phys[[names(dfs[i])]]$Attribute_group)
        dfs[[i]]$Attribute_type <- unique(phys[[names(dfs[i])]]$Attribute_type)
    }
}

## Concatenate the name of the Attribute group with the Attribute
## (just attributes with logical type), e.g. aerophilicity:aerobic
data_ready <- data_ready |>
    map(~ {
        attr_type <- unique(.x$Attribute_type)
        if (attr_type == 'logical')  {
            .x$Attribute <- paste0(.x$Attribute_group, ':', .x$Attribute)
        }
        .x
    })

message('>>>>>>> Creaing output data ', Sys.time(), ' <<<<<<')

## Some final edits to the output data.frames
## In the code below, '__' was used to mark the attributes specific
## to a certain attribute. Example: aerobic__Attribute_source
output <- vector('list', length(dfs))
for (i in seq_along(output)) {
    data <- dfs[[i]]
    names(output)[i] <- names(dfs)[i]
    attr_type <- unique(data$Attribute_type)
    if (attr_type == 'logical') {
        ## Here, attribute names are obtained form the columns.
        ## Attribute names in colums was the format in the propagation step.
        common_names <- grep('__', colnames(data), value = TRUE, invert = TRUE)
        unique_names <- grep('__', colnames(data), value = TRUE)
        attr_names <- unique(sub('__.*$', '', unique_names))
        ## This map function call repeat the process for each attribute
        ## e.g. aerophilicity:aerobic.
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
            ## Get only asr and inh which is the main output from the
            ## propagation step
            k <- k[which(k$Evidence %in% c('asr', 'inh')),]
            k$Frequency <- scores2Freq(k$Score)
            # k <- k[which(k$Frequency != 'unknown'),]
            k
        })
        names(data) <- attr_names
        data <- bind_rows(data)
        ## The attribute group is the same for all attributes.
        ## e.g. aerophilicity for aerophilicity:aerobic and aerophiliciyt:anaerobic
        attr_grp <- unique(data_ready[[names(dfs)[i]]][['Attribute_group']])
        data[['Attribute_group']] <- attr_grp
        data[['Attribute_type']] <- 'logical'
        ## Here, we join the data before the propagation (inh, exp, tas, nas
        ## with the data after the propagation (asr, inh)
        data <- bind_rows(data_ready[[names(dfs)[i]]], data)
        output[[i]] <- data
    } else if (attr_type == 'range') {
        cols <- c(
            'NCBI_ID', 'Taxon_name', 'Attribute',
            'Attribute_value_min', 'Attribute_value_max',
            'Evidence', 'Score', 'Rank'
        )
        ## In this case, the attribute is the same for all data
        attr_name <- unique(sub('__.*$', '', grep('__', names(data), value = TRUE)))
        colnames(data) <- sub('.*__', '', colnames(data))
        data[['Attribute']] <- sub('_', ' ', attr_name)
        data <- data[,cols]
        data <- data[which(!is.na(data$Evidence)),]
        ## Filter asr and inh from the propagation output
        data <- data[which(data$Evidence %in% c('asr', 'inh')),]
        data$Frequency <- scores2Freq(data$Score)
        # data <- data[which(data$Frequency != 'unknown'),]
        attr_grp <- unique(data_ready[[names(dfs)[i]]][['Attribute_group']])
        data[['Attribute_group']] <- attr_grp
        data[['Attribute_type']] <- 'range'
        data <- bind_rows(data_ready[[names(dfs)[i]]], data)
        output[[i]] <- data
    }
}


message('>>>>>>> Exporting first dump file ', Sys.time(), ' <<<<<<')

## Code for a first dump file, containing only numeric attributes.
# numeric_attributes <- keep(output, ~ unique(.x$Attribute_type) == 'range')
numeric_attributes <- keep(output, ~ unique(.x$Attribute_type) == 'range')
full_dump_numeric <- reduce(numeric_attributes, bind_rows)
full_dump_numeric$NCBI_ID <- sub('^[dpcofgst]__', '', full_dump_numeric$NCBI_ID)
full_dump_numeric$Attribute <- gsub(' ', '_', full_dump_numeric$Attribute)
fname <- paste0("full_dump_numeric.csv.bz2")
unlink(fname)
con <- bzfile(fname, "w")
write.csv(full_dump_numeric, file = con, quote = TRUE, row.names = FALSE)
close(con)

## Add code here for automatically selecting thresholds, filtering, and
## changing numeric attributes to categorical/logical
message('>>>>>>> Exporting second dump file ', Sys.time(), ' <<<<<<')

categorical_attributes <- keep(output, ~ unique(.x$Attribute_type) == 'logical')
numeric_attributes_with_thr <- keep(
    .x = numeric_attributes,
    .p = ~ .hasSpecialThresholds(unique(.x$Attribute_group))
)
new_categorical_attributes <- map(
    .x = numeric_attributes_with_thr,
    .f =  ~ {
        thresholds <-.getSpecialThresholds(unique(.x$Attribute_group))
        rangeToLogicalThr(.x, thresholds)
    }
)
cat_list <- c(categorical_attributes, new_categorical_attributes)
full_dump_cat <- reduce(cat_list, bind_rows)
full_dump_cat$NCBI_ID <- sub('^[dpcofgst]__', '', full_dump_cat$NCBI_ID)
full_dump_cat$Attribute <- gsub(' ', '_', full_dump_cat$Attribute)
fname2 <- paste0("full_dump_categorical.csv.bz2")
unlink(fname2)
con <- bzfile(fname2, "w")
write.csv(full_dump_cat, file = con, quote = TRUE, row.names = FALSE)
close(con)

message('>>>>>>> Finished exporting dump files ', Sys.time(), ' <<<<<<')

## Code for creating signatures and exporting signatures
## Add headers
## TODO
