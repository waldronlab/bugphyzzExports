## Script to create bugphyzz exports dump files and signature files

library(logr)
library(bugphyzz) # BiocManager::install('waldronlab/bugphyzz', force = TRUE)
library(taxPPro) # BiocManager::install('sdgamboa/taxPPro', force = TRUE)
library(purrr)
library(rlang)
library(dplyr)
library(data.tree)
library(bugphyzzExports) # Install locally
library(BiocParallel)
library(tidyr)

logfile <- "log_file"
lf <- log_open(logfile, logdir = FALSE, compact = TRUE, show_notes = FALSE)

n_threads <- parallel::detectCores()
if (n_threads > 16) {
    n_threads <- round(n_threads * 0.6)
}
msg <- paste0('Using ', n_threads, ' cores.')
log_print(msg, blank_after = TRUE)

## For now, do not inlcude these physiologies/Attribute groups.
## They need more curation.
exlclude_phys <- c(
    'country', 'geographic location',
    'habitat', 'isolation site',
    'metabolite production', 'metabolite utilization',
    'halophily'
)

## These are true binary attributes.
## We need to keep TRUE and FALSE values and merge them with the
## attribute name
binaries <- c(
    "acetate producing",
    "animal pathogen",
    "antimicrobial sensitivity",
    "biofilm forming",
    "butyrate producing",
    "extreme environment",
    "health associated",
    "hydrogen gas producing",
    "lactate producing",
    "motility",
    "pathogenicity human",
    "plant pathogenicity",
    "sphingolipid producing",
    "spore formation"
)

phys_names <- showPhys()
phys_names <- phys_names[which(!phys_names %in% exlclude_phys)]

msg <- paste0('"', paste0(phys_names, collapse = ', '), '"')
msg_len <- length(phys_names)
msg <- paste('Importing', msg_len, 'physiologies from bugphyzz:', msg, '--', Sys.time())
log_print(msg, blank_after = TRUE)

phys <- physiologies(phys_names, full_source = FALSE)
phys <- map(phys, ~ {
    attr_grp <- unique(.x$Attribute_group)
    if (attr_grp %in% binaries) {
        .x$Attribute <- paste0(.x$Attribute, '--', .x$Attribute_value)
        .x$Attribute_value <- TRUE
    }
    if (unique(.x$Attribute_type) == 'logical') {
        .x <- filter(.x, Attribute_value == TRUE)
    }
    if ('Unit' %in% colnames(.x)) {
        unit <- .x$Unit
        unit <- unique(unit[!is.na(unit)])
        .x$Unit <- unit
    }
    return(.x)
})

## Convert all attributes to categorical before propagation
## This would allow to use the same AUR-ROC method for comparison
categorical <- keep(phys, ~ unique(.x$Attribute_type) == 'logical')
range <- keep(phys, ~ unique(.x$Attribute_type == 'range'))
range <- range[which(names(range) %in% names(THRESHOLDS()))]
range_cat <- map2(range, names(range), ~ rangeToLogicalThr(.x, THRESHOLDS()[[.y]]))
categorical <- c(categorical, range_cat)


## The next chunk checks that only attributes with valid values are included.
## Those attributes with invalid values are reported in the log file.
## <Attribute>--TRUE and <Attribute>--FALSE are added for binary attributes.
log_print('Check that all attributes are valid. Invalid values will be printed and dropped from the full dump file:', blank_after = TRUE)
fname <- system.file('extdata/attributes.tsv', package = 'bugphyzz')
valid_attributes <- unique(read.table(fname, header = TRUE, sep = '\t')$attribute)
more_valid_attributes <- map(categorical, ~ {
    attr_grp <- unique(.x$Attribute_group)
    if (attr_grp %in% binaries) {
        binary_attr <- unique(.x$Attribute)
        return(binary_attr)
    }
}) |>
    discard(is.null) |>
    unlist(recursive = TRUE, use.names = FALSE)
valid_attributes <- unique(c(valid_attributes, more_valid_attributes))

data <- map(categorical, ~ {
    attr_names <- unique(.x$Attribute)
    attr_grp <- unique(.x$Attribute_group)
    lgl <- sum(!attr_names %in% valid_attributes)
    if (lgl > 0) {
        invalid_values <- filter(.x, !Attribute %in% valid_attributes)
        invalid_values <- invalid_values |>
            select(Attribute_group, Attribute) |>
            unique() |>
            as_tibble()
        log_print(paste0('Invalid values for ', attr_grp, ': '))
        log_print(invalid_values, blank_after = TRUE)
    }
    output <- filter(.x, Attribute %in% valid_attributes)
    return(output)
}) |>
    discard(~ !nrow(.x))
data_discarded <- names(categorical)[which(!names(categorical) %in% names(data))]
if (length(data_discarded) > 0) {
    data_discarded <- paste0(', paste0(data_discarded, collapse = ', '), ')
    msg <- paste0(
        "The following physiologies were discarded because they didn't have",
        " any valid Attribute: ", data_discarded, '.'
    )
log_print(msg, blank_after = TRUE)
}

## This chunk is to make sure that al NCBI_IDs have valid taxon names
data <- bplapply(data, BPPARAM = MulticoreParam(workers = n_threads), FUN = function(x) {
    set1 <- filter(x, !is.na(NCBI_ID))
    set2 <- filter(x, is.na(NCBI_ID))
    set1$Rank <- checkRank(set1$NCBI_ID)
    set1$Taxon_name <- checkTaxonName(set1$NCBI_ID)
    bind_rows(set1, set2)
})

## For those numeric attributes converted to categorical attributes,
## propagate them with their ranges and units.
data <- map(data, ~ {
    if ('Attribute_range' %in% colnames(.x)) {
        .x$Attribute <- paste0(.x$Attribute, ' ', .x$Attribute_range)
        .x$Attribute <- sub('\\)$', '', .x$Attribute)
        .x$Attribute <- paste0(.x$Attribute, ' ', .x$Unit, ')')
        .x$Attribute <- sub(' \\)', ')', .x$Attribute)
    }
    return(.x)
})


data_ready <- bplapply(data, BPPARAM = MulticoreParam(workers = n_threads), FUN = function(x) {
    tryCatch(
        error = function(e) e,
        {
            output <- x |>
                prepareDataForPropagation() |>
                mergeOriginalAndEarlyASR() |>
                group_by(NCBI_ID) |>
                mutate(Score = Score / sum(Score)) |>
                ungroup() |>
                distinct()
            return(output)
        }
    )
})

lgl <- map_lgl(data_ready, is_error)
if (any(lgl)) {
    discard_names <- names(data_ready)[which(lgl)]
    msg <- paste0(
        "The following physiologies will be discared for some error during",
        " propagation: ",
        paste0("'", paste0(discard_names, collapse = ', '), ".'")
    )
    log_print(msg, blank_after = TRUE)
    data_ready <- discard(data_ready, is_error)
}

## Propagation with taxPPro
data('tree_list')
tree <- as.Node(tree_list)
propagated <- bplapply(X = data_ready, BPPARAM = MulticoreParam(workers = n_threads), FUN = function(x) {
    msg <- unique(x$Attribute_group)
    msg <- paste0('Propagating ', msg, '...')
    log_print(msg)
    input_tbl <- x |>
        select(NCBI_ID, Attribute, Score, Evidence) |>
        distinct() |>
        complete(NCBI_ID, Attribute, fill = list(Score = 0, Evidence = ''))
    l <- split(input_tbl, factor(input_tbl$NCBI_ID))
    tree$Do(function(node) {
        if (!is.null(l[[node$name]])) {
            node[['table']] <- l[[node$name]]
        }
    })
    tree$Do(asr, traversal = 'post-order')
    tree$Do(inh, traversal = 'pre-order')

    ## Get NCBI_IDs with propagation results
    data_tree_tbl <- tree$Get(function(node) node[['table']], simplify = FALSE) |>
        purrr::discard(~ all(is.na(.x))) |>
        dplyr::bind_rows() |>
        dplyr::relocate(NCBI_ID) |>
        dplyr::filter(Evidence %in% c('', 'asr', 'inh') | is.na(Evidence))

    ## Combine data from propagation and original annotations
    ## >>Taxon name, etc are missing here<<
    data_with_values <- bind_rows(data_tree_tbl, x)

    ## Add missing values for the tree (this is maybe just necessary for
    ## displaying stats
    all_node_names <- tree$Get(function(node) node$name, simplify = TRUE)
    missing_node_names <- all_node_names[which(!all_node_names %in% unique(data_with_values$NCBI_ID))]
    if (length(missing_node_names > 0)) {
        attrs <- unique(x$Attribute)
        empty_df <- data.frame(
            NCBI_ID = sort(rep(missing_node_names, length(attrs))),
            Attribute = rep(attrs, length(missing_node_names)),
            Score = 0,
            Evidence = NA
        )
        final_table <- bind_rows(data_with_values, empty_df)
    } else {
        final_table <- data_with_values
    }

    attr_grp <- unique(x$Attribute_group)
    attr_type <- unique(x$Attribute_type)

    # check_id <- function(id) {
    #     tryCatch(
    #         error = function(e) NA,
    #         taxizedb::taxid2name(id, db = 'ncbi')
    #     )
    # }

    final_table <- final_table |>
        filter(NCBI_ID != 'ArcBac') |>
        mutate(
            Rank = case_when(
                grepl('^d__', NCBI_ID) ~ 'domain',
                grepl('^p__', NCBI_ID) ~ 'phylum',
                grepl('^c__', NCBI_ID) ~ 'class',
                grepl('^o__', NCBI_ID) ~ 'order',
                grepl('^f__', NCBI_ID) ~ 'family',
                grepl('^g__', NCBI_ID) ~ 'genus',
                grepl('^s__', NCBI_ID) ~ 'species',
                grepl('^t__', NCBI_ID) ~ 'strain'
            )
        ) |>
        mutate(NCBI_ID = sub('^[dpcofgst]__', '', NCBI_ID)) |>
        mutate(Taxon_name = ifelse(is.na(Taxon_name), checkTaxonName(NCBI_ID), Taxon_name)) |>
        filter(!is.na(NCBI_ID) & !is.na(Taxon_name)) |>
        mutate(Frequency = taxPPro:::scores2Freq(Score)) |>
        mutate(Attribute_value = TRUE) |>
        mutate(Attribute_group = attr_grp) |>
        mutate(Attribute_type = attr_type)

    tree$Do(function(node) {
        node[['table']] <- NULL
    })

    return(final_table)
})

full_dump_with_0 <- bind_rows(propagated)
full_dump_with_0$Attribute_value <- NULL
full_dump_with_0$Parent_name <- NULL
full_dump_with_0$Parent_rank <- NULL
full_dump_with_0$Parent_NCBI_ID <- NULL
full_dump_with_0$Strain <- NULL
full_dump_with_0$Genome_ID <- NULL
full_dump_with_0$Accession_ID <- NULL
full_dump_with_0 <- unique(full_dump_with_0)

## Some code for dividing the Attribute column into
## Attribute range
full_dump_with_0 <- full_dump_with_0 |>
    dplyr::mutate(
        Attribute_range = ifelse(
            test = grepl('\\(', Attribute),
            yes = sub('^.*(\\(.*)$', '\\1', Attribute),
            no = NA
        )
    ) |>
    dplyr::mutate(Attribute = sub('\\(.*$', '', Attribute)) |>
    dplyr::mutate(Attribute = stringr::str_squish(Attribute))

data.table::fwrite(
    x = full_dump_with_0, file = 'full_dump_with_0.csv', quote = TRUE, sep = ",",
    na = NA, row.names = FALSE, nThread = n_threads
)

pthreads <- paste0('-p', as.character(n_threads))
system2('pbzip2', args = list(pthreads, '-f', 'full_dump_with_0.csv'))

propagated <- map(propagated, ~ {
    total_scores <- .x |>
        group_by(Attribute_group, NCBI_ID) |>
        reframe(Total_score = sum(Score))
    taxids_above_0 <- total_scores |>
        filter(Total_score > 0) |>
        pull(NCBI_ID) |>
        unique()
    output <- .x |>
        filter(NCBI_ID %in% taxids_above_0)
    return(output)
})

full_dump <- bind_rows(propagated)
full_dump$Attribute_value <- NULL
full_dump_with_0$Parent_name <- NULL
full_dump_with_0$Parent_rank <- NULL
full_dump_with_0$Parent_NCBI_ID <- NULL
full_dump$Strain <- NULL
full_dump$Genome_ID <- NULL
full_dump$Accession_ID <- NULL
full_dump <- unique(full_dump)

## Same than above for fixint the Attribute column
full_dump <- full_dump |>
    dplyr::mutate(
        Attribute_range = ifelse(
            test = grepl('\\(', Attribute),
            yes = sub('^.*(\\(.*)$', '\\1', Attribute),
            no = NA
        )
    ) |>
    dplyr::mutate(Attribute = sub('\\(.*$', '', Attribute)) |>
    dplyr::mutate(Attribute = stringr::str_squish(Attribute))

data.table::fwrite(
    x = full_dump, file = 'full_dump.csv', quote = TRUE, sep = ",",
    na = NA, row.names = FALSE, nThread = n_threads
)
system2('pbzip2', args = list(pthreads, '-f', 'full_dump.csv'))

log_close()
