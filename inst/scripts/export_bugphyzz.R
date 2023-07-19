## Script to create bugphyzz exports dump files and signature files

library(logr)
library(bugphyzz) # BiocManager::install('waldronlab/bugphyzz', force = TRUE)
library(taxPPro)
library(purrr)
library(rlang)
library(dplyr)
library(data.tree)
library(bugphyzzExports) # BiocManager::install('sdgamboa/taxPPro', force = TRUE)
library(BiocParallel)
library(tidyr)

logfile <- "log_file"
lf <- log_open(logfile, logdir = FALSE, compact = TRUE, show_notes = FALSE)

exlclude_phys <- c(
    'country', 'geographic location',
    'habitat', 'isolation site',
    'metabolite production', 'metabolite utilization',
    'halophily' # I think this should be separated in optimal and not optimal, or just optimal so it's similar to optimal ph
)

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
    if (attr_grp %in% binaries) { # binaries is defined above
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

categorical <- keep(phys, ~ unique(.x$Attribute_type) == 'logical')
categorical$aerophilicity <- homogenizeAerophilicityAttributeNames(
    categorical$aerophilicity
)
range <- keep(phys, ~ unique(.x$Attribute_type == 'range'))
range <- range[which(names(range) %in% names(THRESHOLDS()))]
range_cat <- map2(range, names(range), ~ rangeToLogicalThr(.x, THRESHOLDS()[[.y]]))
categorical <- c(categorical, range_cat)

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

data <- bplapply(data, BPPARAM = MulticoreParam(workers = 60), FUN = function(x) {
    set1 <- filter(x, !is.na(NCBI_ID))
    set2 <- filter(x, is.na(NCBI_ID))
    set1$Rank <- checkRank(set1$NCBI_ID)
    set1$Taxon_name <- checkTaxonName(set1$NCBI_ID)
    bind_rows(set1, set2)
})

data <- map(data, ~ {
    if ('Attribute_range' %in% colnames(.x)) {
        .x$Attribute <- paste0(.x$Attribute, ' ', .x$Attribute_range)
        .x$Attribute <- sub('\\)$', '', .x$Attribute)
        .x$Attribute <- paste0(.x$Attribute, ' ', .x$Unit, ')')
        .x$Attribute <- sub(' \\)', ')', .x$Attribute)
    }
    return(.x)
})

data_ready <- bplapply(data, BPPARAM = MulticoreParam(workers = 60), FUN = function(x) {
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
            # if ('Unit' %in% colnames(x)) {
            #     unit <- x$Unit
            #     unit <- unique(unit[!is.na(unit)])
            #     if (length(unit) > 1) {
            #         attr_grp <- unique(x$Attribute_group)
            #         warning('More than 1 unit for', attr_grp, call. = FALSE)
            #     }
            #     output$Unit <- unit
            # }
            return(output)
        }
    )
})
data_ready <- discard(data_ready, is_error)

data('tree_list')
tree <- as.Node(tree_list)

propagated <- bplapply(X = data_ready, BPPARAM = MulticoreParam(workers = 60), FUN = function(x) {

    # if ('Unit' %in% colnames(x)) {
    #     unit <- x$Unit
    #     unit <- unique(unit[!is.na(unit)])
    #     if (length(unit) > 1) {
    #         attr_grp <- unique(x$Attribute_group)
    #         warning('More than 1 unit for', attr_grp, call. = FALSE)
    #     }
    # } else {
    #     unit <- NULL
    # }

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

    check_id <- function(id) {
        tryCatch(
            error = function(e) NA,
            taxizedb::taxid2name(id, db = 'ncbi')
        )
    }

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
        mutate(Taxon_name = ifelse(is.na(Taxon_name), check_id(NCBI_ID), Taxon_name)) |>
        filter(!is.na(NCBI_ID) & !is.na(Taxon_name)) |>
        mutate(Frequency = taxPPro:::scores2Freq(Score)) |>
        mutate(Attribute_value = TRUE) |>
        mutate(Attribute_group = attr_grp) |>
        mutate(Attribute_type = attr_type)

    # if (!is.null(unit)) {
    #     final_table$Unit <- unit
    # }

    tree$Do(function(node) {
        node[['table']] <- NULL
    })

    return(final_table)
})

full_dump_with_0 <- bind_rows(propagated)
full_dump_with_0$Attribute_value <- NULL
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
    na = NA, row.names = FALSE, nThread = 60
)
system2('pbzip2', args = list('-p60', '-f', 'full_dump_with_0.csv'))

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
# readr::write_csv(
#     x = full_dump, file = "full_dump.csv.bz2", quote = 'needed', num_threads = 16
# )

# map(fu)
#
# total_scores <- full_dump |>
#     group_by(Attribute_group, NCBI_ID) |>
#     reframe(Total_score = sum(Score))
#
# taxids_above_0 <- total_scores |>
#     filter(Total_score > 0)
#
# full_dump <- full_dump |>
#     filter(NCBI_ID %in% taxids_above_0)

data.table::fwrite(
    x = full_dump, file = 'full_dump.csv', quote = TRUE, sep = ",",
    na = NA, row.names = FALSE, nThread = 60
)
system2('pbzip2', args = list('-p60', '-f', 'full_dump.csv'))

log_close()
