
## Script for importing attribute data from bugphyzz,
## performing cleaning and propagation,
## and exporting the data as a single tsv text file and
## a set of .gmt text files.

## Setup #######################################################################
suppressMessages({
    library(logr)
    library(bugphyzz)
    library(taxPPro)
    library(data.tree)
    library(phytools)
    library(dplyr)
    library(purrr)
    library(tidyr)
    library(ggplot2)
    library(ape)
    library(bugphyzzExports)
})
logfile <- "log_file"
lf <- log_open(logfile, logdir = FALSE, compact = TRUE, show_notes = FALSE)

## Import data #################################################################
phys_names <- c(

    ## multistate-intersection
    'aerophilicity',
    # 'gram stain',
    # 'biosafety level',
    # 'COGEM pathogenicity rating',
    # 'shape',
    # 'spore shape',
    # 'arrangement',
    # 'hemolysis', # didn't run for this (It must be run independently for cross-validation)

    ## multistate-union
    'habitat',
    # 'disease association',
    # 'antimicrobial resistance',
    # 'halophily', ## Curation must be reviewed

    ## multistate-uninion (but not propagation for these)
    # 'isolation site', Do not include. Curation must be reviewed. Compare with habitat.
    # 'growth medium', Do not include. Curation must be reviewed.
    # 'country', Do not include.
    # 'geographic location', Do not include.

    # 'metabolite production', Curation must be reviewed before inclusion.
    # 'metabolite utilization' Curation must be reviewed before inclusion.

    ## binary
    'plant pathogenicity',
    # 'acetate producing',
    # 'sphingolipid producing',
    # 'lactate producing',
    # 'butyrate producing',
    # 'hydrogen gas producing',
    # 'pathogenicity human',
    # 'motility',
    # 'biofilm forming',
    # 'extreme environment',
    # 'animal pathogen',
    # 'antimicrobial sensitivity',
    # 'spore formation', # didn' run for this (run independently for cross-validation)
    # 'health associated', # didn't run for this (run independently for cross-validation)

    ## numeric/range
    'growth temperature'
    # 'optimal ph',
    # 'width',
    # 'length',
    # 'genome size',
    # 'coding genes',
    # 'mutation rate per site per generation',
    # 'mutation rate per site per year'
)
msg <- paste0(
    'Importing ', length(phys_names), ' physiologies for propagation: ',
    paste0(phys_names, collapse = ', '), '.'
)
log_print(msg, blank_after = TRUE)
bugphyzz_data <- physiologies(phys_names)
v_order <- sort(map_int(bugphyzz_data, nrow))
bugphyzz_data <- bugphyzz_data[names(v_order)]

## Convert range/numeric physiologies to categorical ###########################
msg <- paste0(
    'Searching for attributes of type range. They will be converted to type ',
    'multistate-intersection based on thresholds.'
)
log_print(msg, blank_after = TRUE)
phys <- vector('list', length(bugphyzz_data))
for (i in seq_along(phys)) {
    attribute_type <- bugphyzz_data[[i]]$Attribute_type |>
        {\(y) y[!is.na(y)]}() |>
        unique()
    dat_name <- names(bugphyzz_data)[i]
    names(phys)[i] <- dat_name
    if (attribute_type == 'range' && dat_name %in% names(THRESHOLDS())) {
        msg <- paste0(
            dat_name, " is of type range and we have a threshold for it.",
            ' Converting ', dat_name, ' to multistate-intersection.'
        )
        log_print(msg)
        res <- rangeToLogicalThr(bugphyzz_data[[i]], THRESHOLDS()[[dat_name]])
        res$Attribute_type <- 'multistate-intersection'
        phys[[i]] <- res

    } else if (attribute_type == 'range' && !dat_name %in% names(THRESHOLDS())) {
        msg <- paste0(
            dat_name, " is of type range, but we don't have a threshold for it.",
            " Skipping ", dat_name, '.'
        )
        log_print(msg)
        next

    } else {
        phys[[i]] <- bugphyzz_data[[i]]

    }
}
for (i in seq_along(phys)) {
   physName <-  names(phys)[i]
   if (is.null(phys[[i]])) {
       msg <- paste0(physName, ' will be discarded now. Not in thresholds.')
       log_print(msg)
   }
}
log_print("", blank_after = TRUE)
phys <- discard(phys, is.null)

## Check valid attributes ######################################################
msg <- paste(
    'Check that all attributes are valid.',
    'Invalid values will not be exported.',
    'Invalid values can be checked on the log file.'
)
log_print(msg, blank_after = TRUE)
filterAttributes <- function(dat) {
    ## This function is for filtering valid attributes based on the
    ## attributes.tsv file in the bugphyzz R package.
    fpath <- file.path('extdata', 'attributes.tsv')
    attributes_tsv <- system.file(fpath, package = 'bugphyzz')
    attributes_data <- unique(read.table(attributes_tsv, header = TRUE, sep = '\t'))
    ag <- dat |>
        pull(Attribute_group) |>
        {\(y) y[!is.na(y)]}() |>
        unique()
    current_attrs <- dat |>
        pull(Attribute) |>
        {\(y) y[!is.na(y)]}() |>
        unique()
    rgx <- paste0('\\b', ag, '\\b')
    valid_attributes <- attributes_data |>
        filter(grepl(rgx, attribute_group)) |>
        pull(attribute)

    if (!length(valid_attributes)){
        msg <- paste0(
            ag, ' not in valid attribute groups. It needs to be added.'
        )
        log(msg)
        return(NULL)
    }

    lgl_v <- !current_attrs %in% valid_attributes
    if (any(lgl_v)) {
        invalid_values <- current_attrs[which(lgl_v)]
        msg <- paste0(
            'These values are not valid for ', ag, ': ',
            paste0(invalid_values, collapse = ', '), '.'
        )
        log_print(msg)
    }
    output <- dat |>
        filter(Attribute %in% valid_attributes)
    if (!length(output)) {
        msg <- paste0('No valid attributes for ', ag, '. Dropping it.')
        log_print(msg)
        output <- NULL
    }
    return(output)
}
phys <- map(phys, filterAttributes)
for (i in seq_along(phys)) {
    if (!nrow(phys[[i]])) {
        msg <- paste0(
            names(phys)[i], ' discarde due to the lack of valid attributes.',
            ' Please check the log file.'
        )
    }
}
log_print("", blank_after = TRUE)
phys <- discard(phys, is.null)

## Preparing data for propagation ##############################################
msg <- ('Preparing data for propagation...')
log_print('', blank_after = TRUE)
log_print(msg, blank_after = TRUE)
tim <- system.time({
    phys_data_ready <- vector('list', length(phys))
    taxidWarnings <- vector('list', length(phys))
    for (i in seq_along(phys_data_ready)) {
        name <- names(phys)[i]
        msg <- paste0('Preparing ', name, '.')
        log_print(msg)
        names(phys_data_ready)[i] <- name
        names(taxidWarnings)[i] <- name
        wngs <- list()
        suppressWarnings({
            withCallingHandlers({
                dat <- getDataReady(filterData(phys[[i]]))
                if (length(dat) > 0)
                    phys_data_ready[[i]] <- dat
            },
            warning = function(w) {
                if (grepl('taxizedb', w$message)) {
                    msg <- sub('.*unrank.*: (\\d+.*)$', '\\1', w$message)
                    wngs <<- c(wngs, list(msg))
                }
            })
        })
        if (length(wngs) > 0)
            taxidWarnings[[i]] <- wngs
    }
    phys_data_ready <- list_flatten(phys_data_ready)
    ## list_flatten (line above) ensures that data from attributes of type
    ## multistate-union are separated into individual data.frames per attribute
})


for (i in seq_along(phys_data_ready)) {
   physName <-  names(phys_data_ready)[i]
   if (is.null(phys_data_ready[[i]])) {
       msg <- paste0(phyName, ' will be discarded now. No taxids.')
       log_print(msg)
   }
}
log_print("", blank_after = TRUE)
phys_data_ready <- discard(phys_data_ready, is.null)

phys_data_ready <- map(phys_data_ready, ~ {
    attribute_type <- .x |>
        pull(Attribute_type) |>
        {\(y) y[!is.na(y)]}() |>
        unique()
    if (attribute_type %in% c('binary', 'multistate-union')) {
        ## This step is used to include FALSE values,
        ## which are necessary for the ASR step below
        ## In the case of multistate-intersection, FALSE values are inferred because they're mutally exclusive (need to elaborate more here).
        return(completeBinaryData(.x))
    }
    return(.x)
})

log_print('', blank_after = TRUE)
log_print('Total time preparing data for propagation was: ')
log_print(tim, blank_after = TRUE)

taxidWarnings <- discard(taxidWarnings, is.null)
if (!is.null(taxidWarnings)) {
    msg <- 'Some NCBI IDs (taxids) might need to be updated:'
    log_print(msg, blank_after = TRUE)
    log_print(taxidWarnings, blank_after = TRUE)
}

msg <- paste0('Preparing tree data (NCBI and LTP).')
log_print(msg)
tim <- system.time({
    data('tree_list')
    ncbi_tree <- as.Node(tree_list)
    ltp <- ltp()
    tree <- reorder(ltp$tree, 'postorder')
    tip_data <- ltp$tip_data
    node_data <- ltp$node_data

})
log_print(tim, blank_after = TRUE)

start_time <- Sys.time()
msg <- paste0('Performing propagation. It started at ', start_time, '.')
log_print(msg, blank_after = TRUE)

output <- vector('list', length(phys_data_ready))
for (i in seq_along(phys_data_ready)) {
    time1 <- Sys.time()
    dat <- phys_data_ready[[i]]

    attribute_group <- dat$Attribute_group |>
        {\(y) y[!is.na(y)]}() |>
        unique()
    attribute_type <- dat$Attribute_type |>
        {\(y) y[!is.na(y)]}() |>
        unique()
    attribute_nms <- dat$Attribute |>
        {\(y) y[!is.na(y)]}() |>
        unique()

    names(output)[i] <- names(phys_data_ready)[i]

    if (attribute_type == 'multistate-union') {
        attrNMS <- unique(sub('--(TRUE|FALSE)$', '', attribute_nms))
        attrGroupMsg <- paste0(attribute_group, '-', attrNMS)
    } else {
        attrGroupMsg <- attribute_group
    }

    dat_n_tax <- length(unique(dat$NCBI_ID))
    msg <- paste0(
        attrGroupMsg, ' has ', format(dat_n_tax, big.mark = ','), ' taxa.'
    )
    log_print(msg, blank_after = TRUE)

    msg <- paste0(
        'Mapping source annotations to the NCBI tree for ', attrGroupMsg, '.'
    )
    log_print(msg)

    node_list <- split(dat, factor(dat$NCBI_ID))
    tim <- system.time({
        ncbi_tree$Do(function(node) {
            if (node$name %in% names(node_list))
                node$attribute_tbl <- node_list[[node$name]]
        })
    })
    log_print(tim, blank_after = TRUE)

    msg <- paste0(
        'Performing taxonomic pooling for ',
        attrGroupMsg, '.'
    )
    log_print(msg)
    tim <- system.time({
        ncbi_tree$Do(
           function(node) {
                taxPool(
                    node = node,
                    grp = attribute_group,
                    typ = attribute_type)
            },
            traversal = 'post-order'
        )
    })
    log_print(tim, blank_after = TRUE)

    msg <- paste0(
        'Performing inheritance (1) for ',
        attrGroupMsg, '.'
    )
    log_print(msg)

    tim <- system.time({
        ncbi_tree$Do(inh1, traversal = 'pre-order')
    })
    log_print(tim, blank_after = TRUE)

    new_dat <- ncbi_tree$Get(
        'attribute_tbl', filterFun = function(node) {
            grepl('^[gst]__', node$name)
        }
    ) |>
        discard(~ all(is.na(.x))) |>
        bind_rows() |>
        arrange(NCBI_ID, Attribute) |>
        filter(!NCBI_ID %in% dat$NCBI_ID) |>
        bind_rows(dat)

    new_taxids <- new_dat |>
        pull(taxid) |>
        {\(y) y[!is.na(y)]}()

    per <- mean(tip_data$taxid %in% new_taxids) * 100
    if (per < 10) {
        msg <- paste0(
            'Not enough data for ASR. Skipping ASR and inhetiance (2) for ', attrGroupMsg,
            '. Stopped after the first round of propagation.'
        )
        log_print(msg, blank_after = TRUE)

        output[[i]] <- new_dat ## Here, I include data after taxpool and inhertiance 1, ASR and inheritance 2 are skipped

        msg <- paste0('Cleaning nodes for ', attrGroupMsg, '.')
        log_print(msg)
        tim <- system.time({
            ncbi_tree$Do(cleanNode)
        })
        log_print(tim, blank_after = TRUE)

        time2 <- Sys.time()
        time3 <- round(difftime(time2, time1, units = 'min'))
        nrow_fr <- nrow(new_dat)
        msg <- paste0(
            'Number of rows for ', attrGroupMsg, ' were ' ,
            format(nrow_fr, big.mark = ','), '.',
            ' It took ', time3[[1]], ' mins.'
        )
        log_print(msg, blank_after = TRUE)
        log_print('', blank_after = TRUE)

        next
    }

    tip_data_annotated <- left_join(
        tip_data,
        select(new_dat, taxid, Attribute, Score),
        by = 'taxid'
    )

    annotated_tips <- tip_data_annotated |>
        select(tip_label, Attribute, Score) |>
        filter(!is.na(Attribute)) |>
        pivot_wider(
            names_from = 'Attribute', values_from = 'Score', values_fill = 0
        ) |>
        tibble::column_to_rownames(var = 'tip_label') |>
        as.matrix()

    if (attribute_type %in% c('binary', 'multistate-union')) {
        no_annotated_tips <- tip_data |>
            filter(!tip_label %in% rownames(annotated_tips)) |>
            select(tip_label) |>
            mutate(
                Attribute = factor(
                    attribute_nms[[1]], levels = attribute_nms
                )
            ) |>
            complete(tip_label, Attribute) |>
            mutate(Score = ifelse(grepl('--FALSE$', Attribute), 1, 0)) |>
            pivot_wider(names_from = 'Attribute', values_from = 'Score') |>
            tibble::column_to_rownames(var = 'tip_label') |>
            as.matrix() |>
            {\(y) y[,colnames(annotated_tips)]}()
    } else if (attribute_type == 'multistate-intersection') {
        no_annotated_tip_names <- tip_data |>
            filter(!tip_label %in% rownames(annotated_tips)) |>
            pull(tip_label)
        fill_value <- 1 / length(attribute_nms)
        vct <- rep(
            fill_value,
            length(no_annotated_tip_names) * length(attribute_nms)
        )
        no_annotated_tips <- matrix(
            data = vct,
            nrow = length(no_annotated_tip_names),
            ncol = length(attribute_nms),
            dimnames = list(no_annotated_tip_names, attribute_nms)
        )
    }

    input_matrix <- rbind(annotated_tips, no_annotated_tips)
    input_matrix <- input_matrix[tree$tip.label,]

    msg <- paste0(
        'Performing ASR for (round 2 of propagation) ', attrGroupMsg, '.'
    )
    log_print(msg)
    tim <- system.time({
        fit <- fitMk(
            tree = tree, x = input_matrix, model = 'ER',
            pi = 'fitzjohn', lik.func = 'pruning', logscale = TRUE
        )
        asr <- ancr(object = fit, tips = TRUE)
    })
    log_print(tim, blank_after = TRUE)

    res <- asr$ace
    node_rows <- length(tree$tip.label) + 1:tree$Nnode
    rownames(res)[node_rows] <- tree$node.label
    res <- res[tree$node.label,]
    res_df <- res |>
        as.data.frame() |>
        tibble::rownames_to_column(var = 'node_label') |>
        filter(!grepl('^n\\d+', node_label))

    ## Get annotations for nodes
    node_data_annotated <- ltp$node_data |>
        filter(node_label %in% unique(res_df$node_label)) |>
        select(node_label, taxid, Taxon_name, Rank)

    nodes_annotated <- node_data_annotated |>
        left_join(res_df, by = 'node_label') |>
        mutate(Rank = ifelse(Rank == 'superkingdom', 'kingdom', Rank)) |>
        mutate(
            NCBI_ID = case_when(
                Rank == 'kingdom' ~ paste0('k__', taxid),
                Rank == 'phylum' ~ paste0('p__', taxid),
                Rank == 'class' ~ paste0('c__', taxid),
                Rank == 'order' ~ paste0('o__', taxid),
                Rank == 'family' ~ paste0('f__', taxid),
                Rank == 'genus' ~ paste0('g__', taxid),
                Rank == 'species' ~ paste0('s__', taxid),
                Rank == 'strain' ~ paste0('t__', taxid)
            )
        ) |>
        filter(
            Rank %in% c(
                'kingdom', 'phylum', 'class', 'order', 'family', 'genus',
                'species', 'strain'
            )
        ) |>
        mutate(Evidence = 'asr') |>
        relocate(NCBI_ID, taxid, Taxon_name, Rank, Evidence) |>
        pivot_longer(
            cols = 7:last_col(), names_to = 'Attribute', values_to = 'Score'
        ) |>
        mutate(
            Attribute_source = NA,
            Confidence_in_curation = NA,
            Attribute_group = attribute_group,
            Attribute_type = attribute_type,
            Frequency = case_when(
                Score == 1 ~ 'always',
                Score > 0.9 ~ 'usually',
                Score >= 0.5 ~ 'sometimes',
                Score > 0 & Score < 0.5 ~ 'rarely',
                Score == 0 ~ 'never'
            )
        ) |>
        select(-node_label)

    new_taxa_for_ncbi_tree <- nodes_annotated |>
        relocate(NCBI_ID, Rank, Attribute, Score, Evidence)

    new_taxa_for_ncbi_tree_list <- split(
        new_taxa_for_ncbi_tree, factor(new_taxa_for_ncbi_tree$NCBI_ID)
    )

    msg <- paste0(
        'Mapping annotations for inheritabce 2 for ', attrGroupMsg,
        '.'
    )
    log_print(msg)
    tim <- system.time({
        ncbi_tree$Do(function(node) {
            cond1 <- node$name %in% names(new_taxa_for_ncbi_tree_list)
            cond2 <- is.null(node$attribute_tbl) || all(is.na(node$attribute_tbl))
            if (cond1 && cond2) {
                node$attribute_tbl <- new_taxa_for_ncbi_tree_list[[node$name]]
            }
        })
    })
    log_print(tim, blank_after = TRUE)

    msg <- paste0(
        'Performing inheritance (2)  for ', attrGroupMsg
    )
    log_print(msg)
    tim <- system.time({
        ncbi_tree$Do(inh2, traversal = 'pre-order')
    })
    log_print(tim, blank_after = TRUE)

    result <- ncbi_tree$Get(
        attribute = 'attribute_tbl', simplify = FALSE,
        filterFun = function(node) {
            node$name != 'ArcBac' && !is.null(node$attribute_tbl)
        }
    ) |>
        bind_rows() |>
        discard(~ all(is.na(.x)))
    min_thr <- 1 / length(unique(dat$Attribute))

    msg <- paste0(
        'Minimum threshold for positives in ', attrGroupMsg, ' was ',
        min_thr, '.'
    )
    log_print(msg, blank_after = TRUE)

    add_taxa_1 <- dat |>
        filter(!NCBI_ID %in% unique(result$NCBI_ID)) |>
        discard(~ all(is.na(.x)))
    add_taxa_2 <- new_taxa_for_ncbi_tree |>
        filter(!NCBI_ID %in% unique(result$NCBI_ID)) |>
        discard(~ all(is.na(.x)))
    final_result <- bind_rows(list(result, add_taxa_1, add_taxa_2)) |>
        filter(Score > min_thr)

    final_result_size <- lobstr::obj_size(final_result)
    msg <- paste0(
        'Size of propagated data for ', attrGroupMsg, ' is ',
        gdata::humanReadable(final_result_size, standard = 'SI'), '.'
    )
    log_print(msg, blank_after = TRUE)

    output[[i]] <- final_result


    msg <- paste0('Cleaning nodes for ', attrGroupMsg, '.')
    log_print(msg)
    tim <- system.time({
        ncbi_tree$Do(cleanNode)
    })
    log_print(tim, blank_after = TRUE)

    time2 <- Sys.time()
    time3 <- round(difftime(time2, time1, units = 'min'))
    nrow_fr <- nrow(final_result)
    msg <- paste0(
        'Number of rows for ', attrGroupMsg, ' were ' ,
        format(nrow_fr, big.mark = ','), '.',
        ' It took ', time3[[1]], ' mins.'
    )
    log_print(msg, blank_after = TRUE)
    log_print('', blank_after = TRUE)
}
end_time <- Sys.time()
elapsed_time <- round(difftime(end_time, start_time, units = 'min'))

msg <- paste0(
    'Propagation ended at ', elapsed_time,
    '. Total elapsed time for propagtion for ', length(phys_data_ready),
    ' physiologies was ', elapsed_time[[1]], ' min.'
)
log_print(msg, blank_after = TRUE)

## Exporting annotations as a single tsv file ####
final_obj <- bind_rows(output)
final_obj_size <- lobstr::obj_size(final_obj)

msg <- paste0(
    'Size of final object is ',
    gdata::humanReadable(final_obj_size, standard = 'SI')
)
log_print(msg, blank_after = TRUE)

msg <- paste0('Writing final output file.')
log_print(msg, blank_after = TRUE)
final_obj_fname <- paste0('bugphyzz_export_', Sys.Date(), '.tsv')
write.table(
    x = final_obj, file = final_obj_fname, sep = '\t', row.names = FALSE
)

fsize <- gdata::humanReadable(file.size(final_obj_fname), standard = "SI")
msg <- paste0(
    'The size of the tsv file is ', fsize, '. Output file name is ',
    final_obj_fname, '.'
)
log_print(msg, blank_after = TRUE)

si <- sessioninfo::session_info()
log_print(si, blank_after = TRUE)
log_close()
