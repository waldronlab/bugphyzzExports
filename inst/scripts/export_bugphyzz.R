
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

phys_names <- c(
    ## multistate-intersection
    'aerophilicity',

    ## multistate-union
    # 'antimicrobial resistance',

    ## binary
    'acetate producing',

    ## numeric/range
    'growth temperature'
)

msg <- paste0(
    'Importing ', length(phys_names), ' physiologies for propagation: ',
    paste0(phys_names, collapse = ', '), '.'
)
log_print(msg, blank_after = TRUE)
bugphyzz_data <- physiologies(phys_names)
v <- sort(map_int(bugphyzz_data, nrow))
bugphyz_data <- bugphyzz_data[names(v)]

msg <- paste0(
    'Seaching for attributes of type range. They will be converted to type ',
    'multistate-intersection based on thresholds.'
)
log_print(msg, blank_after = TRUE)

phys <- vector('list', length(bugphyzz_data))
for (i in seq_along(phys)) {
    at <- unique(bugphyzz_data[[i]]$Attribute_type)
    dat_name <- names(bugphyzz_data)[i]
    names(phys)[i] <- dat_name
    if (at == 'range' && dat_name %in% names(THRESHOLDS())) {
        msg <- paste0(
            dat_name, " is of type range and we have a threshold for it.",
            ' Converting ', dat_name, ' to multistate-intersection.'
        )
        log_print(msg)
        res <- rangeToLogicalThr(bugphyzz_data[[i]], THRESHOLDS()[[dat_name]])
        res$Attribute_type <- 'multistate-intersection'
        phys[[i]] <- res

    } else if (at == 'range' && !dat_name %in% names(THRESHOLDS())) {
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
## TODO add message letting know when data is being eliminated
phys <- discard(phys, is.null)

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
})

## TODO add tidyr::complete for binary data
phys_data_ready <- map(phys_data_ready, ~ {
    attribute_type <- .x |>
        pull(Attribute_type) |>
        {\(y) y[!is.na(y)]}() |>
        unique()

    if (attribute_type %in% c('binary', 'multistate-union')) {
        return(completeBinaryData(.x))
    }
    return(.x)
})

log_print('', blank_after = TRUE)
log_print('Total time preparing data for propagation was: ')
log_print(tim, blank_after = TRUE)

taxidWarnings <- discard(taxidWarnings, is.null)
if (!is.null(taxidWarnings)) {
    msg <- 'Some NCBI IDs (taxids) need to be updated:'
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

    # tx <- grep('_taxid$', colnames(tip_data), value = TRUE)
    # nodes <- flatten(map(tx, ~ split(tip_data, factor(tip_data[[.x]]))))
    # nodes <- map(nodes, ~ .x[['tip_label']])
    # node_names <- map_int(nodes, ~ getMRCATaxPPro(tree, .x))
    # node_names <- node_names[!is.na(node_names)]
    # nodes_df <- data.frame(
    #     node = unname(node_names),
    #     node_label = names(node_names)
    # ) |>
    #     group_by(node) |>
    #     mutate(node_label = paste0(unique(node_label), collapse = '+')) |>
    #     ungroup() |>
    #     distinct()
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

    names(output)[i] <- attribute_group

    dat_n_tax <- length(unique(dat$NCBI_ID))
    msg <- paste0(
        attribute_group, ' has ', format(dat_n_tax, big.mark = ','), ' taxa.'
    )
    log_print(msg, blank_after = TRUE)

    msg <- paste0(
        'Mapping source annotations to the NCBI tree for ', attribute_group, '.'
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
        'Performing taxonomic pooling (round 1 of propagation) for ',
        attribute_type, '.'
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
        'Performing inhertiance1 (round 1 of propagation) for ',
        attribute_group, '.'
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
    if (per < 1) {
        msg <- paste0(
            'Not enough data for ASR. Skipping ASR for ', attribute_group,
            '. Stopped after the first round of propagation.'
        )
        log_print(msg, blank_after = TRUE)

        output[[i]] <- new_dat

        msg <- paste0('Cleaning nodes for ', attribute_type, '.')
        log_print(msg)
        tim <- system.time({
            ncbi_tree$Do(cleanNode)
        })
        log_print(tim, blank_after = TRUE)

        time2 <- Sys.time()
        time3 <- round(difftime(time2, time1, units = 'min'))
        nrow_fr <- nrow(new_dat)
        msg <- paste0(
            'Number of rows for ', attribute_type, ' were ' ,
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

    # pruned_tree <- ape::keep.tip(tree, tip = rownames(annotated_tips))
    # pruned_tree <- reorder(pruned_tree, 'postorder')
    # pruned_tip_data <- tip_data |>
    #     filter(tip_label %in% pruned_tree$tip.label)
    # pruned_node_data <- data.frame(
    #     node = length(pruned_tree$tip.label) + 1:pruned_tree$Nnode
    # )

    # tx <- grep('_taxid$', colnames(pruned_tip_data), value = TRUE)
    # nodes <- tx |>
    #     map(~ split(pruned_tip_data, factor(pruned_tip_data[[.x]]))) |>
    #     flatten() |>
    #     map(~ .x[['tip_label']])
    # node_names <- map_int(nodes, ~ getMRCATaxPPro(pruned_tree, .x))
    # node_names <- node_names[!is.na(node_names)]
    # nodes_df <- data.frame(
    #     node = unname(node_names),
    #     node_label = names(node_names)
    # ) |>
    #     group_by(node) |>
    #     mutate(node_label = paste0(unique(node_label), collapse = '+')) |>
    #     ungroup() |>
    #     distinct() |>
    #     arrange(node)
    # pruned_node_data <- left_join(pruned_node_data, nodes_df, by = 'node') |>
    #     mutate(
    #         node_label = ifelse(
    #             is.na(node_label), paste0('n', as.character(node)), node_label
    #         )
    #     )
    # pruned_tree$node.label <- pruned_node_data$node_label

    msg <- paste0(
        'Performing ASR for (round 2 of propagation) ', attribute_type, '.'
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
    # nodes_annotated <- res[which(grepl('^\\d+(\\+\\d+)*', rownames(res))),]
    node_data_annotated <- ltp$node_data |>
        filter(node_label %in% unique(res_df$node_label)) |>
        select(node_label, taxid, Taxon_name, Rank)

    nodes_annotated <- node_data_annotated |>
        left_join(res_df, by = 'node_label') |>
        # select(taxid, Taxon_name, Rank) |>
        mutate(Rank = ifelse(Rank == 'superkingdom', 'kingdom', Rank)) |>
        # new_taxa_from_nodes <- nodes_annotated |>
        # mutate(Rank = taxizedb::taxid2rank(taxid)) |>
        # mutate(Rank = ifelse(Rank == 'superkingdom', 'kingdom', Rank)) |>
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
            # taxid = sub('\\w__', '', NCBI_ID),
            # Taxon_name = taxizedb::taxid2name(taxid, db = 'ncbi'),
            Frequency = case_when(
                Score == 1 ~ 'always',
                Score > 0.9 ~ 'usually',
                Score >= 0.5 ~ 'sometimes',
                Score > 0 & Score < 0.5 ~ 'rarely',
                Score == 0 ~ 'never'
            )
        )



    # new_taxa_from_nodes <- nodes_annotated |>
    #     as.data.frame() |>
    #     tibble::rownames_to_column(var = 'NCBI_ID') |>
    #     filter(grepl('^\\d+(\\+\\d+)*', NCBI_ID)) |>
    #     mutate(NCBI_ID = strsplit(NCBI_ID, '\\+')) |>
    #     tidyr::unnest(NCBI_ID) |>
    #     mutate(Rank = taxizedb::taxid2rank(NCBI_ID)) |>
    #     mutate(Rank = ifelse(Rank == 'superkingdom', 'kingdom', Rank)) |>
    #     mutate(
    #         NCBI_ID = case_when(
    #             Rank == 'kingdom' ~ paste0('k__', NCBI_ID),
    #             Rank == 'phylum' ~ paste0('p__', NCBI_ID),
    #             Rank == 'class' ~ paste0('c__', NCBI_ID),
    #             Rank == 'order' ~ paste0('o__', NCBI_ID),
    #             Rank == 'family' ~ paste0('f__', NCBI_ID),
    #             Rank == 'genus' ~ paste0('g__', NCBI_ID),
    #             Rank == 'species' ~ paste0('s__', NCBI_ID),
    #             Rank == 'strain' ~ paste0('t__', NCBI_ID)
    #         )
    #     ) |>
    #     filter(
    #         Rank %in% c(
    #             'kingdom', 'phylum', 'class', 'order', 'family', 'genus',
    #             'species', 'strain'
    #         )
    #     ) |>
    #     mutate(Evidence = 'asr') |>
    #     relocate(NCBI_ID, Rank, Evidence) |>
    #     pivot_longer(
    #         cols = 4:last_col(), names_to = 'Attribute', values_to = 'Score'
    #     ) |>
    #     mutate(
    #         Attribute_source = NA,
    #         Confidence_in_curation = NA,
    #         Attribute_group = Attribute_group_var,
    #         Attribute_type = Attribute_type_var,
    #         # taxid = sub('\\w__', '', NCBI_ID),
    #         Taxon_name = taxizedb::taxid2name(taxid, db = 'ncbi'),
    #         Frequency = case_when(
    #             Score == 1 ~ 'always',
    #             Score > 0.9 ~ 'usually',
    #             Score >= 0.5 ~ 'sometimes',
    #             Score > 0 & Score < 0.5 ~ 'rarely',
    #             Score == 0 ~ 'never'
    #         )
    #     )
    new_taxa_for_ncbi_tree <- nodes_annotated |>
        relocate(NCBI_ID, Rank, Attribute, Score, Evidence)
    new_taxa_for_ncbi_tree_list <- split(
        new_taxa_for_ncbi_tree, factor(new_taxa_for_ncbi_tree$NCBI_ID)
    )

    ## Perform taxonomic pooling and inheritance (propagation round 3)
    ## Mapping new internal nodes to the NCBI taxonomy tree ####
    msg <- paste0(
        'Mapping annotations for third round of propagation for ', attribute_type,
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
        'Performing inheritance (round 3 of propagation) for ', attribute_type
    )
    log_print(msg)
    tim <- system.time({
        ncbi_tree$Do(inh2, traversal = 'pre-order')
    })
    log_print(tim, blank_after = TRUE)

    ## Extracting files ####
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
        'Minimum threshold for positives in ', attribute_group, ' was ',
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
        'Size of propagated data for ', attribute_type, ' is ',
        gdata::humanReadable(final_result_size, standard = 'SI'), '.'
    )
    log_print(msg, blank_after = TRUE)

    output[[i]] <- final_result


    msg <- paste0('Cleaning nodes for ', attribute_type, '.')
    log_print(msg)
    tim <- system.time({
        ncbi_tree$Do(cleanNode)
    })
    log_print(tim, blank_after = TRUE)

    time2 <- Sys.time()
    time3 <- round(difftime(time2, time1, units = 'min'))
    nrow_fr <- nrow(final_result)
    msg <- paste0(
        'Number of rows for ', attribute_type, ' were ' ,
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
    'The size of the tsv file is ', fsize, '.'
)
log_print(msg, blank_after = TRUE)

si <- sessioninfo::session_info()
log_print(si, blank_after = TRUE)
log_close()

