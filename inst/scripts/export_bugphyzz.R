## Script to create bugphyzz exports dump files and signature files in gmt format

suppressMessages({
    library(logr)
    library(bugphyzz)
    library(taxPPro)
    library(rlang)
    library(dplyr)
    library(bugphyzzExports)
    library(tidyr)
    library(tibble)
    library(phytools)
    library(castor)
    library(purrr)
})

logfile <- "log_file"
lf <- log_open(logfile, logdir = FALSE, compact = TRUE, show_notes = FALSE)

attributes_by_type <- list(
    numeric = c( # Attribute_type is range
        "growth temperature",
        "coding genes",
        "genome size",
        "length",
        "width",
        "optimal ph",
        "mutation rate per site per generation",
        "mutation rate per site per year"
    ),
    binary = c( # Attribute_type is binary
        "animal pathogen",
        "antimicrobial sensitivity",
        "biofilm forming",
        "extreme environment",
        "health associated",
        "hydrogen gas producing",
        "lactate producing",
        "motility",
        "plant pathogenicity",
        "spore formation",
        "host-associated"
    ),
    multistate = c( # Attribute_type is multistate-intersection
        "aerophilicity",
        "gram stain",
        "biosafety level",
        "COGEM pathogenicity rating",
        "shape",
        "spore shape",
        "arrangement",
        "hemolysis"
    )
    # do_not_propagate = c(
    #     "habitat"
    # )
)


ltp <- ltp()
tree <- ltp$tree
tip_data <- ltp$tip_data
node_data <- ltp$node_data

tree_data <- bind_rows(
    select(as_tibble(tip_data), label = tip_label, NCBI_ID, Taxon_name, Rank, taxid),
    select(as_tibble(node_data), label = node_label, NCBI_ID, Taxon_name, Rank, taxid)
)

phys_names <- unlist(attributes_by_type, use.names = FALSE)
msg <- paste0('"', paste0(phys_names, collapse = ', '), '"')
msg_len <- length(phys_names)
msg <- paste('Importing', msg_len, 'physiologies from bugphyzz:', msg, '--', Sys.time())
log_print(msg, blank_after = TRUE)
system.time(phys <- physiologies(phys_names, full_source = FALSE))

dat_ready <- vector("list", length(phys))
for (i in seq_along(dat_ready)) {
    message(names(phys)[i])
    dat_ready[[i]] <- getDataReady(filterData(phys[[i]]))
    names(dat_ready)[i] <- names(phys)[i]
}

# backup <- dat_ready
# dat_ready <- dat_ready[c("aerophilicity", "growth temperature", "animal pathogen")]
# dat_ready <- backup

propagated <- vector('list', length(dat_ready))
for (i in seq_along(propagated)) {
    dat <- dat_ready[[i]]
    attr_grp <- dat |>
        pull(Attribute_group) |>
        {\(y) y[!is.na(y)]}() |>
        unique()
    attr_type <- dat |>
        pull(Attribute_type) |>
        {\(y) y[!is.na(y)]}() |>
        unique()
    msg <- paste0("Propagating ", attr_grp, " of type ", attr_type, ".")
    log_print(msg)

    fdat <- left_join(
        tip_data, dat_ready[[i]], by = 'NCBI_ID',
        relationship = "many-to-many"
    )

    if (attr_type %in% c("multistate-intersection", "binary")) {
        ## Propagation multistate ####
        annotated_tips <- fdat |>
            filter(!is.na(Attribute)) |> # NAs in the Attribute column correspond to unnanotated tips
            select(tip_label, Attribute, Score) |>
            pivot_wider(
                names_from = "Attribute", values_from = "Score", values_fill = 0
            ) |>
            column_to_rownames(var = "tip_label") |>
            as.data.frame() |>
            as.matrix()
        no_annotated_tips_names <- fdat |>
            filter(is.na(Attribute)) |>
            pull(tip_label)
        no_annotated_tips <- matrix(
            data = rep(1/ncol(annotated_tips), length(no_annotated_tips_names) * ncol(annotated_tips)),
            ncol = ncol(annotated_tips),
            dimnames = list(rownames = no_annotated_tips_names, colnames = colnames(annotated_tips))
        )
        input_mat <- rbind(annotated_tips, no_annotated_tips)
        input_mat <- input_mat[tree$tip.label, ]

        tim <- system.time({
            fit <- fitMk(tree = tree, x = input_mat, model = "ER", pi = "equal")
            ace <- ancr(fit, tips = TRUE)
        })
        log_print(tim)
        asr_df <- ace$ace |>
            as.data.frame() |>
            rownames_to_column(var = "label") |>
            filter(label %in% tree_data$label) |>
            filter(!label %in% rownames(annotated_tips)) |>
            pivot_longer(
                names_to = "Attribute", values_to = "Score", cols = 2:last_col()
            ) |>
            left_join(tree_data, by = "label") |>
            filter(!is.na(NCBI_ID)) |>
            filter(!NCBI_ID %in% unique(dat$NCBI_ID)) |>
            mutate(
                Evidence = 'asr',
                Confidence_in_curation = NA,
                Attribute_source = NA,
                Frequency = scores2Freq(Score)
            )
        if (attr_type == "multistate-intersection") {
            res <- bind_rows(asr_df, filter(dat, !is.na(Evidence))) |>
                mutate(
                    Attribute_group = attr_grp,
                    Attribute_type = attr_type
                ) |>
                rename(Attribute_value = Attribute) |>
                mutate(Attribute = Attribute_group) |>
                relocate(
                    NCBI_ID, Taxon_name, Rank, Attribute, Attribute_value,
                    Frequency, Score,
                    Attribute_source, Confidence_in_curation
                ) |>
                select(-label)
        } else if (attr_type == "binary") {
            res <- bind_rows(asr_df, filter(dat, !is.na(Evidence))) |>
                mutate(
                    Attribute_group = attr_grp,
                    Attribute_type = attr_type
                ) |>
                separate(
                    col = "Attribute", into = c("Attribute", "Attribute_value"),
                    sep = "--"
                ) |>
                relocate(
                    NCBI_ID, Taxon_name, Rank, Attribute, Attribute_value,
                    Frequency, Score,
                    Attribute_source, Confidence_in_curation
                ) |>
                select(-label)
        }
        propagated[[i]] <- res
        names(propagated)[i] <- names(dat_ready)[i]
    } else if (attr_type %in% c("range")) { # Provided as a range by physiologies in bugphyzz
        ## Propagation numeric ####
        annotated_tips <- fdat |>
            filter(!is.na(Attribute_value)) |>
            select(tip_label, Attribute_value) |>
            column_to_rownames(var = "tip_label") |>
            as.matrix()
        no_annotated_tips <- fdat |>
            filter(is.na(Attribute_value)) |>
            select(tip_label, Attribute_value) |>
            column_to_rownames(var = "tip_label") |>
            as.matrix()
        input_mat <- rbind(annotated_tips, no_annotated_tips)
        input_vct <- input_mat[tree$tip.label,, drop = TRUE]
        tim <- system.time({
            asr <- hsp_squared_change_parsimony(
                tree = tree, tip_states = input_vct, weighted = TRUE,
                check_input = TRUE
            )
        })
        log_print(tim)
        asr_df <- data.frame(
            label = c(tree$tip.label, tree$node.label),
            Attribute_value = asr$states
        ) |>
            filter(!grepl("^n\\d+$", label)) |>
            filter(!label %in% rownames(annotated_tips))
        nsti <- getNsti(tree = tree, annotated_tip_labels = rownames(annotated_tips))
        predicted_dat <- left_join(
            asr_df, nsti, by = c("label" = "tip_label")) |>
            left_join(tree_data, by = "label") |>
            select(-label) |>
            filter(!is.na(NCBI_ID)) |>
            mutate(
                Attribute_source = NA,
                Confidence_in_curation = NA,
                Frequency = "unknown",
                Score = NA,
                Evidence = "asr",
            ) |>
            filter(!NCBI_ID %in% dat$NCBI_ID)
        res <- bind_rows(dat, predicted_dat) |>
            mutate(
                Attribute_group = attr_grp,
                Attribute = attr_grp,
                Attribute_type = attr_type,
            ) |>
            relocate(NCBI_ID, Taxon_name, Rank, Attribute, Attribute_value)
        propagated[[i]] <- res
        names(propagated)[i] <- names(dat_ready)[i]
    }
}






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
## This would allow to use the same AUC-ROC method for comparison
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

data_ready <- map(categorical, ~ {
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
data_discarded <- names(categorical)[which(!names(categorical) %in% names(data_ready))]
if (length(data_discarded) > 0) {
    data_discarded <- paste0(', paste0(data_discarded, collapse = ', '), ')
    msg <- paste0(
        "The following physiologies were discarded because they didn't have",
        " any valid Attribute: ", data_discarded, '.'
    )
log_print(msg, blank_after = TRUE)
}


## This chunk is to make sure that all NCBI_IDs have valid taxon names
for (i in seq_along(data_ready)) {
    x <- data_ready[[i]]
    set1 <- filter(x,!is.na(NCBI_ID))
    set2 <- filter(x, is.na(NCBI_ID))
    set1$Rank <- checkRank(set1$NCBI_ID)
    set1$Taxon_name <- checkTaxonName(set1$NCBI_ID)
    data_ready[[i]] <- bind_rows(set1, set2)
    log_print(paste0("Finished checking valid taxon names for: ", names(data_ready)[i]),
              blank_after = TRUE)
}

## For those numeric attributes converted to categorical attributes,
## propagate them with their ranges and units.
data_ready <- map(data_ready, ~ {
    if ('Attribute_range' %in% colnames(.x)) {
        .x$Attribute <- paste0(.x$Attribute, ' ', .x$Attribute_range)
        .x$Attribute <- sub('\\)$', '', .x$Attribute)
        .x$Attribute <- paste0(.x$Attribute, ' ', .x$Unit, ')')
        .x$Attribute <- sub(' \\)', ')', .x$Attribute)
    }
    return(.x)
})

# Following code chunk produces the following warning 21 times - is this harmless?
# There were 21 warnings (use warnings() to see them)
# > warnings()
# Warning messages:
#     1: There was 1 warning in `dplyr::mutate()`.
# ℹ In argument: `Evidence = forcats::fct_relevel(Evidence, "asr")`.
for (i in seq_along(data_ready)) {
    data_ready[[i]] <- tryCatch(
        output <- data_ready[[i]] |>
            prepareDataForPropagation() |>
            mergeOriginalAndEarlyASR() |>
            group_by(NCBI_ID) |>
            mutate(Score = Score / sum(Score)) |>
            ungroup() |>
            distinct()
    )
    log_print(paste0(
        "Finished preparing data for propagation: ",
        names(data_ready)[i]
    ),
    blank_after = TRUE)
}

lgl <- map_lgl(data_ready, is_error)
if (any(lgl)) {
    discard_names <- names(data_ready)[which(lgl)]
    msg <- paste0(
        "The following physiologies will be discared for some error during",
        " preparation for propagation: ",
        paste0("'", paste0(discard_names, collapse = ', '), ".'")
    )
    log_print(msg, blank_after = TRUE)
    data_ready <- discard(data_ready, is_error)
}

## Propagation with taxPPro
data('tree_list')
tree <- as.Node(tree_list)
propagated <- as.list(character(length(data_ready)))
names(propagated) = names(data_ready)


for (i in seq_along(data_ready)){
    p <- system.time({
    x <- data_ready[[i]]
    msg <- unique(x$Attribute_group)
    msg <- paste('Propagating', msg, 'at', Sys.time())
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

    propagated[[i]] <- final_table
    })
    print(p)
}


rm(categorical, data_ready, data_tree_tbl, data_with_values, empty_df,
   final_table, input_tbl, l, output, phys, range, range_cat, set1, set2,
   tree_list, x, all_node_names, attr_grp, attr_type, attrs, binaries,
   data_discarded, exclude_phys, i, lf, lgl, missing_node_names,
   more_valid_attributes, msg, msg_len, p, phys_names, tree, valid_attributes)

## Create header
## Create a header for both the dump files and the gmt files.
header <- paste0("# bugphyzz ", Sys.Date(),
                 ", License: Creative Commons Attribution 4.0 International",
                 ", URL: https://waldronlab.io/bugphyzz\n")
cat(header, file = 'full_dump_with_0.csv')

dropcols <- c("Attribute_value", "Parent_name", "Parent_rank", "Parent_NCBI_ID",
    "Strain", "Genome_ID", "Accession_ID")

for (i in seq_along(propagated)) {
    log_print(paste("Dumping", names(propagated)[i], "to file with zeros"), blank_after = TRUE)
    propagated[[i]] <-
        propagated[[i]][, !colnames(propagated[[i]]) %in% dropcols]
    ## Some code for dividing the Attribute column into
    ## Attribute range
    propagated[[i]]$Attribute_range = ifelse(
        test = grepl('\\(', propagated[[i]]$Attribute),
        yes = sub('^.*(\\(.*)$', '\\1', propagated[[i]]$Attribute),
        no = NA
    )
    propagated[[i]]$Attribute = sub('\\(.*$', '', propagated[[i]]$Attribute)
    propagated[[i]]$Attribute = stringr::str_squish(propagated[[i]]$Attribute)
    write.table(
        x = propagated[[i]],
        file = 'full_dump_with_0.csv',
        quote = TRUE,
        sep = ",",
        row.names = FALSE,
        append = TRUE,
        col.names = identical(i, 1L)
    )
}

system2('pbzip2', args = list('-f', 'full_dump_with_0.csv'))

cat(header, file = 'full_dump.csv')
for (i in seq_along(propagated)) {
    log_print(paste("Dumping", names(propagated)[i], "to file without zeros"), blank_after = TRUE)

    total_scores <- propagated[[i]] |> group_by(Attribute_group, NCBI_ID) |>
            reframe(Total_score = sum(Score))
    taxids_above_0 <- total_scores |>
        filter(Total_score > 0) |>
        pull(NCBI_ID) |>
        unique()
    propagated[[i]] <- propagated[[i]] |>
        filter(NCBI_ID %in% taxids_above_0)
    rm(taxids_above_0, total_scores)
    write.table(
        x = propagated[[i]], file = 'full_dump.csv', quote = TRUE, sep = ",",
        row.names = FALSE,
        append = TRUE, col.names = identical(i, 1L)
    )
}

system2('pbzip2', args = list('full_dump.csv'))

## Export gmt files
log_print('Writing GMT files...')
ranks <- c('genus', 'strain', 'species', 'mixed')
tax_id_types <- c('Taxon_name', 'NCBI_ID')


# helper function to add a header line to an already written dump or GMT file
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
            'bugphyzz-', ranks[i], '-', tax_id_types[j], '.gmt'
        )
        for (k in seq_along(propagated)){
            log_print(paste("rank:", ranks[i], "/ tax_id_type:", tax_id_types[j], "/ physiology:", names(propagated)[k]), blank_after = TRUE)
            sig <- getBugphyzzSignatures(
                df = propagated[[k]], tax.level = ranks[i], tax.id.type = tax_id_types[j]
            )
            bugsigdbr::writeGMT(sigs = sig, gmt.file = gmt_file, append = TRUE)
        }
        addHeader(header, gmt_file)
        counter <- counter + 1
    }
}

si <- sessioninfo::session_info()
log_print(si, blank_after = TRUE)

log_close()
