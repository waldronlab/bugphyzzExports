## Script to create bugphyzz exports dump files and signature files

library(logr)
library(bugphyzz)
library(taxPPro)
library(purrr)
library(rlang)
library(dplyr)
library(data.tree)
library(bugphyzzExports)
library(BiocParallel)
library(tidyr)

logfile <- "log_file"
lf <- log_open(logfile, logdir = FALSE, compact = TRUE, show_notes = FALSE)

phys_names <- c(

    ## Categorical <<<<
    "acetate producing",
    "aerophilicity",
    "animal pathogen",
    "antimicrobial resistance",
    "antimicrobial sensitivity",
    "arrangement",
    "biofilm forming",
    "biosafety level",
    "butyrate producing",
    # # "COGEM pathogenicity rating",
    # # "country",
    "disease association",
    "extreme environment",
    # # "geographic location",
    "gram stain",
    "growth medium",
    # # "habitat",
    "health associated",
    "hemolysis",
    "hydrogen gas producing",
    # # "isolation site",
    "lactate producing",
    # # "metabolite production",
    # # "metabolite utilization",
    "motility",
    "pathogenicity human",
    "plant pathogenicity",
    "shape",
    "sphingolipid producing",
    "spore formation",
    "spore shape",

    ## Ranges <<<<
    "coding genes",
    "genome size",
    "growth temperature",
    "halophily",
    "length",
    "mutation rate per site per generation",
    "mutation rates per site per year",
    "optimal ph",
    "width"
)

msg <- paste0('"', paste0(phys_names, collapse = ', '), '"')
msg_len <- length(phys_names)
msg <- paste('Importing', msg_len, 'physiologies from bugphyzz:', msg, '--', Sys.time())
log_print(msg, blank_after = TRUE)

phys <- physiologies(phys_names, full_source = FALSE)
categorical <- keep(phys, ~ unique(.x$Attribute_type) == 'logical')
categorical$aerophilicity <- homogenizeAerophilicityAttributeNames(
    categorical$aerophilicity
)
range <- keep(phys, ~ unique(.x$Attribute_type == 'range'))
range <- range[which(names(range) %in% names(THRESHOLDS()))]
range_cat <- map2(range, names(range), ~ rangeToLogicalThr(.x, THRESHOLDS()[[.y]]))
categorical <- c(categorical, range_cat)

log_print('Check that all attributes are valid. Invalid values will be printed:', blank_after = TRUE)
fname <- system.file('extdata/attributes.tsv', package = 'bugphyzz')
valid_attributes <- unique(read.table(fname, header = TRUE, sep = '\t')$attribute)
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
    filter(.x, Attribute %in% valid_attributes)
}) |>
    discard(~ !nrow(.x))

data_ready <- bplapply(data, BPPARAM = MulticoreParam(workers = 16), FUN = function(x) {
    tryCatch(
        error = function(e) e,
        {
            x |>
                prepareDataForPropagation() |>
                mergeOriginalAndEarlyASR() |>
                group_by(NCBI_ID) |>
                mutate(Score = Score / sum(Score)) |>
                ungroup() |>
                distinct()
        }
    )
})
data_ready <- discard(data_ready, is_error)

data('tree_list')
tree <- as.Node(tree_list)

nodes <- tree$Get(function(node) node$name)

propagated <- bplapply(X = data_ready, BPPARAM = MulticoreParam(workers = 16), FUN = function(x) {
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
    data_tree_tbl <- tree$Get(function(node) node[['table']], simplify = FALSE) |>
        purrr::discard(~ all(is.na(.x))) |>
        dplyr::bind_rows() |>
        dplyr::relocate(NCBI_ID)

    attrs <- unique(x$Attribute)
    all_node_names <- tree$Get(function(node) node$name, simplify = TRUE)
    all_node_names <- all_node_names[which(!all_node_names %in% unique(data_tree_tbl$NCBI_ID))]

    if (length(all_node_names > 0)) {
        empty_df <- data.frame(
            NCBI_ID = sort(rep(all_node_names, length(attrs))),
            Attribute = rep(attrs, length(all_node_names)),
            Score = 0,
            Evidence = NA
        )
        inferred_values <- bind_rows(data_tree_tbl, empty_df)
    } else {
        inferred_values <- data_tree_tbl
    }

    other_ids <- x$NCBI_ID[which(!phys$NCBI_ID %in% inferred_values$NCBI_ID)]
    other_ids <- unique(other_ids)
    other_phys <- filter(x, NCBI_ID %in% other_ids)
    final_table <- bind_rows(inferred_values, other_phys)

    tree$Do(function(node) {
        node[['table']] <- NULL
    })

    return(final_table)
})

full_dump <- bind_rows(propagated)
full_dump$NCBI_ID <- sub('^[dpcofgst]__', '', full_dump$NCBI_ID)
full_dump$Attribute <- gsub(' ', '_', full_dump$Attribute)

readr::write_csv(x = full_dump, file = "full_dump.csv.bz2", quote = 'needed', num_threads = 16)
log_close()
