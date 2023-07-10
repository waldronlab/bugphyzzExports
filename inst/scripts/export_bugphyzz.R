## Script to create bugphyzz exports dump files and signature files

library(bugphyzz)
library(taxPPro)
library(purrr)
library(rlang)
library(dplyr)
library(data.tree)
library(bugphyzzExports)

phys_names <- c(

    ## Categorical <<<<

    "acetate producing",
    "aerophilicity",
    "animal pathogen",
    "antimicrobial resistance",
    "antimicrobial sensitivity",
    # "arrangement",
    # "biofilm forming",
    # "biosafety level",
    "butyrate producing",
    # "COGEM pathogenicity rating",
    # "country",
    # "disease association",
    # "extreme environment",
    # "geographic location",
    "gram stain",
    # "growth medium",
    # "habitat",
    "health associated",
    # "hemolysis",
    "hydrogen gas producing",
    # "isolation site",
    "lactate producing",
    # "metabolite production",
    # "metabolite utilization",
    # "motility",
    "pathogenicity human",
    "plant pathogenicity",
    "shape",
    "sphingolipid producing",
    "spore formation",
    "spore shape",

    ## Ranges <<<<

    # "coding genes",
    # "genome size",
    "growth temperature",
    "halophily",
    "length",
    # "mutation rate per site per generation",
    # "mutation rates per site per year",
    "optimal ph",
    "width"
)
phys <- physiologies(phys_names, full_source = FALSE) |>
    map(removeAccessionAndGenomeID)

categorical <- keep(phys, ~ unique(.x$Attribute_type == 'logical'))
categorical$aerophilicity <- homogenizeAerophilicityAttributeNames(
    categorical$aerophilicity
)

range <- keep(phys, ~ unique(.x$Attribute_type == 'range'))
range <- range[which(names(range) %in% names(THRESHOLDS()))]
range_cat <- map2(range, names(range), ~ rangeToLogicalThr(.x, THRESHOLDS()[[.y]]))
categorical <- c(categorical, range_cat)

## Make sure that only valid attribute values are included.
fname <- system.file('extdata/attributes.tsv', package = 'bugphyzz')
attributes <- read.table(fname, header = TRUE, sep = '\t')


data <- map(categorical, ~ filter(.x, Attribute %in% unique(attributes$attribute)))
data <- keep(data, ~ nrow(.x) > 0)

data_ready <- vector('list', length(data))
for (i in seq_along(data_ready)) {
    message('Preparing ', names(data)[i], '.')
    names(data_ready)[i] <- names(data)[i]
    data_ready[[i]] <- tryCatch(
        error = function(e) e,
        {
            data[[i]] |>
                prepareDataForPropagation() |>
                mergeOriginalAndEarlyASR() |>
                group_by(NCBI_ID) |>
                mutate(Score = Score / sum(Score)) |>
                ungroup() |>
                distinct()
        }
    )
}
data_ready <- discard(data_ready, is_error)

message('>>>>>>> Propagating data ', Sys.time(), ' <<<<<<')
data('tree_list')
tree <- as.Node(tree_list)
propagated <- vector('list', length(data_ready))
for (i in seq_along(propagated)) {
    input_tbl <- data_ready[[i]] |>
        select(NCBI_ID, Attribute, Score, Evidence) |>
        distinct() |>
        tidyr::complete(
            NCBI_ID, Attribute, fill = list(Score = 0, Evidence = '')
        )
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
        relocate(NCBI_ID)

    attrs <- unique(data[[i]]$Attribute)
    all_node_names <- tree$Get(function(node) node$name, simplify = FALSE) |>
        unlist(recursive = TRUE, use.names = FALSE) |>
        unique()
    all_node_names <- all_node_names[which(!all_node_names %in% unique(data_tree_tbl$NCBI_ID))]
    empty_df <- data.frame(
        NCBI_ID = sort(rep(all_node_names, length(attrs))),
        Attribute = rep(attrs, length(all_node_names)),
        Score = 0,
        Evidence = NA
    )
    inferred_values <- bind_rows(data_tree_tbl, empty_df)

    other_ids <- data[[i]]$NCBI_ID[which(!phys$NCBI_ID %in% inferred_values$NCBI_ID)]
    other_ids <- unique(other_ids)
    other_phys <- filter(data[[i]], NCBI_ID %in% other_ids)
    final_table <- bind_rows(inferred_values, other_phys)
    propagated[[i]] <- final_table

    tree$Do(function(node) {
        node[['table']] <- NULL
    })
}
names(propagated) <- names(data_ready)

full_dump <- bind_rows(propagated)
full_dump$NCBI_ID <- sub('^[dpcofgst]__', '', full_dump$NCBI_ID)
full_dump$Attribute <- gsub(' ', '_', full_dump$Attribute)
full_dump_fname <- paste0("full_dump.csv.bz2")
unlink(full_dump_fname)
con <- bzfile(full_dump_fname, "w")
write.csv(x = full_dump, file = con, quote = TRUE, row.names = FALSE)
close(con)
