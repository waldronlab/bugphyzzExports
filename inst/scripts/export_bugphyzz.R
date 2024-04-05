
## This script contains the code for generating
## the text files imported by the `bugphyzz::importBugphyzz` function.
## Check the bugphyzz package on waldronlab/bugphyzz or Bioconductor.

suppressMessages({
    library(logr)
    library(bugphyzz)
    library(taxPPro)
    library(rlang)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(phytools)
    library(castor)
    library(purrr)
})

logfile <- "log_file"
lf <- log_open(logfile, logdir = FALSE, compact = TRUE, show_notes = FALSE)

attributes_by_type <- list(
    numeric = c(
        "growth temperature"
        # "coding genes",
        # "genome size",
        # "length",
        # "width",
        # "optimal ph",
        # "mutation rate per site per generation",
        # "mutation rate per site per year"
    ),
    binary = c(
        # "animal pathogen",
        # "antimicrobial sensitivity",
        # "biofilm forming",
        # "extreme environment",
        # "health associated",
        # "hydrogen gas producing",
        # "lactate producing",
        "motility"
        # "plant pathogenicity",
        # "spore formation",
        # "host-associated",
        # 'sphingolipid producing',
        # 'butyrate producing'
    ),
    multistate = c(
        "aerophilicity"
        # "gram stain",
        # "biosafety level",
        # "COGEM pathogenicity rating",
        # "shape",
        # "spore shape",
        # "arrangement",
        # "hemolysis"
    )
)

ltp <- ltp()
tree <- ltp$tree
tip_data <- ltp$tip_data
node_data <- ltp$node_data

tree_data <- bind_rows(
    select(
        as_tibble(tip_data), label = tip_label, NCBI_ID, Taxon_name, Rank,
        taxid
    ),
    select(
        as_tibble(node_data), label = node_label, NCBI_ID, Taxon_name, Rank,
        taxid
    )
)

phys_names <- unlist(attributes_by_type, use.names = FALSE)
msg <- paste0('"', paste0(phys_names, collapse = ', '), '"')
msg_len <- length(phys_names)
msg <- paste('Importing', msg_len, 'physiologies from bugphyzz:', msg)
log_print(msg, blank_after = TRUE)
tim <- system.time({
    phys <- physiologies(phys_names, full_source = FALSE)
})
log_print(tim, blank_after = TRUE)

dat_ready <- vector("list", length(phys))
dat_ready_w <- vector("list", length(phys))

for (i in seq_along(dat_ready)) {
    ws <- list()
    withCallingHandlers(
       warning = function(w) ws[[length(ws) + 1]] <<- w$message,
       expr = dat_ready[[i]] <-  getDataReady(filterData(phys[[i]]))
    )
    dat_ready_w[[i]] <- ws
    names(dat_ready_w)[i] <- names(phys)[i]
    names(dat_ready)[i] <- names(phys)[i]
}
wng_msgs <- as.character(list_flatten(dat_ready_w))

if (length(wng_msgs) > 0) {
    need_update_ids <- grep("taxon IDs", wng_msgs, value = TRUE) |>
        sub("^.*d: ", "", x = _) |>
        paste0(collapse = ", ") |>
        strsplit(", ") |>
        flatten_chr() |>
        unique() |>
        sort()
    invalid_attrs <- grep("^Invalid attributes", wng_msgs, value = TRUE)

    log_print("The following NCBI IDs might need to be updated:")
    for (i in seq_along(need_update_ids)) {
        log_print(need_update_ids[i])
    }

    log_print("")
    log_print("The following attributes are not valid and have been dropped:")
    for (i in seq_along(invalid_attrs)) {
        log_print(sub(": .*", ":", invalid_attrs[i]))
        ias <- strsplit(sub("^.*: ", "", invalid_attrs[i]),"---")[[1]]
        for (j in seq_along(ias)) {
            log_print(ias[j])
        }
        log_print("")
    }
}

nodes_with_taxid <- grep("^n", tree$node.label, value = TRUE, invert = TRUE)

propagated <- vector('list', length(dat_ready))
for (i in seq_along(propagated)) {
    names(propagated)[i] <- names(dat_ready)[i]
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
        tip_data, dat, by = 'NCBI_ID',
        relationship = "many-to-many"
    )

    tim <- system.time({
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

            ltp_per <- floor(nrow(annotated_tips) / Ntip(tree)  * 100)
            if (ltp_per < 1) {
                msg <- paste0("Not enough data for ASR for ", attr_grp, ". Skipping ASR.")
                log_print(msg)

                if (attr_type == "multistate-union") {
                    dat <- dat |>
                        rename(Attribute_value = Attribute) |>
                        mutate(Attribute = Attribute_group) |>
                        relocate(
                            NCBI_ID, Taxon_name, Rank, Attribute, Attribute_value,
                            Frequency, Score,
                            Attribute_source, Confidence_in_curation
                        )
                } else if (attr_type == "binary") {
                    dat <- dat |>
                        separate(
                        col = "Attribute", into = c("Attribute", "Attribute_value"),
                        sep = "--"
                    ) |>
                        relocate(
                            NCBI_ID, Taxon_name, Rank, Attribute, Attribute_value,
                            Frequency, Score,
                            Attribute_source, Confidence_in_curation
                        )

                }
                propagated[[i]] <- dat
                next
            }

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

            fit <- fitMk(tree = tree, x = input_mat, model = "ER", pi = "equal")
            ace <- ancr(fit, tips = TRUE)
            output_mat <- ace$ace
            rownames(output_mat) <- c(tree$tip.label, tree$node.label)

            asr_df <- output_mat |>
                as.data.frame() |>
                rownames_to_column(var = "label") |>
                filter(label %in% tree_data$label) |>
                filter(!label %in% rownames(annotated_tips)) |>
                left_join(tree_data, by = "label", relationship = "many-to-many") |>
                group_by(NCBI_ID) |>
                slice_max(order_by = label, n = 1, with_ties = FALSE) |>
                ungroup() |>
                pivot_longer(
                    names_to = "Attribute", values_to = "Score", cols = colnames(output_mat)
                ) |>
                filter(!is.na(NCBI_ID)) |>
                filter(!NCBI_ID %in% unique(dat$NCBI_ID)) |>
                mutate(
                    Evidence = 'asr',
                    Confidence_in_curation = NA,
                    Attribute_source = NA,
                    Frequency = scores2Freq(Score)
                )

                # left_join(tree_data, by = "label", relationship = "many-to-many") |>

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

        } else if (attr_type %in% c("range")) {
            annotated_tips <- fdat |>
                filter(!is.na(Attribute_value)) |>
                select(tip_label, Attribute_value) |>
                column_to_rownames(var = "tip_label") |>
                as.matrix()

            ltp_per <- round(nrow(annotated_tips) / Ntip(tree)  * 100)
            if (ltp_per < 1) {
                msg <- paste0("Not enough data for ASR for ", attr_grp, ". Skipping ASR.")
                log_print(msg)
                propagated[[i]] <- dat |>
                    mutate(Attribute = Attribute_group) |>
                    relocate(
                        NCBI_ID, Taxon_name, Rank,
                        Attribute, Attribute_value,
                        Evidence, Frequency,
                        Attribute_source, Confidence_in_curation
                    )
                next
            }

            no_annotated_tips <- fdat |>
                filter(is.na(Attribute_value)) |>
                select(tip_label, Attribute_value) |>
                column_to_rownames(var = "tip_label") |>
                as.matrix()

            input_mat <- rbind(annotated_tips, no_annotated_tips)
            input_vct <- input_mat[tree$tip.label,, drop = TRUE]
            asr <- hsp_squared_change_parsimony(
                tree = tree, tip_states = input_vct, weighted = TRUE,
                check_input = TRUE
            )
            asr_df <- data.frame(
                label = c(tree$tip.label, tree$node.label),
                Attribute_value = asr$states
            ) |>
                # filter(!grepl("^n\\d+$", label)) |>
                filter(label %in% tree_data$label) |>
                filter(!label %in% rownames(annotated_tips))

            nsti <- getNsti(tree = tree, annotated_tip_labels = rownames(annotated_tips), nodes = nodes_with_taxid)

            predicted_dat <- left_join(
                asr_df, nsti, by = c("label" = "label")) |>
                left_join(tree_data, by = "label", relationship = "many-to-many") |>
                select(-label) |>
                filter(!is.na(NCBI_ID)) |>
                mutate(
                    Attribute_source = NA,
                    Confidence_in_curation = NA,
                    Frequency = nsti2Freq(nsti),
                    Score = freq2Scores2(Frequency),
                    # Frequency = "always",
                    # Score = 1,
                    Evidence = "asr"
                ) |>
                filter(!NCBI_ID %in% dat$NCBI_ID) |>
                group_by(NCBI_ID) |>
                slice_max(order_by = Attribute_value, n = 1, with_ties = FALSE) |>
                ungroup()

            res <- bind_rows(dat, predicted_dat) |>
                mutate(
                    Attribute_group = attr_grp,
                    Attribute = attr_grp,
                    Attribute_type = attr_type,
                ) |>
                relocate(
                    NCBI_ID, Taxon_name, Rank,
                    Attribute, Attribute_value,
                    Evidence, Frequency,
                    Attribute_source, Confidence_in_curation
                )

            propagated[[i]] <- res
        }
    })

    log_print(tim)
}

## The following physiologies were not propagated. Some reasons:
## Too few annotations (< 1%)
## Some are not well defined (might need further review)
## or redundant (free-living and host-associated|FALSE)

## TODO re-factor the code below to reduce repetition
# 'disease association',
# 'antimicrobial resistance'

## habitat
h <- physiologies("habitat")[[1]]
hl <- split(h, h$Attribute) |>
    purrr::map(~ mutate(.x, Attribute_type = "binary"))
hready <- purrr::map(hl, ~ getDataReady(filterData(.x))) |>
    bind_rows() |>
    separate(col = "Attribute", into = c("Attribute", "Attribute_value"), sep = "--") |>
    mutate(
        Attribute_value = Attribute,
        Attribute = Attribute_group,
        Attribute_type = "multistate-union"
    ) |>
    filter(!is.na(Evidence)) |>
    relocate(
        NCBI_ID, Taxon_name, Rank,
        Attribute, Attribute_value,
        Evidence, Frequency, Score,
        Attribute_source, Confidence_in_curation,
        Attribute_type
    )
propagated[["habitat"]] <- hready

## disease association
da <- physiologies("disease association")[[1]]
dal <- split(da, da$Attribute) |>
    purrr::map(~ mutate(.x, Attribute_type = "binary"))
daready <- purrr::map(dal, ~ getDataReady(filterData(.x))) |>
    bind_rows() |>
    separate(col = "Attribute", into = c("Attribute", "Attribute_value"), sep = "--") |>
    mutate(
        Attribute_value = Attribute,
        Attribute = Attribute_group,
        Attribute_type = "multistate-union"
    ) |>
    filter(!is.na(Evidence)) |>
    relocate(
        NCBI_ID, Taxon_name, Rank,
        Attribute, Attribute_value,
        Evidence, Frequency, Score,
        Attribute_source, Confidence_in_curation,
        Attribute_type
    )
propagated[["disease association"]] <- daready

## antimicrobial resistance
## TODO The changes below should be done on the bugphyzzWrangling repo and
## upload directly to the spreadsheet curation
ar <- physiologies("antimicrobial resistance")[[1]] |>
    mutate(Attribute = stringr::str_squish(tolower(Attribute))) |>
    mutate(Attribute = sub("not resistant", "sensitive", Attribute)) |>
    mutate(Attribute = sub("carbapenem-resistan", "resistance to", Attribute)) |>
    mutate(Attribute = sub("resistant", "resistance", Attribute)) |>
    mutate(Attribute = sub("tot", "to", Attribute)) |>
    mutate(Attribute_value = ifelse(grepl("sensitive", Attribute), FALSE, Attribute_value)) |>
    mutate(Attribute = sub("sensitive", "resistance", Attribute)) |>
    mutate(Attribute_type = "binary")

arl <- split(ar, ar$Attribute)
arready <- purrr::map(arl, ~ getDataReady(filterData(.x))) |>
    bind_rows() |>
    # separate(col = "Attribute", into = c("Attribute", "Attribute_value"), sep = "--") |>
    filter(!is.na(Evidence)) |>
    mutate(
        Attribute_value = Attribute,
        Attribute = Attribute_group,
        Attribute_type = "multistate-union"
    ) |>
    filter(!is.na(Evidence)) |>
    relocate(
        NCBI_ID, Taxon_name, Rank,
        Attribute, Attribute_value,
        Evidence, Frequency, Score,
        Attribute_source, Confidence_in_curation,
        Attribute_type
    ) |>
    mutate(
        Attribute_value = case_when(
            grepl("--TRUE$", Attribute_value) ~ sub("--TRUE$", "", Attribute_value),
            grepl("--FALSE$", Attribute_value) ~ sub("\\w+ to (.*)--FALSE$", "sensitive to \\1", Attribute_value)
        )
    )


propagated[["antimicrobial resistance"]] <- arready
final_set <- propagated |>
    purrr::map(~ {
        .x |>
            select(-taxid) |>
            mutate(NCBI_ID = sub("^\\w__", "", NCBI_ID)) |>
            filter(!is.na(Evidence)) |>
            select(-Attribute_group) |>
            relocate(
                NCBI_ID, Taxon_name, Rank,
                Attribute, Attribute_value,
                Evidence, Frequency, Score,
                Attribute_source, Confidence_in_curation,
                Attribute_type
            )
    })
multistate_data <- bind_rows(final_set[c(attributes_by_type$multistate, "habitat", "antimicrobial resistance", "disease association")])
binary_data <- bind_rows(final_set[attributes_by_type$binary])
numeric_data <- bind_rows(final_set[attributes_by_type$numeric])
numeric_data <- numeric_data |>
    mutate(Attribute_type = "numeric")

## Create a header for both the dump files and the gmt files.
header <- paste0("# bugphyzz ", Sys.Date(),
                 ", License: Creative Commons Attribution 4.0 International",
                 ", URL: https://waldronlab.io/bugphyzz\n")

cat(header, file = 'bugphyzz_multistate.csv')
write.table(
    x = multistate_data, file = "bugphyzz_multistate.csv", quote = TRUE,
    sep = ",", row.names = FALSE, append = TRUE, col.names = TRUE
)

cat(header, file = 'bugphyzz_binary.csv')
write.table(
    x = binary_data, file = "bugphyzz_binary.csv", quote = TRUE,
    sep = ",", row.names = FALSE, append = TRUE, col.names = TRUE
)

cat(header, file = 'bugphyzz_numeric.csv')
write.table(
    x = numeric_data, file = "bugphyzz_numeric.csv", quote = TRUE,
    sep = ",", row.names = FALSE, append = TRUE, col.names = TRUE
)


all_data <- list(multistate_data, binary_data, numeric_data)
all_data_list <- purrr::map(all_data, ~ split(.x, .x$Attribute)) |>
    purrr::list_flatten()

## Export gmt files
log_print('Writing GMT files...')
ranks <- c('genus', 'strain', 'species', 'mixed')
tax_id_types <- c('Taxon_name', 'NCBI_ID')

# helper function to add a header line to an already written dump or GMT file
addHeader <- function(header, out.file) {
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
        for (k in seq_along(all_data_list)){
            log_print(paste0("rank: ", ranks[i], "; tax_id_type: ", tax_id_types[j], "; physiology: ", names(propagated)[k]), blank_after = FALSE)
            sig <- makeSignatures(
                dat = all_data_list[[k]], tax_id_type = tax_id_types[j],
                tax_level = ranks[i]
            )
            if (!length(sig)) {
                next
            }
            bugsigdbr::writeGMT(sigs = sig, gmt.file = gmt_file, append = TRUE)
        }
        addHeader(header, gmt_file)
        counter <- counter + 1
    }
}

si <- sessioninfo::session_info()
log_print(si, blank_after = TRUE)
log_close()
