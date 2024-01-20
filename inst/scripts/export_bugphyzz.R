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
        "host-associated",
        'sphingolipid producing',
        'butyrate producing'
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
system.time({
    phys <- physiologies(phys_names, full_source = FALSE)
})

dat_ready <- vector("list", length(phys))
for (i in seq_along(dat_ready)) {
    message(names(phys)[i])
    dat_ready[[i]] <- getDataReady(filterData(phys[[i]]))
    names(dat_ready)[i] <- names(phys)[i]
}

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
                pivot_longer(
                    names_to = "Attribute", values_to = "Score", cols = 2:last_col()
                ) |>
                left_join(tree_data, by = "label", relationship = "many-to-many") |>
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

            nsti <- getNsti(tree = tree, annotated_tip_labels = rownames(annotated_tips))

            predicted_dat <- left_join(
                asr_df, nsti, by = c("label" = "tip_label")) |>
                left_join(tree_data, by = "label", relationship = "many-to-many") |>
                select(-label) |>
                filter(!is.na(NCBI_ID)) |>
                mutate(
                    Attribute_source = NA,
                    Confidence_in_curation = NA,
                    Frequency = "unknown",
                    Score = NA,
                    Evidence = "asr"
                ) |>
                filter(!NCBI_ID %in% dat$NCBI_ID)

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
## Some are not well defined (might need further review) or redundant (free-linvn and host-associated|FALSE)

## TODO refactor the code below to reduce repetition
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
ar <- physiologies("antimicrobial resistance")[[1]]
arl <- split(ar, ar$Attribute) |>
    purrr::map(~ mutate(.x, Attribute_type = "binary"))
arready <- purrr::map(arl, ~ getDataReady(filterData(.x))) |>
    bind_rows() |>
    separate(col = "Attribute", into = c("Attribute", "Attribute_value"), sep = "--") |>
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
            log_print(paste("rank:", ranks[i], "/ tax_id_type:", tax_id_types[j], "/ physiology:", names(propagated)[k]), blank_after = TRUE)
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
