#' Subset a data.frame by threshold
#'
#' @param df A data.frame.
#' @param thr A named list of thresholds.
#'
#' @return A data.frame
#' @export
#'
subsetByThreshold <- function(df, thr) {
    lower <- thr[['lower']]
    upper <- thr[['upper']]
    if (is.na(lower)) {
        output <- df[which(df$Attribute_value_max < upper),]
        output$Attribute_range <- paste0('(<', upper, ')')
    } else if (is.na(upper)) {
        output <- df[which(df$Attribute_value_min > lower),]
        output$Attribute_range <- paste0('(>', lower, ')')
    } else {
        output <- df[which(df$Attribute_value_min >= lower & df$Attribute_value_max < upper),]
        output$Attribute_range <- paste0('(>=', lower, ', <', upper, ')')
    }
    return(output)
}

#' Convert a range data.frame to logical based on thresholds
#'
#' \code{rangeToLogicalThr}
#'
#' @param df A data.frame.
#' @param thresholds A list of thresholds.
#'
#' @return A data.frame.
#'
#' @export
#'
rangeToLogicalThr <- function(df, thresholds) {
    subsets <- vector('list', length(thresholds))
    thr_names <- names(thresholds)
    attr_grp <- unique(df$Attribute_group)
    for (i in seq_along(thresholds)) {
        names(subsets)[i] <- thr_names[i]
        x <- subsetByThreshold(df, thresholds[[i]])
        x$Attribute <- thr_names[i]
        x$Attribute_value_min <- NULL
        x$Attribute_value_max <- NULL
        x$Attribute_value <- TRUE
        x$Attribute_type <- 'logical'
        subsets[[i]] <- x
    }
    output <- dplyr::bind_rows(subsets) |>
        dplyr::distinct()
        # dplyr::mutate(Attribute = paste0(.data$Attribute, ' ', .data$Attribute_range))
    return(output)
}

#' THRESHOLDS
#'
#' \code{THRESHOLDS} this function returns a list of thresholds for numeric
#' attributes. This woud allow to convert them to categorical/logical in
#' order to use them in the propagation workflow.
#'
#' @return A named nested list.
#' @export
#'
THRESHOLDS <- function() {
    list(
        `coding genes` = list(
            `very small` = c(lower = NA, upper = 473),
            small = c(lower = 474, upper = 600),
            average = c(lower = 601, upper = 6000),
            `very large` = c(lower = 6001, upper = NA)
        ),
        `fatty acid` = list(
            `minimally present` = c(lower = NA, upper = 0.99),
            `lower concentration` = c(lower = 0.991, upper = 23.2),
            concentrated = c(lower = 23.21, upper = 55.5),
            `very concentrated` = c(lower = 55.51, upper = NA)
        ),
        `genome size` = list(
            small = c(lower = NA, upper = 490885),
            average = c(lower = 490886, upper = 998123),
            large = c(lower = 998124, upper = 6997434),
            `very large` = c(lower = 6997435, upper = NA)
        ),
        `growth temperature` = list(
            psychrophile = c(lower = NA, upper = 24.9),
            mesophile = c(lower = 25, upper = 45),
            thermophile = c(lower = 46, upper = 60),
            hyperthermophile = c(lower = 61, upper = NA)
        ),
        length = list(
            small = c(lower = NA, upper = 3.8),
            average = c(lower = 3.9, upper = 22),
            large = c(lower = 23, upper = 60),
            `very large` = c(lower = 61, upper = NA)
        ),
        `mutation rate per site per generation` = list(
            slow = c(lower = NA, upper = 2.92),
            medium = c(lower = 2.93, upper = 16),
            fast = c(lower = 17, upper = NA)
        ),
        `mutation rate per site per year` = list(
            slow = c(lower = NA, upper = 7.5),
            medium = c(lower = 7.6, upper = 20),
            `medium fast` = c(lower = 21, upper = 54.2),
            fast = c(lower = 54.3, upper = NA)
        ),
        `optimal ph` = list(
            acidic = c(lower = NA, upper = 6), # 5 --> 6
            neutral = c(lower = 6, upper = 8), # upper 7 --> 8
            alkaline = c(lower = 8, upper = 9.76), # upper 9.75 --> 9.76
            `very alkaline` = c(lower = 9.76, upper = NA)
        ),
        width = list(
            small = c(lower = NA, upper = 0.9),
            average = c(lower = 0.91, upper = 3.5),
            large = c(lower = 3.51, upper = 12),
            `very large` = c(lower = 13, upper = NA)
        )
    )
}

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
