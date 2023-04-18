#' Subset a data.frame by threshold
#'
#' @param df A data.frame.
#' @param thr A named list of thresholds.
#'
#' @return A data.frame
#' @export
#'
subsetByThreshold <- function(df, thr) {
    lower <- thr$lower
    upper <- thr$upper
    if (is.na(lower)) {
        output <- df[which(df$Attribute_value_max < upper),]
    } else if (is.na(upper)) {
        output <- df[which(df$Attribute_value_min >= lower),]
    } else {
        output <- df[which(df$Attribute_value_min >= lower & df$Attribute_value_max < upper),]
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
        x$Attribute <- paste0(attr_grp, ':', thr_names[i])
        x$Attribute_value_min <- NULL
        x$Attribute_value_max <- NULL
        x$Attribute_value <- TRUE
        subsets[[i]] <- x
    }
    output <- dplyr::bind_rows(subsets)
    return(output)
}
