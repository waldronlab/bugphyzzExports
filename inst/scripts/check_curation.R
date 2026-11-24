#!/usr/bin/env Rscript
## Check curation files (sources/sheets/*.csv) against bugphyzz's rules (#115).
## Run from the repository root:
##
##   Rscript inst/scripts/check_curation.R sources/sheets/shape.csv ...
##   Rscript inst/scripts/check_curation.R --all
##   Rscript inst/scripts/check_curation.R --base origin/devel [files]
##
## --all checks every file in sources/sheets/. --base <ref> reports only the
## problems on lines added or changed since the merge base with <ref>, so
## rows that were already there do not block a pull request; without files it
## checks the files of sources/sheets/ changed since then. Problems that
## belong to no line (unreadable file, unknown physiology) are always
## reported; header problems are on line 1.
##
## Prints one `file:line: message` per problem (line = where the record
## starts) and exits 1 if there is any. Needs bugphyzz installed: the rules
## come from its curation_template.tsv, attributes.tsv, attribute_synonyms.tsv,
## attribute_sources.tsv and spreadsheet_links.tsv, through its own functions.
## Tested by test_check_curation.R.

bp <- function(name) get(name, envir = asNamespace("bugphyzz"))

problem <- function(file, line, kind, message, end = line) {
    n <- length(message)
    data.frame(
        file = rep_len(file, n), line = as.integer(line), end = as.integer(end),
        kind = rep_len(kind, n), message = message, stringsAsFactors = FALSE
    )
}
no_problem <- function() problem(character(0), integer(0), "", character(0))

## Physiology and attribute type of a file, from its name ("gram_stain.csv" is
## "gram stain"), as .importSpreadsheets() takes them from
## spreadsheet_links.tsv. `type` is the name curation_template.tsv uses in
## attribute_types: logical for the discrete sheets, numeric or range.
## The fatty acid sheet is wide and read by .fattyAcidComposition().
sheet_type <- function(path) {
    phys <- gsub("_", " ", sub("\\.csv$", "", basename(path)))
    links <- utils::read.table(
        system.file("extdata", "spreadsheet_links.tsv", package = "bugphyzz"),
        header = TRUE, sep = "\t"
    )
    i <- match(phys, links$physiology)
    if (!is.na(i)) {
        at <- links$attribute_type[i]
        type <- if (at %in% bp(".discreteAttributeTypes")()) "logical" else at
        return(list(physiology = phys, type = type, wide = FALSE))
    }
    if (phys %in% bp(".customLinks")()$physiology)
        return(list(physiology = phys, type = "numeric", wide = TRUE))
    NULL
}

## Read a CSV the way read.csv() does, keeping every value as text and the
## line each record starts on. A quoted field can hold a line break, so a
## record can span lines. Returns the header, the records (padded or cut to
## the header's width), their first and last lines and field counts.
read_records <- function(path) {
    n <- utils::count.fields(
        path, sep = ",", quote = "\"", comment.char = "",
        blank.lines.skip = FALSE
    )
    if (!length(n)) stop("file is empty")
    last <- which(!is.na(n))
    ## count.fields() gives NA for every line of a record but its last
    if (!length(last) || max(last) < length(n) ||
        length(n) != length(readLines(path, warn = FALSE)))
        stop("a quoted field is never closed")
    first <- c(1L, utils::head(last, -1) + 1L)
    keep <- n[last] > 0                 # read.csv() skips blank lines
    first <- first[keep]; last <- last[keep]; n <- n[last]
    d <- utils::read.table(
        path, sep = ",", quote = "\"", comment.char = "", header = FALSE,
        fill = TRUE, col.names = paste0("V", seq_len(max(n))),
        colClasses = "character", na.strings = character(0),
        encoding = "UTF-8", stringsAsFactors = FALSE
    )
    if (nrow(d) != length(n)) stop("records and lines could not be matched")
    header <- unlist(d[1, seq_len(n[1])], use.names = FALSE)
    d <- d[-1, seq_len(n[1]), drop = FALSE]
    names(d) <- make.unique(header)
    list(
        header = header, data = d, first = first[-1], last = last[-1],
        nfields = n[-1]
    )
}

## Required columns a file must have, in order. curation_template.tsv
## describes the table after import; three of its required columns are made
## by the import and so are not in the files: Evidence and
## Confidence_in_curation are joined from attribute_sources.tsv by
## Attribute_source, and a range sheet has Attribute_value ("1-2", "<5"),
## which .modifyRange() splits into Attribute_value_min and _max. The wide
## fatty acid sheet has its values in one column per fatty acid.
file_columns <- function(type, wide = FALSE) {
    req <- bp(".requiredColumns")(type)
    req[req == "Attribute_value_max"] <- "Attribute_value"
    req <- setdiff(
        req, c("Attribute_value_min", "Evidence", "Confidence_in_curation")
    )
    if (wide) req <- setdiff(req, "Attribute_value")
    req
}

check_columns <- function(file, header, req) {
    miss <- setdiff(req, header)
    if (length(miss))
        return(problem(
            file, 1, "columns", paste0("required column ", miss, " is missing")
        ))
    if (!identical(header[seq_along(req)], req))
        return(problem(file, 1, "columns", paste0(
            "required columns are out of order: the file starts with ",
            paste(header[seq_along(req)], collapse = ", "), "; expected ",
            paste(req, collapse = ", ")
        )))
    NULL
}

squish <- function(x) stringr::str_squish(x)
## read.csv() reads the text NA as missing
missing_value <- function(x) !nzchar(x) | x == "NA"
missing_label <- function(col, x) paste0(col, ifelse(nzchar(x), " is NA", " is empty"))

## Rows whose values fail the valid_values regex of curation_template.tsv.
## Values are compared as physiologies() leaves them: whitespace squished,
## Frequency in lower case. An empty value fails only in a required column.
## Attribute and Attribute_value have their own checks below.
check_regex <- function(d, type) {
    tpl <- bp(".template")(d)
    tpl <- tpl[grepl(type, tpl$attribute_types), , drop = FALSE]
    tpl <- tpl[!tpl$column_name %in% c("Attribute", "Attribute_value"), ]
    out <- list()
    for (i in seq_len(nrow(tpl))) {
        col <- tpl$column_name[i]
        x <- squish(d[[col]])
        if (col %in% c("Frequency", "Evidence", "Confidence_in_curation"))
            x <- tolower(x)
        empty <- missing_value(x)
        bad <- !grepl(tpl$valid_values[i], ifelse(empty, "", x))
        if (tpl$requiredness[i] != "required") bad <- bad & !empty
        if (any(bad))
            out[[col]] <- data.frame(row = which(bad), message = ifelse(
                empty[bad], missing_label(col, x[bad]),
                paste0(col, " \"", x[bad], "\" does not match ",
                       tpl$valid_values[i])
            ))
    }
    do.call(rbind, out)
}

## Attribute_source must be in attribute_sources.tsv, or the row gets no
## Evidence and no Confidence_in_curation (.addSourceInfo()).
check_source <- function(d) {
    x <- squish(d$Attribute_source)
    u <- unique(x[!missing_value(x)])     # empty is reported by check_regex
    if (!length(u)) return(NULL)
    known <- bp(".addSourceInfo")(data.frame(Attribute_source = u))
    unknown <- known$Attribute_source[is.na(known$Evidence)]
    bad <- x %in% unknown
    if (!any(bad)) return(NULL)
    data.frame(row = which(bad), message = paste0(
        "Attribute_source \"", x[bad],
        "\" is not in bugphyzz's attribute_sources.tsv"
    ))
}

## Attribute_value as .processSheet() needs it: TRUE or FALSE in a discrete
## sheet (anything else is dropped), a number in a numeric sheet
## (.numericToRange() calls as.double(); the template regex applies too), and
## in a range sheet a value .modifyRange() keeps and can split.
check_value <- function(d, type) {
    x <- squish(d$Attribute_value)
    empty <- missing_value(x)
    if (type == "logical") {
        ok <- x %in% c("TRUE", "FALSE")
        what <- "is not TRUE or FALSE"
    } else if (type == "numeric") {
        tpl <- bp(".template")(d)
        rgx <- tpl$valid_values[tpl$column_name == "Attribute_value"]
        ok <- grepl(rgx, x) & !is.na(suppressWarnings(as.double(x)))
        what <- "is not a number"
    } else {
        u <- unique(x[!empty])
        r <- suppressWarnings(bp(".modifyRange")(
            data.frame(.value = u, Attribute_value = u)
        ))
        r <- r[!is.na(r$Attribute_value_min) & !is.na(r$Attribute_value_max), ]
        ok <- x %in% r$.value
        what <- "is not a number or a range such as 1-2, <5 or >5"
    }
    bad <- !ok
    if (!any(bad)) return(NULL)
    data.frame(row = which(bad), message = ifelse(
        empty[bad], missing_label("Attribute_value", x[bad]),
        paste0("Attribute_value \"", x[bad], "\" ", what)
    ))
}

## Attribute must be in attributes.tsv under the physiology's attribute_group
## (.validAttributes(), the match the export and taxPPro make), after
## .applySynonyms() has mapped it with attribute_synonyms.tsv. Synonyms apply
## to rows with Attribute_value TRUE only; a term mapped to nothing is dropped
## on purpose and is not a problem.
check_attribute <- function(d, physiology, type) {
    x <- squish(d$Attribute)
    empty <- missing_value(x)
    valid <- bp(".validAttributes")(physiology)
    ok <- x %in% valid
    if (type == "logical") {
        isTrue <- squish(d$Attribute_value) == "TRUE"
        u <- unique(x[isTrue & !empty])
        if (length(u)) {
            syn <- bp(".applySynonyms")(data.frame(
                Attribute_group = physiology, Attribute = u, .term = u,
                Attribute_value = TRUE
            ))
            mapped <- split(syn$Attribute, factor(syn$.term, levels = u))
            okTerm <- vapply(mapped, function(a) all(a %in% valid), logical(1))
            ok[isTrue] <- okTerm[x[isTrue]] %in% TRUE
        }
    }
    bad <- !ok | empty
    if (!any(bad)) return(NULL)
    known <- x %in% bp(".attributes")()
    data.frame(row = which(bad), message = ifelse(
        empty[bad], missing_label("Attribute", x[bad]),
        paste0(
            "Attribute \"", x[bad], "\" is not in bugphyzz's attributes.tsv",
            ifelse(known[bad], paste0(" under \"", physiology, "\""), "")
        )
    ))
}

## Every problem of one file: a data frame with file, line, end (last line of
## the record), kind (file, parse, columns, value, vocabulary) and message.
check_file <- function(path) {
    if (!file.exists(path))
        return(problem(path, NA, "file", "file not found"))
    st <- sheet_type(path)
    if (is.null(st))
        return(problem(path, NA, "file", paste0(
            "no physiology named \"",
            gsub("_", " ", sub("\\.csv$", "", basename(path))),
            "\" in bugphyzz's spreadsheet_links.tsv"
        )))
    r <- tryCatch(suppressWarnings(read_records(path)), error = function(e) e)
    if (inherits(r, "error"))
        return(problem(path, NA, "parse", paste0(
            "cannot be read as CSV: ", conditionMessage(r)
        )))
    out <- list()
    req <- file_columns(st$type, st$wide)
    out$columns <- check_columns(path, r$header, req)
    d <- r$data
    width <- length(r$header)
    wrong <- r$nfields != width
    out$parse <- problem(
        path, r$first[wrong], "parse",
        sprintf("%d fields, expected %d", r$nfields[wrong], width),
        r$last[wrong]
    )
    rows <- function(p, kind) {
        if (is.null(p)) return(NULL)
        p <- p[!wrong[p$row], , drop = FALSE]   # reported as parse problems
        problem(path, r$first[p$row], kind, p$message, r$last[p$row])
    }
    out$regex <- rows(check_regex(d, st$type), "value")
    if ("Attribute_source" %in% names(d))
        out$source <- rows(check_source(d), "value")
    if (!st$wide && all(c("Attribute", "Attribute_value") %in% names(d))) {
        out$value <- rows(check_value(d, st$type), "value")
        out$attribute <- rows(
            check_attribute(d, st$physiology, st$type), "vocabulary"
        )
    }
    out <- do.call(rbind, c(list(no_problem()), unname(out)))
    out <- out[order(out$line, out$kind), , drop = FALSE]
    rownames(out) <- NULL
    out
}

## Line numbers on the new side of a `git diff -U0`
diff_added_lines <- function(diff) {
    h <- regmatches(diff, regexec("^@@ -\\S+ \\+([0-9]+)(,([0-9]+))? @@", diff))
    h <- h[lengths(h) > 0]
    unlist(lapply(h, function(m) {
        n <- if (nzchar(m[4])) as.integer(m[4]) else 1L
        as.integer(m[2]) + seq_len(n) - 1L
    }))
}

git <- function(dir, args) {
    out <- suppressWarnings(
        system2("git", c("-C", shQuote(dir), args), stdout = TRUE, stderr = TRUE)
    )
    if (!is.null(attr(out, "status")))
        stop("git ", paste(args, collapse = " "), ": ",
             paste(out, collapse = "\n"), call. = FALSE)
    out
}

## Lines of `path` (working tree) added or changed since the merge base of
## `base` and HEAD
changed_lines <- function(path, base) {
    diff_added_lines(git(dirname(path), c(
        "diff", "-U0", "--merge-base", shQuote(base), "--",
        shQuote(basename(path))
    )))
}

## Keep the problems of records that touch a changed line, and those with no
## line
filter_changed <- function(problems, lines) {
    hit <- mapply(
        function(a, b) is.na(a) || any(lines >= a & lines <= b),
        problems$line, problems$end
    )
    problems[as.logical(hit), , drop = FALSE]
}

changed_files <- function(base, dir = "sources/sheets") {
    f <- git(dir, c(
        "diff", "--name-only", "--diff-filter=d", "--relative", "--merge-base",
        shQuote(base), "--", "."
    ))
    file.path(dir, f[grepl("^[^/]+\\.csv$", f)])
}

check_curation <- function(files, base = NULL) {
    out <- lapply(files, function(f) {
        p <- check_file(f)
        if (!is.null(base) && nrow(p)) p <- filter_changed(p, changed_lines(f, base))
        p
    })
    do.call(rbind, out)
}

report <- function(problems, nfiles) {
    p <- problems
    loc <- ifelse(is.na(p$line), p$file, paste0(p$file, ":", p$line))
    if (nrow(p)) cat(paste0(loc, ": ", p$message), sep = "\n")
    if (identical(Sys.getenv("GITHUB_ACTIONS"), "true") && nrow(p)) {
        esc <- function(x) gsub("\n", "%0A", gsub("\r", "%0D", gsub("%", "%25", x)))
        prop <- function(x) gsub(",", "%2C", gsub(":", "%3A", esc(x)))
        cat(paste0(
            "::error file=", prop(p$file),
            ifelse(is.na(p$line), "", paste0(",line=", p$line)),
            "::", esc(p$message)
        ), sep = "\n")
    }
    kinds <- table(factor(
        p$kind, levels = c("file", "parse", "columns", "value", "vocabulary")
    ))
    cat(
        nrow(p), " problem(s) in ", length(unique(p$file)), " of ", nfiles,
        " file(s) checked",
        if (nrow(p)) paste0(
            " (", paste0(names(kinds), " ", kinds, collapse = ", "), ")"
        ),
        "\n", sep = ""
    )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
    usage <- "usage: check_curation.R [--base <git ref>] (--all | file.csv ...)"
    base <- NULL
    i <- match("--base", args)
    if (!is.na(i)) {
        if (i == length(args)) stop(usage, call. = FALSE)
        base <- args[i + 1]
        args <- args[-c(i, i + 1)]
    }
    all <- "--all" %in% args
    files <- setdiff(args, "--all")
    if (any(grepl("^--", files))) stop(usage, call. = FALSE)
    if (all) {
        files <- list.files("sources/sheets", "\\.csv$", full.names = TRUE)
    } else if (!length(files)) {
        if (is.null(base)) stop(usage, call. = FALSE)
        files <- changed_files(base)
    }
    p <- check_curation(files, base)
    if (is.null(p)) p <- no_problem()
    report(p, length(files))
    if (nrow(p)) quit(status = 1)
}

if (sys.nframe() == 0L) main()
