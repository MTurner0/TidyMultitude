#' @title Subset columns of `SummarizedExperiment`s by rows of `colData`.
#'
#' @description Subsets `colData`, retaining only rows that produce a value of
#'   `TRUE` for all conditions. The columns of `assay` data are then subset to
#'   match.
#'
#' @param .data A `SummarizedExperiment` or `MultiAssayExperiment`.
#' 
#' @param ... Expressions that return
#'   a logical value, and are defined in terms of the variables in `.data`. If
#'   multiple expressions are included, they are combined with the `&` operator.
#'   Only rows for which all conditions evaluate to `TRUE` are kept.
#'
#' @export
filter_colData <- function(.data, ...) {
  UseMethod("filter_colData")
}

#' @rdname filter_colData
#' @importFrom dplyr filter
#' @export
filter_colData.MultiAssayExperiment <- function(.data, ...) {
  filter_args <- as.list(substitute(list(...)))
  filtered_mae <- .data
  for(i in seq_along(experiments(.data))) {
    filtered_mae[[i]] <- tidy_colData_helper(.data[[i]], filter, filter_args)
  }
  filtered_mae
}

#' @rdname filter_colData
#' @importFrom dplyr filter
#' @export
filter_colData.SummarizedExperiment <- function(.data, ...) {
  filter_args <- as.list(substitute(list(...)))
  tidy_colData_helper(.data, filter, filter_args)
}

#' @title Subset columns of `SummarizedExperiment`s by rows of `rowData`.
#'
#' @param .data A `SummarizedExperiment` or `MultiAssayExperiment`.
#' 
#' @param ... Expressions that return
#'   a logical value, and are defined in terms of the variables in `.data`. If
#'   multiple expressions are included, they are combined with the `&` operator.
#'   Only rows for which all conditions evaluate to `TRUE` are kept.
#'
#' @export
filter_rowData <- function(.data, ...) {
  UseMethod("filter_rowData")
}

#' @rdname filter_rowData
#' @importFrom dplyr filter
#' @export
filter_rowData.MultiAssayExperiment <- function(.data, ...) {
  filter_args <- as.list(substitute(list(...)))
  filtered_mae <- .data
  for(i in seq_along(experiments(.data))) {
    filtered_mae[[i]] <- tidy_rowData_helper(.data[[i]], filter, filter_args)
  }
  filtered_mae
}

#' @rdname filter_rowData
#' @importFrom dplyr filter
#' @export
filter_rowData.SummarizedExperiment <- function(.data, ...) {
  filter_args <- as.list(substitute(list(...)))
  tidy_rowData_helper(.data, filter, filter_args)
}