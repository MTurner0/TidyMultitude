#' @title Arrange columns of a `MultiAssayExperiment` or `SummarizedExperiment`
#'
#' @description Allow columns of assays within a `MultiAssayExperiment` or
#' `SummarizedExperiment` to be rearranged by the observations of their
#' `colData`.
#'
#' @param .data Either a `MultiAssayExperiment` or `SummarizedExperiment`
#'   object.
#' @param ... Columns in `colData` to arrange by.
#'
#' @returns An object of the same class as `.data`.
#'
#' @export
arrange_colData <- function(.data, ...) {
  UseMethod("arrange_colData")
}

#' @rdname arrange_colData
#' @importFrom dplyr arrange
#' @export
arrange_colData.MultiAssayExperiment <- function(.data, ...) {
  arrange_args <- as.list(substitute(list(...)))
  reordered_mae <- .data
  for (i in seq_along(experiments(.data))) {

    # Arrange colData for each SE in the MAE
    reordered_mae[[i]] <- tidy_colData_helper(.data[[i]], arrange, arrange_args)
  }
  reordered_mae
}

#' @rdname arrange_colData
#' @importFrom dplyr arrange
#' @export
arrange_colData.SummarizedExperiment <- function(.data, ...) {
  arrange_args <- as.list(substitute(list(...)))
  tidy_colData_helper(.data, arrange, arrange_args)
}

#' @title Arrange rows of a `MultiAssayExperiment` or `SummarizedExperiment`
#'
#' @description Allow rows of assays within a `MultiAssayExperiment` or
#'   `SummarizedExperiment` to be rearranged by the observations of their
#'   `rowData`. Unlike `colData`, `rowData` is experiment-specific, so only a
#'   single experiment can be operated on at a time.
#'
#' @param .data Either a `MultiAssayExperiment` or `SummarizedExperiment`
#'   object.
#' @param `experiment` If `.data` is a `MultiAssayExperiment`, the experiment to
#'   operate upon must be specified.
#' @param ... Columns in `rowData` to arrange by.
#'
#' @return An object of the same class as `.data`.
#'
#' @export
arrange_rowData <- function(.data, experiment, ...) {
  UseMethod("arrange_rowData")
}

#' @rdname arrange_rowData
#' @importFrom dplyr arrange
#' @export
arrange_rowData.MultiAssayExperiment <- function(.data, experiment, ...) {
  exp_name <- paste0(substitute(experiment))
  if (length(exp_name) > 1 || !(exp_name %in% names(.data))) {
    stop("Invalid experiment name.")
  }
  arrange_args <- as.list(substitute(list(...)))
  .data[[exp_name]] <- tidy_rowData_helper(
    .data[[exp_name]], arrange, arrange_args
  )
  .data
}

#' @rdname arrange_rowData
#' @importFrom dplyr arrange
#' @export
arrange_rowData.SummarizedExperiment <- function(.data, ...) {
  arrange_args <- as.list(substitute(list(...)))
  tidy_rowData_helper(.data, arrange, arrange_args)
}
