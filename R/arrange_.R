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
#' `SummarizedExperiment` to be rearranged by the observations of their
#' `rowData`.
#'
#' @param .data Either a `MultiAssayExperiment` or `SummarizedExperiment`
#'   object.
#' @param ... Columns in `rowData` to arrange by.
#'
#' @export
arrange_rowData <- function(.data, ...) {
  UseMethod("arrange_rowData")
}

#' @rdname arrange_rowData
#' @importFrom dplyr arrange
#' @export
arrange_rowData.MultiAssayExperiment <- function(.data, ...) {
  arrange_args <- as.list(substitute(list(...)))
  reordered_mae <- .data
  for (i in seq_along(experiments(.data))) {

    # Arrange colData for each SE in the MAE
    reordered_mae[[i]] <- tidy_rowData_helper(.data[[i]], arrange, arrange_args)
  }
  reordered_mae
}

#' @rdname arrange_rowData
#' @importFrom dplyr arrange
#' @export
arrange_rowData.SummarizedExperiment <- function(.data, ...) {
  arrange_args <- as.list(substitute(list(...)))
  tidy_rowData_helper(.data, arrange, arrange_args)
}
