#' @title Arrange columns of a `MultiAssayExperiment` or `SummarizedExperiment`
#'
#' @description Allow columns of assays within a `MultiAssayExperiment` or
#' `SummarizedExperiment` to be rearranged by the observations of their
#' `colData`.
#'
#' @param .data Either a `MultiAssayExperiment` or `SummarizedExperiment`
#'   object.
#'
#' @export
arrange_colData <- function(.data, ...) {
  UseMethod("arrange_colData")
}

#' @rdname arrange_colData
#' @export 
arrange_colData.MultiAssayExperiment <- function(.data, ...) {
  arrange_args <- as.list(substitute(list(...)))
  reordered_mae <- .data
  for(i in seq_along(experiments(.data))) {
    
    # Arrange colData for each SE in the MAE
    reordered_mae[[i]] <- TidyMultitude:::tidy_colData_helper(.data[[i]], dplyr::arrange, arrange_args)
  }
  reordered_mae
}

#' @rdname arrange_colData
#' @export
arrange_colData.SummarizedExperiment <- function(.data, ...) {
  arrange_args <- as.list(substitute(list(...)))
  TidyMultitude:::tidy_colData_helper(.data, dplyr::arrange, arrange_args)
}

#' @title Arrange rows of a `MultiAssayExperiment` or `SummarizedExperiment`
#'
#' @param .data Either a `MultiAssayExperiment` or `SummarizedExperiment`
#'   object.
#'
#' @export
arrange_rowData <- function(.data, ...) {
  UseMethod("arrange_rowData")
}

#' @rdname arrange_rowData
#' @export 
arrange_rowData.MultiAssayExperiment <- function(.data, ...) {
  arrange_args <- as.list(substitute(list(...)))
  reordered_mae <- .data
  for(i in seq_along(experiments(.data))) {
    
    # Arrange colData for each SE in the MAE
    reordered_mae[[i]] <- TidyMultitude:::tidy_rowData_helper(.data[[i]], dplyr::arrange, arrange_args)
  }
  reordered_mae
}

#' @rdname arrange_rowData
#' @export
arrange_rowData.SummarizedExperiment <- function(.data, ...) {
  arrange_args <- as.list(substitute(list(...)))
  TidyMultitude:::tidy_rowData_helper(.data, dplyr::arrange, arrange_args)
}