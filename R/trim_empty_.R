#' @title Remove rows containing all zeroes
#'
#' @description Sequencing data often produces large but sparse matrices. To
#'   save space and time, it may be worthwhile to remove empty columns and/or
#'   rows. Note that a row will only be removed if it is empty for all assays in
#'   a `SummarizedExperiment`.
#'
#' @param .data Either a `MultiAssayExperiment` or a `SummarizedExperiment`.
#' @param experiment If `.data` is a `MultiAssayExperiment`, choose which
#'   experiment to remove missing values from.
#' @param counts If `TRUE`, empty rows can be checked using `rowSums(x) > 0`,
#'   which is faster than `any(x != 0)`.
#' 
#' @return An object of the same class as `.data`.
#'
#' @export
trim_empty_rows <- function(.data, experiment, counts = TRUE) {
  UseMethod("trim_empty_rows")
}

#' @rdname trim_empty_rows
#' @export
trim_empty_rows.MultiAssayExperiment <- function(
    .data, experiment, counts = TRUE
    ) {
  exp_name <- paste0(substitute(experiment))
  if (!(exp_name %in% names(.data))) {
    stop("Invalid experiment name.")
  }
  .data[[exp_name]] <- trim_empty_rows(.data[[exp_name]], counts = counts)
  return(.data)
}

#' @rdname trim_empty_rows
#' @importFrom purrr map
#' @importFrom SummarizedExperiment colData rowData
#' @export
trim_empty_rows.SummarizedExperiment <- function(.data, counts = TRUE) {
  # Create a vector identifying rows
  if (counts) {
    nonempty_indices <- (sapply(as.list(assays(.data)), rowSums) > 0) %>%
      apply(., 1, any)
  } else {

    # Define internal functions
    matrix_nonzero <- function(x) {
      apply(x, 1, nonzero)
    }
    nonzero <- function(x) {
      any(x != 0)
    }

    nonempty_indices <- sapply(as.list(assays(.data)), matrix_nonzero) %>%
      apply(., 1, any)
  }

  new_rowdata <- rowData(.data)[nonempty_indices, ]

  new_assay_list <- map(
    .x = as.list(assays(.data)),
    ~ .x[nonempty_indices, ]
    )

  names(new_assay_list) <- assays(.data) %>%
    names()

  SummarizedExperiment(
    assays = new_assay_list,
    colData = colData(.data),
    rowData = new_rowdata
    )
}
