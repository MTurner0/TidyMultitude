#' @title Remove rows containing all zeroes
#' 
#' @description Sequencing data often produces large but sparse matrices. To save space and time, it may be worthwhile to remove empty columns and/or rows.
#' 
#' @param .data Either a `MultiAssayExperiment` or a `SummarizedExperiment`.
#' @param counts If `TRUE`, empty rows can be checked using `sum(x) > 1`, which is faster than `any(x != 0)`.
#' 
#' @export
trim_empty_rows <- function(.data, experiment, counts = TRUE) {
  UseMethod("trim_empty_rows")
}

#' @rdname trim_empty_rows
#' @export
trim_empty_rows.MultiAssayExperiment <- function(.data, experiment, counts = TRUE) {
  exp_name <- deparse(substitute(experiment))
  .data[[exp_name]] <- trim_empty_rows.SummarizedExperiment(.data[[exp_name]],
                                                            counts = counts)
  return(.data)
}

#' @rdname trim_empty_rows
#' @export
trim_empty_rows.SummarizedExperiment <- function(.data, counts = TRUE) {
  # Create a vector indicating which rows of the matrix and rowData are nonzero
  if (counts) {
    nonempty_indices <- apply(.data@assays@data@listData[[1]], 1, sum) > 0
  } else {
    nonzero <- function(x) { any(x != 0)  }
    nonempty_indices <- apply(.data@assays@data@listData[[1]], 1, nonzero)
  }
  new_rowdata <- rowData(.data)[nonempty_indices, ]
  new_assay <- .data@assays@data@listData[[1]][nonempty_indices, ]
  new_assay_list <- new_assay %>% list(.)
  names(new_assay_list) <- names(.data@assays)
  trimmed_experiment <- SummarizedExperiment(assays = new_assay_list,
                                             colData = .data@colData,
                                             rowData = new_rowdata)
  return(trimmed_experiment)
}