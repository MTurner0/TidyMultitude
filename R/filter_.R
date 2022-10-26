#' @title Subset columns of `SummarizedExperiment`s by rows of `colData`.
#'
#' @description Subsets `colData`, retaining only rows that produce a value of
#'   `TRUE` for all conditions. The columns of `assay` data are then subset to
#'   match. For a `MultiAssayExperiment`, this subsetting will be applied to all
#'   of its experiments.
#'
#' @param .data A `SummarizedExperiment` or `MultiAssayExperiment`.
#'
#' @param ... Expressions that return a logical value, and are defined in terms
#'   of the variables in `.data`. If multiple expressions are included, they are
#'   combined with the `&` operator. Only rows for which all conditions evaluate
#'   to `TRUE` are kept.
#'
#' @return An object of the same class as `.data`.
#'
#' @examples
#' # For all experiments, keep only sample data collected at the first visit.
#' filter_colData(mae, Visit == 1)
#'
#' @export
filter_colData <- function(.data, ...) {
  UseMethod("filter_colData")
}

#' @rdname filter_colData
#' @importFrom dplyr filter
#' @export
filter_colData.MultiAssayExperiment <- function(.data, ...) {

  if (!has_local_colData(.data)) {
    stop("Experiments don't have local colData. Try running map_colData.")
  }

  filter_args <- as.list(substitute(list(...)))
  filtered_mae <- .data
  for (i in seq_along(experiments(.data))) {
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
#' @description Subsets `rowData`, retaining only rows that produce a value of
#'   `TRUE` for all conditions. The rows of `assay` data are then subset to
#'   match. Unlike `colData`, `rowData` is experiment-specific, so only a single
#'   experiment can be operated on at a time.
#'
#' @param .data A `SummarizedExperiment` or `MultiAssayExperiment`.
#'
#' @param `experiment` If `.data` is a `MultiAssayExperiment`, the experiment to
#'   operate upon must be specified.
#'
#' @param ... Expressions that return a logical value, and are defined in terms
#'   of the variables in `.data`. If multiple expressions are included, they are
#'   combined with the `&` operator. Only rows for which all conditions evaluate
#'   to `TRUE` are kept.
#'
#' @return An object of the same class as `.data`.
#' 
#' @examples
#' mae %>%
#'   filter_rowData(phy16S, Class == "Bacilli") %>%
#'   filter_rowData(cyto, cytokine %in% c("IFN-g", "TNF-a"))
#'
#' @export
filter_rowData <- function(.data, experiment, ...) {
  UseMethod("filter_rowData")
}

#' @rdname filter_rowData
#' @importFrom dplyr filter
#' @export
filter_rowData.MultiAssayExperiment <- function(.data, experiment, ...) {
  exp_name <- paste0(substitute(experiment))
  if (length(exp_name) > 1 || !(exp_name %in% names(.data))) {
    stop("Invalid experiment name.")
  }
  filter_args <- as.list(substitute(list(...)))
  .data[[exp_name]] <- tidy_rowData_helper(
    .data[[exp_name]], filter, filter_args
    )
  .data
}

#' @rdname filter_rowData
#' @importFrom dplyr filter
#' @export
filter_rowData.SummarizedExperiment <- function(.data, ...) {
  filter_args <- as.list(substitute(list(...)))
  tidy_rowData_helper(.data, filter, filter_args)
}
