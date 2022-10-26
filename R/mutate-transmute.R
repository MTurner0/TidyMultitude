#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Documentation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Create, modify, and delete experiments or assays
#'
#' @description When operating on `MultiAssayExperiment`s, `mutate()` adds new
#'   experiments and preserves existing ones; `transmute()` adds new experiments
#'   and drops existing ones. New experiments overwrite existing experiments of
#'   the same name.
#'
#'   When operating on experiments, `mutate()` adds new assays and preserves
#'   existing ones; `transmute()` adds new assays and drops existing ones. New
#'   assays overwrite existing assays of the same name.
#'
#' @param `.data` A `MultiAssayExperiment` or `SummarizedExperiment`.
#'
#' @param `experiment` If `.data` is a `MultiAssayExperiment`, then this
#'   argument selects an experiment on which to perform assay manipulations.
#'
#' @param `...` <[dplyr::dplyr_data_masking()]>
#'   Name-value pairs. The name gives the name of the experiment or assay in the
#'   output. The value can be: \itemize{ \item When operating on
#'   `MultiAssayExperiment`s, an experiment compatible with the
#'   [MultiAssayExperiment::ExperimentList-class]
#'   container. \item When operating  on `SummarizedExperiment`s, an matrix with
#'   the dimensions of the `SummarizedExperiment`.}
#'
#' @seealso [dplyr::mutate()] [dplyr::transmute()]
#' @name mutate
#'   
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom dplyr mutate
#' @export
mutate.MultiAssayExperiment <- function(
    .data, experiment, ..., .keep = c("all", "none")
    ) {
  if (missing(.keep)) {
    .keep <- "all"
  }

  mutate_quosures <- dplyr:::dplyr_quosures(...)

  # Operate on SEs within MAE
  if (missing(experiment)) {
    return(quosure_helper(.data, mutate_quosures, .keep))
  } else { # Operate on assays within SE
    experiment_name <- paste0(substitute(experiment))
    .data[[experiment_name]] <- quosure_helper(
      .data[[experiment_name]], mutate_quosures, .keep
      )
    return(.data)
  }
}

#' @importFrom dplyr mutate
#' @export
mutate.SummarizedExperiment <- function(.data, ..., .keep = c("all", "none")) {
  if (missing(.keep)) {
    .keep <- "all"
  }

  mutate_quosures <- dplyr:::dplyr_quosures(...)
  quosure_helper(.data, mutate_quosures, .keep)
}

#' @importFrom dplyr transmute
#' @export
transmute.MultiAssayExperiment <- function(.data, experiment, ...) {
  mutate_quosures <- dplyr:::dplyr_quosures(...)

  # Operate on SEs within MAE
  if (missing(experiment)) {
    return(quosure_helper(.data, mutate_quosures, .keep = "none"))
  } else { # Operate on assays within SE
    experiment_name <- paste0(substitute(experiment))
    .data[[experiment_name]] <- quosure_helper(
      .data[[experiment_name]], mutate_quosures, .keep = "none"
    )
    return(.data)
  }

}

#' @importFrom dplyr transmute
#' @export
transmute.SummarizedExperiment <- function(.data, ...) {
  mutate_quosures <- dplyr:::dplyr_quosures(...)
  quosure_helper(.data, mutate_quosures, .keep = "none")
}

#' @noRd
#' @export
dplyr::mutate

#' @noRd
#' @export
dplyr::transmute
