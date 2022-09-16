#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Documentation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Create, modify, and delete assays
#'
#' @description `mutate()` adds new assays and preserves existing ones.
#'   `transmute()` adds new assays and drops existing ones. New assays overwrite
#'   existing assays of the same name. 
#'
#' @seealso \code{\link[dplyr:mutate]{dplyr::mutate}}
#' @rdname mutate
#' @name mutate
#' 
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom dplyr mutate
#' @export
mutate.MultiAssayExperiment <- function(.data, experiment, ...) {
  mutate_quosures <- dplyr:::dplyr_quosures(...)
  experiment_name <- paste0(substitute(experiment))
  .data[[experiment_name]] <- TidyMultitude:::quosure_helper(.data[[experiment_name]], 
                                             mutate_quosures)
  .data
}

#' @importFrom dplyr mutate
#' @export
mutate.SummarizedExperiment <- function(.data, ...) {
  mutate_quosures <- dplyr:::dplyr_quosures(...)
  TidyMultitude:::quosure_helper(.data, quosure_list = mutate_quosures)
}

#' @importFrom dplyr transmute
#' @export
transmute.MultiAssayExperiment <- function(.data, experiment, ...) {
  mutate_quosures <- dplyr:::dplyr_quosures(...)
  experiment_name <- paste0(substitute(experiment))
  .data[[experiment_name]] <- TidyMultitude:::quosure_helper(.data[[experiment_name]], 
                                             mutate_quosures,
                                             drop_unused = TRUE)
  .data
}

#' @importFrom dplyr transmute
#' @export
transmute.SummarizedExperiment <- function(.data, ...) {
  mutate_quosures <- dplyr:::dplyr_quosures(...)
  TidyMultitude:::quosure_helper(.data, quosure_list = mutate_quosures, drop_unused = TRUE)
}

#' @export
dplyr::mutate

#' @export
dplyr::transmute