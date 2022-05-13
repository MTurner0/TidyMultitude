#' @title Create, modify, and delete assays
#'
#' @description `mutate()` adds new assays and preserves existing ones.
#'   `transmute()` adds new assays and drops existing ones. New assays overwrite
#'   existing assays of the same name. 
#'
#' @export
mutate <- function(.data, ...) {
  UseMethod("mutate")
}

#' @rdname mutate
#' @export
mutate.MultiAssayExperiment <- function(.data, experiment, ...) {
  mutate_quosures <- dplyr:::dplyr_quosures(...)
  experiment_name <- paste0(substitute(experiment))
  .data[[experiment_name]] <- quosure_helper(.data[[experiment_name]], 
                                             mutate_quosures)
  .data
}

#' @rdname mutate
#' @export
mutate.SummarizedExperiment <- function(.data, ...) {
  mutate_quosures <- dplyr:::dplyr_quosures(...)
  quosure_helper(.data, quosure_list = mutate_quosures)
}

#' @rdname mutate
#' @export
transmute <- function(.data, experiment, ...) {
  UseMethod("transmute")
}

#' @rdname mutate
#' @export
transmute.MultiAssayExperiment <- function(.data, experiment, ...) {
  mutate_quosures <- dplyr:::dplyr_quosures(...)
  experiment_name <- paste0(substitute(experiment))
  .data[[experiment_name]] <- quosure_helper(.data[[experiment_name]], 
                                             mutate_quosures,
                                             drop_unused = TRUE)
  .data
}

#' @rdname mutate
#' @export
transmute.SummarizedExperiment <- function(.data, ...) {
  mutate_quosures <- dplyr:::dplyr_quosures(...)
  quosure_helper(.data, quosure_list = mutate_quosures, drop_unused = TRUE)
}