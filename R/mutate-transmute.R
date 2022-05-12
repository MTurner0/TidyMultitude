#' @title Create, modify, and delete assays
#'
#' @description `mutate()` adds new assays and preserves existing ones. New
#' assays overwrite existing assays of the same name.
#'
#' @export
mutate <- function(.data, ...) {
  UseMethod("mutate")
}

#' @export
mutate.MultiAssayExperiment <- function(.data, experiment, ...) {
  mutate_quosures <- dplyr:::dplyr_quosures(...)
  experiment_name <- paste0(substitute(experiment))
  .data[[experiment_name]] <- quosure_helper(.data[[experiment_name]], 
                                             mutate_quosures)
  .data
}

#' @export
mutate.SummarizedExperiment <- function(.data, ...) {
  mutate_quosures <- dplyr:::dplyr_quosures(...)
  quosure_helper(.data, quosure_list = mutate_quosures)
}