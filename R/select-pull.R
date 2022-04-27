#' @title Select `SummarizedExperiment`s from a `MultiAssayExperiment`
#' 
#' @description 
#' 
#' @param .data A `MultiAssayExperiment`.
#' @param ... The names of the `SummarizedExperiment`s to be kept.
#' 
#' @return A `MultiAssayExperiment` containing the selected `SummarizedExperiment`s.
#' 
#' @export
select <- function(.data, ...) {
  UseMethod("select")
}

#' @rdname select
#' @export
select.MultiAssayExperiment <- function(.data, ...) {
  # Convert `...` into a vector of strings of experiment names
  experiment_vector <- rlang::quos(...) %>% map(rlang::quo_text) %>% unlist()
  return(.data[, , experiment_vector]) 
}