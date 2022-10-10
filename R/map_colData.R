#' @title Map global `colData` to local
#'
#' @description Ensures that each Experiment within a MultiAssayExperiment has 
#' `colData` relevant to its samples.
#'
#' @param .data A MultiAssayExperiment object.
#'   
#' @importFrom MultiAssayExperiment getWithColData
#' @export
map_colData <- function(.data) {
  UseMethod("map_colData")
}

#' @rdname map_colData
#' @export
map_colData.MultiAssayExperiment <- function(.data) {
  if(has_local_colData(.data)) {
    return(.data)
  }
  exp_names <- experiments(.data) %>% names()
  new_exp_list <- list()
  for (i in seq_along(exp_names)) {
    suppressWarnings(new_exp_list[exp_names[i]] <- getWithColData(.data, i))
  }
  MultiAssayExperiment(experiments = new_exp_list)
}