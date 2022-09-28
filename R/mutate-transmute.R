#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Documentation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Create, modify, and delete assays or `SummarizedExperiments`
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
mutate.MultiAssayExperiment <- function(.data, experiment, ..., .keep = c("all", "none")) {
  if (missing(.keep)) {
    .keep <- "all"
  }
  
  mutate_quosures <- dplyr:::dplyr_quosures(...)
  
  # Operate on SEs within MAE
  if (missing(experiment)) {
    return(quosure_helper(.data, mutate_quosures, .keep))
  }
  
  # Operate on assays within SE
  else {
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
  if(missing(.keep)) {
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
  }
  
  # Operate on assays within SE
  else {
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

#' @export
dplyr::mutate

#' @export
dplyr::transmute