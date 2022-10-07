#' @title Keep biological units based on matching factors in `colData`
#'
#' @description `MultiAssayExperiment::intersectColumns()` only works if the
#'   experiments have the same IDs within the `sampleMap`. This function will
#'   use `colData` features to match biological units (currently only supports
#'   two experiments).
#'
#' @param .data A MultiAssayExperiment object.
#'
#' @param experiment1 A SummarizedExperiment within `.data`.
#'
#' @param experiment2 Another SummarizedExperiment within `.data`.
#'
#' @param by A vector of `colData` features (i.e. column names) to match units
#'   by.
#'   
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
#' @export
intersect_colData <- function(.data, experiment1, experiment2, by) {
  
  # If the experiments are not specified, use the first two experiments by
  # default
  if (missing(experiment1) & missing(experiment2)) {
    exp1_name <- names(.data)[1]
    exp2_name <- names(.data)[2]
  } else {
    exp1_name <- paste0(substitute(experiment1))
    exp2_name <- paste0(substitute(experiment2))
  }
  
  mutated_colData1 <- colData(.data[[exp1_name]]) %>% 
    dplyr::as_tibble() %>% 
    mutate(id_helper_mrujhmAqKlLj9cJ = 1:nrow(.))
  
  mutated_colData2 <- colData(.data[[exp2_name]]) %>% 
    dplyr::as_tibble() %>% 
    mutate(id_helper_QFSIxvNAMbN8ltI = 1:nrow(.))
  
  # Make new experiment 1
  exp1_indices <- merge(mutated_colData1, mutated_colData2, by = by) %>% 
    select(id_helper_mrujhmAqKlLj9cJ) %>% pull()
  
  exp1_assays <- purrr::map(.x = as.list(assays(.data[[exp1_name]])),
                            ~ .x[, exp1_indices])
  names(exp1_assays) <- assays(.data[[exp1_name]]) %>% 
    names()
  
  new_exp1 <- SummarizedExperiment(assays = exp1_assays,
                                   colData = colData(.data[[exp1_name]])[exp1_indices, ],
                                   rowData = rowData(.data[[exp1_name]]))
  
  # Make new experiment 2
  exp2_indices <- merge(mutated_colData1, mutated_colData2, by = by) %>% 
    select(id_helper_QFSIxvNAMbN8ltI) %>% pull()
  
  exp2_assays <- purrr::map(.x = as.list(assays(.data[[exp2_name]])),
                            ~ .x[, exp2_indices])
  names(exp2_assays) <- assays(.data[[exp2_name]]) %>% 
    names()
  
  new_exp2 <- SummarizedExperiment(assays = exp2_assays,
                                   colData = colData(.data[[exp2_name]])[exp2_indices, ],
                                   rowData = rowData(.data[[exp2_name]]))
  
  # Build MAE
  exp_list <- list(new_exp1, new_exp2)
  names(exp_list) <- c(exp1_name, exp2_name)
  
  MultiAssayExperiment(experiments = exp_list)
}