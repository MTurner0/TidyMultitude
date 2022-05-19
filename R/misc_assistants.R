#' @title Scaling and Centering of Rows in Matrix-like Objects
#'
#' @description Extends `base::scale` to center and/or scale the rows of a
#'   numeric matrix.
#'
#' @param x a numeric matrix(like object), such as an assay.
#' @param center a logical value indicating whether rows should be centered.
#' @param scale a logical value indicating whether rows should be scaled.
#'
#' @export
scale_rowwise <- function(x, center = TRUE, scale = TRUE) {
  
  # Will only center if specified
  if (center) {
    x <- sweep(x, 1, apply(x, 1, mean))
  }
  
  # Will only scale if specified
  if (scale) {
    scales <- apply(x, 1, sd)
    for(i in 1:nrow(x)) { x[i, ] <- x[i, ]/scales[i] }
  }
  return(x)
}

#' @export
tidy_colData_helper <- function(.data, FUN, list_of_args) {
  
  # First element of list_of_args will be the word "list" -- replace with
  # colData
  list_of_args[[1]] <- colData(.data) %>% 
    as_tibble() %>% 
    
    # Add columns of indices that will be used to subset assay columns
    # Use a name that is unlikely to appear in colData
    dplyr::mutate(id_helper_QjWTNFtWmSBc8XS = 1:nrow(.))
  
  # Transform colData with specified function
  modded_colData <- do.call(FUN, args = list_of_args)
  
  # Subset columns of assay data by rows of colData
  modded_assay_list <- purrr::map(.x = as.list(assays(.data)),
                               ~ .x[, modded_colData$id_helper_QjWTNFtWmSBc8XS])
  names(modded_assay_list) <- assays(.data) %>% 
    names()
  
  SummarizedExperiment(assays = modded_assay_list, 
                       # Remove indexing column
                       colData = modded_colData %>% 
                         dplyr::select(-id_helper_QjWTNFtWmSBc8XS),
                       rowData = rowData(.data))
}

#' @export
quosure_helper <- function(.data, quosure_list, drop_unused = FALSE) {
  for (i in seq_along(quosure_list)) {
    assays(.data)[[names(quosure_list)[i]]] <- rlang::eval_tidy(quosure_list[[i]],
                                                         data = assays(.data) %>% 
                                                           as.list())
  }
  # Provides the transmute functionality
  if (drop_unused) {
    assays(.data) <- assays(.data)[names(quosure_list)]
  }
  .data
}

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
    as_tibble() %>% 
    dplyr::mutate(id_helper_mrujhmAqKlLj9cJ = 1:n())
  
  mutated_colData2 <- colData(.data[[exp2_name]]) %>% 
    as_tibble() %>% 
    dplyr::mutate(id_helper_QFSIxvNAMbN8ltI = 1:n())
  
  # Make new experiment 1
  exp1_indices <- merge(mutated_colData1, mutated_colData2, by = by) %>% 
    dplyr::select(id_helper_mrujhmAqKlLj9cJ) %>% dplyr::pull()
  
  exp1_assays <- purrr::map(.x = as.list(assays(.data[[exp1_name]])),
                            ~ .x[, exp1_indices])
  names(exp1_assays) <- assays(.data[[exp1_name]]) %>% 
    names()
  
  new_exp1 <- SummarizedExperiment(assays = exp1_assays,
                                   colData = colData(.data[[exp1_name]])[exp1_indices, ],
                                   rowData = rowData(.data[[exp1_name]]))
  
  # Make new experiment 2
  exp2_indices <- merge(mutated_colData1, mutated_colData2, by = by) %>% 
    dplyr::select(id_helper_QFSIxvNAMbN8ltI) %>% dplyr::pull()
  
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
