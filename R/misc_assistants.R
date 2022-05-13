#' @title Scaling and Centering of Rows in Matrix-like Objects
#' 
#' @description Extends `base::scale` to center and/or scale the rows of a numeric matrix.
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
quosure_helper <- function(.data, quosure_list) {
  for (i in seq_along(quosure_list)) {
    assays(.data)[[names(quosure_list)[i]]] <- rlang::eval_tidy(quosure_list[[i]],
                                                         data = assays(.data) %>% 
                                                           as.list())
  }
  .data
}