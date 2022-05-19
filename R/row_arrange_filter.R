
#Filter Function
filter_rowData <- function(.data, ...) {
  UseMethod("filter_rowData")
}


filter_rowData.MultiAssayExperiment <- function(.data, ...) {
  filter_args <- as.list(substitute(list(...)))
  filtered_mae <- .data
  for(i in seq_along(experiments(.data))) {
    filtered_mae[[i]] <- tidy_rowData_helper(.data[[i]], filter, filter_args)
  }
  filtered_mae
}


filter_rowData.SummarizedExperiment <- function(.data, ...) {
  filter_args <- as.list(substitute(list(...)))
  tidy_rowData_helper(.data, filter, filter_args)
}

#arrange function
arrange_rowData <- function(.data, ...) {
  UseMethod("arrange_rowData")
}


arrange_rowData.MultiAssayExperiment <- function(.data, ...) {
  arrange_args <- as.list(substitute(list(...)))
  reordered_mae <- .data
  for(i in seq_along(experiments(.data))) {
    
    # Arrange colData for each SE in the MAE
    reordered_mae[[i]] <- tidy_rowData_helper(.data[[i]], arrange, arrange_args)
  }
  reordered_mae
}


arrange_rowData.SummarizedExperiment <- function(.data, ...) {
  arrange_args <- as.list(substitute(list(...)))
  tidy_rowData_helper(.data, arrange, arrange_args)
}



tidy_rowData_helper <- function(.data, FUN, list_of_args) {
  
  # First element of list_of_args will be the word "list" -- replace with
  # rowData
  list_of_args[[1]] <- rowData(.data) %>% 
    as_tibble() %>% 
    
    # Add columns of indices that will be used to subset assay columns
    # Use a name that is unlikely to appear in colData
    dplyr::mutate(id_helper_QjWTNFtWmSBc8XS = 1:nrow(.))
  
  # Transform rowData with specified function
  modded_rowData <- do.call(FUN, args = list_of_args)
  
  # Subset columns of assay data by rows of colData
  modded_assay_list <- purrr::map(.x = as.list(assays(.data)),
                                  ~ .x[modded_rowData$id_helper_QjWTNFtWmSBc8XS, ])
  names(modded_assay_list) <- assays(.data) %>% 
    names()
  
  SummarizedExperiment(assays = modded_assay_list, 
                       # Remove indexing column
                       rowData = modded_rowData %>% 
                         dplyr::select(-id_helper_QjWTNFtWmSBc8XS),
                       colData = colData(.data))
}