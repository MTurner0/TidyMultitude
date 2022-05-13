#Explore the possibility of the using/extending tibble SEs 

#tibble SE functionality
library(tidySummarizedExperiment)
library(dplyr)
library(MultiAssayExperiment)
#For needed data
library(HMP2Data)

#Make SE object
data("momspiCyto_mtx")
data("momspiCyto_samp")
momspiCyto <- SummarizedExperiment(assays = list(cyto_conc = momspiCyto_mtx),
                                   colData = momspiCyto_samp,
                                   rowData = data.frame(cytokine = 
                                                          rownames(momspiCyto_mtx)))

data("momspi16S_mtx")
data("momspi16S_samp")
data("momspi16S_tax")
momspi16S <- SummarizedExperiment(assays = list(counts = momspi16S_mtx),
                                  colData = momspi16S_samp,
                                  rowData = momspi16S_tax)

combo <- MultiAssayExperiment(experiments = list(phy16S = momspi16S,
                                                 cyto = momspiCyto))

#ORDER COLUMNS BY SUBJECT_ID AND VISIT

beta2 <- function(.data, ...) {
  arrange_args <- as.list(substitute(list(...)))
  reordered_mae <- .data
  for(i in seq_along(experiments(.data))) {
    reordered_mae[[i]] <- tidy_colData_helper(.data[[i]], arrange, arrange_args)
  }
  reordered_mae
}

# A variant of arrange_columns.SE that is called from inside the arrange_columns.MAE function
# Must handle data masking differently than arrange_columns.SE, though the two could theoretically be merged

beta <- function(.data, ...) {
  arrange_args <- as.list(substitute(list(...)))
  tidy_colData_helper(.data, arrange, arrange_args)
}

beta_helper <- function(.data, list_of_args) {
  # First element of list_of_args will be the word "list" -- replace with
  # colData
  list_of_args[[1]] <- colData(.data) %>% 
    as_tibble() %>% 
    # Add columns of indices that will be used to subset assay columns
    mutate(inds = 1:nrow(.))
  arranged_colData <- do.call(arrange, args = list_of_args)
  # Subset columns of assay data by rows of colData
  arranged_assay_list <- assay(.data)[, arranged_colData$inds] %>% 
    list(.)
  names(arranged_assay_list) <- names(.data@assays)
  arranged_se <- SummarizedExperiment(assays = arranged_assay_list,
                                      #Remove inds column
                                      colData = arranged_colData %>% 
                                        dplyr::select(-inds),
                                      rowData = rowData(.data))
  return(arranged_se)
}

tidy_colData_helper <- function(.data, FUN, list_of_args) {
  # First element of list_of_args will be the word "list" -- replace with
  # colData
  list_of_args[[1]] <- colData(.data) %>% 
    as_tibble() %>% 
    # Add columns of indices that will be used to subset assay columns
    # Use a name that is unlikely to appear in colData
    mutate(QjWTNFtWmSBc8XS = 1:nrow(.))
  modded_colData <- do.call(FUN, args = list_of_args)
  # Subset columns of assay data by rows of colData
  modded_assay_list <- assay(.data)[, modded_colData$QjWTNFtWmSBc8XS] %>% 
    list(.)
  names(modded_assay_list) <- assays(.data) %>% 
    names()
  SummarizedExperiment(assays = modded_assay_list, 
                       # Remove indexing column
                       colData = modded_colData %>% 
                         dplyr::select(-QjWTNFtWmSBc8XS),
                       rowData = rowData(.data))
}

arrange_columns_helper <- function(.data, by1, by2, var1, var2) {
  assay_name <- names(.data@assays)
  new_assay <- .data@assays@data@listData[[1]][order(.data@colData[[var1]], 
                                                     .data@colData[[var2]])]
  new_assay_list <- new_assay %>% 
    matrix(., nrow = nrow(rowData(.data)),
           ncol = nrow(.data@colData)) %>% list(.)
  names(new_assay_list) <- assay_name
  new_coldata <- .data@colData %>% as.data.frame() %>% arrange({{ by1 }}, {{ by2 }})
  reordered_experiment <- SummarizedExperiment(assays = new_assay_list,
                                               colData = new_coldata,
                                               rowData = rowData(.data))
  return(reordered_experiment)
}


arrange_columns.SummarizedExperiment <- function(.data, by1, by2) {
  var1 <- deparse(substitute(by1))
  var2 <- deparse(substitute(by2))
  assay_name <- names(.data@assays)
  new_assay <- .data@assays@data@listData[[1]][order(.data@colData[[var1]], 
                                                     .data@colData[[var2]])]
  new_assay_list <- new_assay %>% 
    matrix(., nrow = nrow(rowData(.data)),
            ncol = nrow(.data@colData)) %>% list(.)
  names(new_assay_list) <- assay_name
  new_coldata <- .data@colData %>% as.data.frame() %>% arrange({{ by1 }}, {{ by2 }})
  reordered_experiment <- SummarizedExperiment(assays = new_assay_list,
                                             colData = new_coldata,
                                             rowData = rowData(.data))
  return(reordered_experiment)
}


# DROP OTUs WITH ZERO OBSERVATIONS

trim_empty_rows <- function(.data, experiment) {
  UseMethod("trim_empty_rows")
}

trim_empty_rows.MultiAssayExperiment <- function(.data, experiment) {
  exp_name <- deparse(substitute(experiment))
  .data[[exp_name]] <- trim_empty_rows.SummarizedExperiment(.data[[exp_name]])
  return(.data)
}

delta <- function(.data, counts = TRUE) {
  # Create a vector identifying rows
  if (counts) {
    nonempty_indices <- (sapply(as.list(assays(.data)), rowSums) > 0) %>% 
      apply(., 1, any)
  } else {
    
    # Define internal functions
    matrix_nonzero <- function(x) {
      apply(x, 1, nonzero)
    }
    nonzero <- function(x) {
      any(x != 0)
    }
    
    nonempty_indices <- sapply(as.list(assays(.data)), matrix_nonzero) %>% 
      apply(., 1, any)
  }
  
  new_rowdata <- rowData(.data)[nonempty_indices, ]
  new_assay_list <- purrr::map(.x = as.list(assays(.data)),
                          ~ .x[nonempty_indices, ])
  names(new_assay_list) <- assays(.data) %>% names()
  SummarizedExperiment(assays = new_assay_list,
                       colData = colData(.data),
                       rowData = new_rowdata)
}

matrix_nonzero <- function(x) {
  apply(x, 1, nonzero)
}

nonzero <- function(x) {
  any(x != 0)
}



#Create scale_rowwise function
scale_rowwise <- function(x, center = TRUE, scale = TRUE) {
  if (center) {
    x <- sweep(x, 1, apply(x, 1, mean))
  }
  if (scale) {
    scales <- apply(x, 1, sd)
    for(i in 1:nrow(x)) { x[i, ] <- x[i, ]/scales[i] }
  }
  return(x)
}

#Select a subset of SEs from an MAE
select.MultiAssayExperiment <- function(.data, ...) {
  # Convert `...` into a vector of strings of experiment names
  experiment_vector <- rlang::quos(...) %>% map(rlang::quo_text) %>% unlist()
  return(.data[, , experiment_vector]) 
}

#Pull an SE from an MAE
pull.MultiAssayExperiment <- function(.data, var = -1) {
  if(is.numeric(var) & var < 0) {
    # Convert negative index to positive index
    var <- experiments(.data) %>% length() + 1 + var
  }
  return(.data[[var]])
}

#SE
filter_columns.MultiAssayExperiment <- function(.data, ...) {
  z <- experiments(.data) %>% length()
  filter_args <- as.list(substitute(list(...)))
  filtered_mae <- .data
  for(i in 1:z) {
    filtered_mae[[i]] <- filter_columns_helper(.data[[i]], filter_args)
  }
  return(filtered_mae)
}

filter_columns_depreciated <- function(.data, ...) {
  new_colData <- colData(.data) %>% 
    as_tibble() %>% 
    mutate(inds = 1:nrow(.)) %>% 
    filter(...)
  new_assay_list <- assay(.data)[, new_colData$inds] %>% list(.)
  names(new_assay_list) <- names(.data@assays)
  new_se <- SummarizedExperiment(assays = new_assay_list,
                                 colData = new_colData %>% select(-inds),
                                 rowData = rowData(.data))
  return(new_se)
}

gamma <- function(.data, ...) {
  filter_args <- as.list(substitute(list(...)))
  return(tidy_colData_helper(.data, filter, filter_args))
}

filter_columns_helper <- function(.data, list_of_args) {
  #First element of list_of_args will be the word "list" -- need to replace
  list_of_args[[1]] <- colData(.data) %>% 
    as_tibble() %>% 
    #Add columns of indices that will be used to subset assay columns
    mutate(inds = 1:nrow(.))
  filtered_colData <- do.call(filter, args = list_of_args)
  filtered_assay_list <- assay(.data)[, filtered_colData$inds] %>% list(.)
  names(filtered_assay_list) <- names(.data@assays)
  filtered_se <- SummarizedExperiment(assays = filtered_assay_list,
                                      #Remove inds column
                                      colData = filtered_colData %>% select(-inds),
                                      rowData = rowData(.data))
  return(filtered_se)
}

dummy <- function(...) {
  expr <- as.list(substitute(list(...)))
  return(expr)
}

dummy_2 <- function(hello = c("lorem", "ipsum", "sit", "amet")) {
  return(rlang::arg_match(hello))
}

dummy_inner <- function(.data, list_of_args) {
  #First element of list_of_args will be the word "list" -- need to replace
  list_of_args[[1]] <- colData(.data) %>% 
    as_tibble() %>% 
    #Add columns of indices that will be used to subset assay columns
    mutate(inds = 1:nrow(.))
  filtered_colData <- do.call(filter, args = list_of_args)
  return(filtered_colData)
}