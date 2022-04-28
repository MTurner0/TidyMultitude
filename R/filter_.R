#' @title Subset columns of `SummarizedExperiment`s by rows of `colData`.
#'
#' @description Subsets `colData`, retaining only rows that produce a value of
#'   `TRUE` for all conditions. The columns of `assay` data are then subset to
#'   match.
#'
#' @param .data A `SummarizedExperiment` or `MultiAssayExperiment`.
#' @param ... <\code{\link[dplyr]{dplyr_data_masking}}> Expressions that return
#'   a logical value, and are defined in terms of the variables in `.data`. If
#'   multiple expressions are included, they are combined with the `&` operator.
#'   Only rows for which all conditions evaluate to `TRUE` are kept.
#'
#' @export
filter_colData <- function(.data, ...) {
  UseMethod("filter_colData")
}

#' @rdname filter_colData
#' @export
filter_colData.MultiAssayExperiment <- function(.data, ...) {
  z <- experiments(.data) %>% length()
  filter_args <- as.list(substitute(list(...)))
  filtered_mae <- .data
  for(i in 1:z) {
    filtered_mae[[i]] <- filter_colData_helper(.data[[i]], filter_args)
  }
  return(filtered_mae)
}

#' @rdname filter_colData
#' @export
filter_colData.SummarizedExperiment <- function(.data, ...) {
  filter_args <- as.list(substitute(list(...)))
  return(filter_colData_helper(.data, filter_args))
}

#' @export
filter_colData_helper <- function(.data, list_of_args) {
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
                                      colData = filtered_colData %>% 
                                        dplyr::select(-inds),
                                      rowData = rowData(.data))
  return(filtered_se)
}