#' @title Arrange columns of a `MultiAssayExperiment` or `SummarizedExperiment`
#' 
#' @description 
#' Allow columns of assays within a `MultiAssayExperiment` or `SummarizedExperiment` to be rearranged by the observations of their `colData`.
#' 
#' @param .data Either a `MultiAssayExperiment` or `SummarizedExperiment` object.
#' 
#' @export
arrange_columns <- function(.data, by1, by2) {
  UseMethod("arrange_columns")
}

#' @rdname arrange_columns
#' @export 
arrange_columns.MultiAssayExperiment <- function(.data, by1, by2) {
  var1 <- deparse(substitute(by1))
  var2 <- deparse(substitute(by2))
  z <- length(.data@ExperimentList)
  reordered_mae <- .data
  for (i in 1:z) { 
    # Each experiment within the MAE has arrange_columns called on it
    reordered_mae[[i]] <- arrange_columns_helper(.data[[i]], {{ by1 }}, {{ by2 }}, var1, var2)
  }
  return(reordered_mae)
}

#' @rdname arrange_columns
#' @export
arrange_columns.SummarizedExperiment <- function(.data, by1, by2) {
  # Need text feature names for order()
  var1 <- deparse(substitute(by1))
  var2 <- deparse(substitute(by2))
  reordered_se <- arrange_columns_helper(.data, {{ by1 }}, {{ by2 }}, var1, var2)
  return(reordered_se)
}

#' @rdname arrange_columns
#' @export
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