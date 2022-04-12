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

arrange_columns <- function(.data, by1, by2) {
  UseMethod("arrange_columns")
}

arrange_columns.MultiAssayExperiment <- function(.data, by1, by2) {
  var1 <- deparse(substitute(by1))
  var2 <- deparse(substitute(by2))
  z <- length(.data@ExperimentList)
  copypasta <- .data
  for (i in 1:z) { 
    copypasta[[i]] <- arrange_columns_helper(.data[[i]], {{ by1 }}, {{ by2 }}, var1, var2)
  }
  return(copypasta)
}

# A variant of arrange_columns.SE that is called from inside the arrange_columns.MAE function
# Must handle data masking differently than arrange_columns.SE, though the two could theoretically be merged

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

trim_empty_rows.SummarizedExperiment <- function(.data, counts = TRUE) {
  # Create a vector indicating which rows of the matrix and rowData are nonzero
  if (counts) {
    nonempty_indices <- apply(.data@assays@data@listData[[1]], 1, sum) > 0
  } else {
    nonempty_indices <- apply(.data@assays@data@listData[[1]], 1, nonzero)
  }
  new_rowdata <- rowData(.data)[nonempty_indices, ]
  new_assay <- .data@assays@data@listData[[1]][nonempty_indices, ]
  new_assay_list <- new_assay %>% list(.)
  names(new_assay_list) <- names(.data@assays)
  trimmed_experiment <- SummarizedExperiment(assays = new_assay_list,
                                             colData = .data@colData,
                                             rowData = new_rowdata)
  return(trimmed_experiment)
}


nonzero <- function(x) {
  any(x != 0)
}


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
