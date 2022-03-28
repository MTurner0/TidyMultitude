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

#Order columns by visit_number within subject_id

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
