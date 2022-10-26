#' @title A `MultiAssayExperiment` Object
#'
#' @description 26 patients were randomly assigned to a `treatment` or
#'   `placebo`.
#'
#'   On each of 5 visits, stool samples were obtained for 16S rRNA-seq (of 10
#'   bacterial species). However, 5 of the sequencing results failed quality
#'   control and had to be discarded. The results are stored in the `counts`
#'   assay in the `phy16S` `SummarizedExperiment`.
#'
#'   On visits 1, 3, and 5, the concentrations of 8 cytokines were measured.
#'   These results are stored in the `cyto_conc` assay in the `cyto`
#'   `SummarizedExperiment`.
#'
#'   (The data within this `MultiAssayExperiment` is entirely fabricated. It
#'   exists as a "toy" example of a `MultiAssayExperiment`.)
#'
#' @examples
#' print(mae)
"mae"
