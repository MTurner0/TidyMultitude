#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Documentation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Subset a `MultiAssayExperiment` by experiments.
#'
#' @description Select experiments in a `MultiAssayExperiment`.
#'
#' @param .data A `MultiAssayExperiment`.
#' @param ... The names of the experiments to be kept, unquoted.
#'
#' @return A `MultiAssayExperiment` containing the selected
#'   `SummarizedExperiment`s. The output has the following properties:
#'   \itemize{\item Attributes of selected experiments (e.g. `colData`,
#'   `rowData`, and `assay`s for `SummarizedExperiment` objects) are not
#'   affected. \item Output experiments are a subset of input experiments,
#'   potentially with a different order.}
#'
#' @seealso [dplyr::select()]
#'
#' @examples
#' # Return a MultiAssayExperiment object containing only 
#' # the SummarizedExperiment cyto
#' mae %>% select(cyto)
#' 
#' @name select
#'   
NULL

#' @title Extract a single experiment.
#'
#' @description `pull()` is similar to `[[]]`.
#'
#' @param .data A `MultiAssayExperiment` containing at least one experiment.
#' @param  var A variable specified as: \itemize{ \item a literal experiment
#'   name; \item a positive integer, giving the position of the experiment
#'   counting from the top; \item a negative integer, giving the position of the
#'   experiment counting from the bottom. } The default returns the last
#'   experiment. If there is only one experiment (e.g. if [select()]
#'   has been used to select a single experiment from the
#'   `MultiAssayExperiment`), then this will be chosen.
#'
#' @return An experiment.
#'
#' @seealso [dplyr::pull()]
#'
#' @examples
#' names(mae)
#' # [1] "phy16S" "cyto"
#' 
#' # Both of the lines below would return only the cyto SummarizedExperiment
#' mae %>% pull(cyto)
#' mae %>% pull()
#' 
#' # Both of the lines below would return only the phy16S SummarizedExperiment
#' mae %>% pull(phy16S)
#' mae %>% pull(1)
#'
#' @name pull
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom dplyr select
#' @importFrom rlang quos quo_text
#' @importFrom purrr map
#' @export
select.MultiAssayExperiment <- function(.data, ...) {
  # Convert `...` into a vector of strings of experiment names
  experiment_vector <- quos(...) %>%
    map(quo_text) %>%
    unlist()
  for (i in seq_along(experiment_vector)) {
    if (!(experiment_vector[i] %in% names(.data))) {
      stop(paste0("Invalid experiment name: ", experiment_vector[i]))
    }
  }
  return(suppressWarnings(.data[, , experiment_vector]))
}

#' @importFrom dplyr pull
#' @export
pull.MultiAssayExperiment <- function(.data, var = -1) {
  name <- substitute(var)
  if(is.symbol(name)) {
    var <- paste0(name)
    stopifnot("Invalid experiment name." = var %in% names(.data))
  }
  if(is.numeric(var) && var < 0) {
    
    # Convert negative index to positive index
    var <- experiments(.data) %>% length() + 1 + var
  }
  return(.data[[var]])
}

#' @noRd
#' @export
dplyr::select

#' @noRd
#' @export
dplyr::pull
