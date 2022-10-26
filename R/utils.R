# Extends `base::scale` to center and/or scale the rows of a numeric matrix.
scale_rowwise <- function(x, center = TRUE, scale = TRUE) {

  # Will only center if specified
  if (center) {
    x <- sweep(x, 1, apply(x, 1, mean))
  }

  # Will only scale if specified
  if (scale) {
    scales <- apply(x, 1, sd)
    for (i in seq_len(nrow(x))) {
      x[i, ] <- x[i, ] / scales[i]
      }
  }
  return(x)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Helpers for colData and rowData operations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom dplyr as_tibble n
#' @importFrom purrr map
#' @importFrom SummarizedExperiment colData rowData
tidy_colData_helper <- function(.data, FUN, list_of_args) {

  # First element of list_of_args will be the word "list" -- replace with
  # colData
  list_of_args[[1]] <- colData(.data) %>%
    as_tibble() %>%

    # Add columns of indices that will be used to subset assay columns
    # Use a name that is unlikely to appear in colData
    mutate(id_helper_QjWTNFtWmSBc8XS = seq_len(n()))

  # Transform colData with specified function
  modded_colData <- do.call(FUN, args = list_of_args)

  # Subset columns of assay data by rows of colData
  modded_assay_list <- map(
    .x = as.list(assays(.data)),
    ~ .x[, modded_colData$id_helper_QjWTNFtWmSBc8XS]
    )
  names(modded_assay_list) <- assays(.data) %>%
    names()

  SummarizedExperiment(
    assays = modded_assay_list,
    colData = modded_colData %>% select(-id_helper_QjWTNFtWmSBc8XS),
    rowData = rowData(.data)
    )
}

#' @importFrom dplyr as_tibble n
#' @importFrom purrr map
#' @importFrom SummarizedExperiment colData rowData
tidy_rowData_helper <- function(.data, FUN, list_of_args) {

  # First element of list_of_args will be the word "list" -- replace with
  # rowData
  list_of_args[[1]] <- rowData(.data) %>%
    as_tibble() %>%

    # Add columns of indices that will be used to subset assay columns
    # Use a name that is unlikely to appear in colData
    mutate(id_helper_QjWTNFtWmSBc8XS = seq_len(n()))

  # Transform rowData with specified function
  modded_rowData <- do.call(FUN, args = list_of_args)

  # Subset columns of assay data by rows of colData
  modded_assay_list <- map(
    .x = as.list(assays(.data)),
    ~ .x[modded_rowData$id_helper_QjWTNFtWmSBc8XS, ]
    )

  names(modded_assay_list) <- assays(.data) %>%
    names()

  SummarizedExperiment(
    assays = modded_assay_list,
    rowData = modded_rowData %>% select(-id_helper_QjWTNFtWmSBc8XS),
    colData = colData(.data)
    )
}

# Checks whether a SummarizedExperiment has colData
#' @importFrom SummarizedExperiment colData
has_colData <- function(.data) {
  !(0 %in% (.data %>% colData() %>% dim()))
}

# Checks whether all Experiments in an MAE have local colData
has_local_colData <- function(.data) {
  suppressWarnings(
    all(
      lapply(experiments(.data), has_colData)
      )
    )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Helper for mutate/transmute
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

quosure_helper <- function(.data, quosure_list, .keep = "all") {
  UseMethod("quosure_helper")
}

quosure_helper.MultiAssayExperiment <- function(
    .data, quosure_list, .keep = "all"
    ) {
  exp_list <- experiments(.data)
  for (i in seq_along(quosure_list)) {
    exp_list[[names(quosure_list)[i]]] <- rlang::eval_tidy(
      quosure_list[[i]], data = exp_list %>% as.list()
    )
  }
  # Provides the transmute functionality
  if (.keep == "none") {
    exp_list <- exp_list[names(quosure_list)]
  }
  MultiAssayExperiment(experiments = exp_list)
}

quosure_helper.SummarizedExperiment <- function(
    .data, quosure_list, .keep = "all"
    ) {
  for (i in seq_along(quosure_list)) {
    assays(.data)[[names(quosure_list)[i]]] <- rlang::eval_tidy(
      quosure_list[[i]], data = assays(.data) %>% as.list()
    )
  }
  # Provides the transmute functionality
  if (.keep == "none") {
    assays(.data) <- assays(.data)[names(quosure_list)]
  }
  .data
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Helper for intersect_colData
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr as_tibble n
intersect_ELs <- function(EL1, EL2, by) {

  indices <- merge(colData(EL1[[1]]) %>%
                     as_tibble() %>%
                     mutate(id_helper_mrujhmAqKlLj9cJ = seq_len(dplyr::n())),
                   colData(EL2[[1]]) %>%
                     as_tibble() %>%
                     mutate(id_helper_QFSIxvNAMbN8ltI = seq_len(dplyr::n())),
                   by = by)

  c(
    lapply(
      EL1,
      {
        function(x) x[, pull(indices, id_helper_mrujhmAqKlLj9cJ)]
        }
    ),
    lapply(
      EL2,
      {
        function(x) x[, pull(indices, id_helper_QFSIxvNAMbN8ltI)]
        }
    )
  )
}

MAE_to_EL <- function(.data) {
  lapply(names(.data), {
    function(x) experiments(.data)[x]
    })
}
