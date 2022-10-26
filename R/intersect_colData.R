#' @title Keep biological units based on matching factors in `colData`
#'
#' @description [MultiAssayExperiment::intersectColumns()] only works if the
#'   experiments have the same IDs within the `sampleMap`. This function will
#'   use `colData` features to match biological units.
#'
#' @param .data A `MultiAssayExperiment` object.
#'
#' @param by A vector of `colData` features (i.e. column names) to match units
#'   by. The column names must appear in the `colData` for each experiment.
#'
#' @return A `MultiAssayExperiment` whose experiments all contain `colData` with
#'   matching features.
#'
#' @examples
#' intersect_colData(mae, by = c("Patient", "Visit"))
#'
#' @importFrom purrr reduce
#' @export
intersect_colData <- function(.data, by) {
  
  if (!has_local_colData(.data)) {
    stop("Experiments don't have local colData. Try running map_colData.")
  }

  MultiAssayExperiment(
    experiments = reduce(
      .x = MAE_to_EL(.data),
      { function(x, y) intersect_ELs(x, y, by = by) }
    )
  )
}
