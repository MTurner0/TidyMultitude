library(dplyr)
library(TidyMultitude)

exp <- curatedTCGAData::curatedTCGAData(
  "BRCA", "GISTIC*", "2.0.1", dry.run = FALSE
)

map_global_to_local <- function(.data) {
  exp_names <- experiments(.data) %>% names()
  new_exp_list <- list()
  for (i in seq_along(exp_names)) {
    suppressWarnings(new_exp_list[exp_names[i]] <- getWithColData(.data, i))
  }
  MultiAssayExperiment(experiments = new_exp_list)
}


has_colData <- function(.data) {
  !(0 %in% (.data %>% colData() %>% dim()))
}
