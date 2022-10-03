#MUTATE: process and add assays

# Define my preprocessing function
prep <- function(data, scale = c(TRUE, FALSE)){
  # Default scaling to true
  if(missing(scale)){scale <- TRUE}
  
  # Note that the data is transposed -- I want to perform these steps on rows
  scaled <- scale(t(data), center = TRUE, scale = scale) # Center and scale
  tr <- sum(scaled * scaled)/(dim(scaled)[1]-1) # Fast trace computation
  processed <- scaled/sqrt(tr) # Normalize magnitude
  return(t(processed)) # Rotate back
}

# mutate.SummarizedExperiment(new = foo(old)) Make dumb version first: expects
# those arguments in order rather than having to dissect arg

mutate_simplest <- function(experiment, old, FUN, new) {
  indexer <- deparse(substitute(old))
  new_name <- deparse(substitute(new))
  assays(experiment)[[new_name]] <- assays(experiment)[[indexer]] %>% FUN()
  experiment
}

transmute_simplest <- function(experiment, old, FUN, new) {
  indexer <- deparse(substitute(old))
  new_name <- deparse(substitute(new))
  assays(experiment)[[new_name]] <- assays(experiment)[[indexer]] %>% FUN()
  assays(experiment) <- assays(experiment)[new_name] # Drop old assay
  experiment
}


eval_tidy(test_quosure$counts, data = assays(se) %>% as.list())

mutate.MultiAssayExperiment <- function(.data, experiment, ...) {
  mutate_quosures <- dplyr:::dplyr_quosures(...)
  experiment_name <- paste0(substitute(experiment))
  .data[[experiment_name]] <- quosure_helper(.data[[experiment_name]], 
                                             mutate_quosures)
  .data
}

mutate_se <- function(.data, ...) {
  my_quosures <- dplyr:::dplyr_quosures(...)
  for (i in seq_along(my_quosures)) {
    assays(.data)[[names(my_quosures)[i]]] <- rlang::eval_tidy(my_quosures[[i]], 
                                     data = assays(.data) %>% 
                              as.list())
  }
  .data
}

mutate.SummarizedExperiment <- function(.data, ...) {
  mutate_quosures <- dplyr:::dplyr_quosures(...)
  quosure_helper(.data, quosure_list = mutate_quosures)
}

# TRY ME

quosure_helper <- function(.data, quosure_list, drop_unused = FALSE) {
  for (i in seq_along(quosure_list)) {
    assays(.data)[[names(quosure_list)[i]]] <- eval_tidy(quosure_list[[i]],
                                                        data = assays(.data) %>% 
                                                          as.list())
  }
  if (drop_unused) {
    keep <- (assays(.data) %>% names())[assays(.data) %>% names() %in%
                                          names(quosure_list)]
    assays(.data) <- assays(.data)[keep]
  }
  .data
}

mutate_step1 <- function(experiment, ...) {
  
  new_name <- names(processed)
  call <- deparse(processed[[1]])
  assays(experiment)[[new_name]] <- assays(experiment)[[indexer]] %>% FUN()
  return(experiment)
}

my_quosure_test <- function(.data, ...) {
  res <- dplyr:::dplyr_quosures(...)
  inner_quosure_test(.data, res)
}

inner_quosure_test <- function(.data, args) {
  rlang::eval_tidy(args, data = .data)
}







# Needs to check whether `experiment`

bully <- function(.data, experiment, ...) {
  mutate_quosures <- dplyr:::dplyr_quosures(...)
  
  # Operate on SEs within MAE
  if (missing(experiment)) {
    print(mutate_quosures)
    .data <- dummy(.data, mutate_quosures)
    return(.data)
  }

  # Operate on assays within SE
  else {
    print(mutate_quosures)
    experiment_name <- paste0(substitute(experiment))
    .data[[experiment_name]] <- dummy(.data[[experiment_name]], 
                                      mutate_quosures)
    return(.data) 
  }
}

dummy <- function(.data, quosure_list, drop_unused = FALSE) {
  UseMethod("dummy")
}

dummy.MultiAssayExperiment <- function(.data, quosure_list, drop_unused = FALSE) {
  experiment_list <- experiments(.data)
  for (i in seq_along(quosure_list)) {
    experiment_list[[names(quosure_list)[i]]] <- rlang::eval_tidy(quosure_list[[i]],
                                                                data = experiment_list %>% 
                                                                  as.list())
  }
  # Provides the transmute functionality
  if (drop_unused) {
    experiment_list <- experiment_list[names(quosure_list)]
  }
  MultiAssayExperiment(experiments = experiment_list)
}

dummy.SummarizedExperiment <- function(.data, quosure_list, drop_unused = FALSE) {
  for (i in seq_along(quosure_list)) {
    assays(.data)[[names(quosure_list)[i]]] <- rlang::eval_tidy(quosure_list[[i]],
                                                                data = assays(.data) %>% 
                                                                  as.list())
  }
  # Provides the transmute functionality
  if (drop_unused) {
    assays(.data) <- assays(.data)[names(quosure_list)]
  }
  .data
}


