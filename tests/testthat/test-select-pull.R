test_that("Pull yoinks an experiment", {
  data("mae")
  expect_s4_class(pull(mae), "SummarizedExperiment")
})

test_that("Pull yoinks the correct experiment", {
  data("mae")
  expect_equal(pull(mae), mae[[length(experiments(mae))]])
  expect_equal(pull(mae, 1), mae[[1]])
})

test_that("Select preserves the class", {
  data("mae")
  expect_s4_class(select(mae, cyto), "MultiAssayExperiment")
  expect_lte(length(experiments(select(mae, cyto))),
             length(experiments(mae)))
})
