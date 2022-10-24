
#Testing mutate function
test_that("Subsetting data on Summarized Experiment", {
  data("mae")
  mae1=mutate(mae,subset = cyto[1:5, 1:3])
  expect_equal(dim(mae1[[3]]),dim(mae[[2]][1:5, 1:3]) )
  expect_s4_class(mae1[[3]], "SummarizedExperiment")
})

test_that("Adding a new Assay to Summarized Experiment", {
  data("mae")
  mae2=mae %>% mutate(cyto, plus_one = cyto_conc + 1)
  
  expect_equal(c(SummarizedExperiment::assayNames(mae2[[2]])), c("cyto_conc", "plus_one"))
})

#Testing transmute function

test_that("Subsetting data on Summarized Experiment while dropping rest of summarized experiments", {
  data("mae")
  mae3=transmute(mae,subset = cyto[1:5, 1:3])
  expect_equal(dim(mae3[[1]]),dim(mae[[2]][1:5, 1:3]) )
  expect_equal(mae3[['subset']], mae[[2]][1:5, 1:3])
  expect_s4_class(mae3[[1]], "SummarizedExperiment")
})

test_that("Adding new assay and drops existing assays from Summarized Experiment", {
  data("mae")
  mae4=mae %>% transmute(cyto, plus_one = cyto_conc + 1)
  expect_equal(c(SummarizedExperiment::assayNames(mae4[[2]])), c("plus_one"))
})