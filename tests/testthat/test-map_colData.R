exp <- curatedTCGAData::curatedTCGAData(
  "BRCA", "GISTIC*", "2.0.1", dry.run = FALSE
)

test_that("Local colData initially empty", {
  expect_false(TidyMultitude:::has_local_colData(exp))
})

exp <- map_colData(exp)

test_that("Local colData filled by map_colData", {
  expect_true(TidyMultitude:::has_local_colData(exp))
})



