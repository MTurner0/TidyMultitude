require(curatedTCGAData)

exp <- curatedTCGAData(
  "BRCA", "GISTIC*", "2.0.1", dry.run = FALSE
)

# Vector of colnames to intersect by later
# Has to be pulled from global colData right now
test_by <- exp %>% colData() %>% colnames()

test_that("Local colData initially empty", {
  expect_false(has_local_colData(exp))
})

exp <- map_colData(exp)

test_that("Local colData filled by map_colData", {
  expect_true(has_local_colData(exp))
})

after <- intersect_colData(exp, by = test_by)

test_that("intersect_colData retains experiments", {
  expect_equal(names(exp), names(after))
})
