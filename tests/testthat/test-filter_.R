data("mae")

first_vis <- filter_colData(mae, Visit == 1)

test_that("colData is filtered", {
  checkpoint <- lapply(
    names(first_vis), {
      function(x) all(colData(first_vis[[x]])[["Visit"]] == 1)
    }
  )
  lapply(checkpoint, expect_true)
})
