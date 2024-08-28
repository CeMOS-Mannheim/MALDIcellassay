test_that("plotPeak is ggplot", {
  data("Blank2022res")
  
  p <- plotPeak(Blank2022res, mzIdx = 1, tol = 0.1)
  expect_true("ggplot" %in% class(p))
  expect_equal(length(p$layers), 2)
  expect_true("GeomRect" %in% class(p$layers[[1]]$geom))
  expect_true("GeomLine" %in% class(p$layers[[2]]$geom))
})

test_that("plotPeak stops as intended", {
  data("Blank2022res")
  expect_error(plotPeak(object = 1, mzIdx = 1))
  expect_error(plotPeak(Blank2022res))
  expect_error(plotPeak(Blank2022res, mzIdx = NULL))
  expect_error(plotPeak(Blank2022res, mzIdx = NA))
  expect_error(plotPeak(Blank2022res, mzIdx = 1, tol = 0))
  expect_error(plotPeak(Blank2022res, mzIdx = 1, tol = -1))
})


