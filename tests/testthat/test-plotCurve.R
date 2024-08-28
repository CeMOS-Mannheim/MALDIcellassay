test_that("plotCurve is ggplot", {
  data("Blank2022res")
  
  # works with sd bars
  p <- plotCurves(Blank2022res, mzIdx = 1, errorbars = "sd")
  expect_true("ggplot" %in% class(p))
  expect_equal(length(p$layers), 3)
  expect_true("GeomLine" %in% class(p$layers[[1]]$geom))
  expect_true("GeomPoint" %in% class(p$layers[[2]]$geom))
  expect_true("GeomErrorbar" %in% class(p$layers[[3]]$geom))
  
  # works with sem bars
  p <- plotCurves(Blank2022res, mzIdx = 1, errorbars = "sem")
  expect_true("ggplot" %in% class(p))
  expect_equal(length(p$layers), 3)
  expect_true("GeomLine" %in% class(p$layers[[1]]$geom))
  expect_true("GeomPoint" %in% class(p$layers[[2]]$geom))
  expect_true("GeomErrorbar" %in% class(p$layers[[3]]$geom))
  
  # works with no bars
  p <- plotCurves(Blank2022res, mzIdx = 1, errorbars = "none")
  expect_true("ggplot" %in% class(p))
  expect_equal(length(p$layers), 2)
  expect_true("GeomLine" %in% class(p$layers[[1]]$geom))
  expect_true("GeomPoint" %in% class(p$layers[[2]]$geom))
  
  # works with multiple mzIdx
  l <- plotCurves(Blank2022res, mzIdx = 1:2, errorbars = "sd")
  expect_true("ggplot" %in% class(l[[1]]))
  expect_true("ggplot" %in% class(l[[2]]))
})

test_that("plotCurves stops as intended", {
  data("Blank2022res")
  expect_error(plotCurves(1))
  expect_error(plotCurves(Blank2022res, mzIdx = ""))
  expect_error(plotCurves(Blank2022res, mzIdx = 9999))
  expect_error(plotCurves(Blank2022res, mzIdx = -1))
})

test_that("plotCurves plots all curves if no mzIdx is given", {
  data("Blank2022res")
  expect_equal(length(plotCurves(Blank2022res, mzIdx = NULL)), 23)
})
