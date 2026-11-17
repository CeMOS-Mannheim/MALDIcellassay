test_that("plotPeak is ggplot", {
  data("Blank2022res")

  # use the feature with the highest average intensity so the window
  # definitely contains peaks
  intmat <- getIntensityMatrix(Blank2022res)
  mzIdx <- which.max(colMeans(intmat, na.rm = TRUE))

  p <- plotPeak(Blank2022res, mzIdx = mzIdx, tol = 1)
  expect_true("ggplot" %in% class(p))
  # binning rectangle + profile curves + peak markers
  expect_equal(length(p$layers), 3)
  expect_true("GeomRect" %in% class(p$layers[[1]]$geom))
  expect_true("GeomLine" %in% class(p$layers[[2]]$geom))
  expect_true("GeomPoint" %in% class(p$layers[[3]]$geom))
})

test_that("plotPeak can hide the peak markers", {
  data("Blank2022res")

  intmat <- getIntensityMatrix(Blank2022res)
  mzIdx <- which.max(colMeans(intmat, na.rm = TRUE))

  p <- plotPeak(Blank2022res, mzIdx = mzIdx, tol = 1, showPeaks = FALSE)
  expect_equal(length(p$layers), 2)
  expect_true("GeomLine" %in% class(p$layers[[2]]$geom))
})

test_that("plotPeak draws centroided data as sticks", {
  data("Blank2022res")

  intmat <- getIntensityMatrix(Blank2022res)
  mzIdx <- which.max(colMeans(intmat, na.rm = TRUE))

  p <- plotPeak(Blank2022res, mzIdx = mzIdx, tol = 1, centroided = TRUE)
  # binning rectangle + centroids as sticks + peak markers
  expect_equal(length(p$layers), 3)
  expect_true("GeomSegment" %in% class(p$layers[[2]]$geom))
  expect_true("GeomPoint" %in% class(p$layers[[3]]$geom))

  p_no <- plotPeak(Blank2022res, mzIdx = mzIdx, tol = 1,
                   centroided = TRUE, showPeaks = FALSE)
  expect_equal(length(p_no$layers), 2)
  expect_true("GeomSegment" %in% class(p_no$layers[[2]]$geom))
})

test_that("plotPeak stops as intended", {
  data("Blank2022res")
  expect_error(plotPeak(object = 1, mzIdx = 1))
  expect_error(plotPeak(Blank2022res))
  expect_error(plotPeak(Blank2022res, mzIdx = NULL))
  expect_error(plotPeak(Blank2022res, mzIdx = NA))
  expect_error(plotPeak(Blank2022res, mzIdx = 1, tol = 0))
  expect_error(plotPeak(Blank2022res, mzIdx = 1, tol = -1))
  expect_error(plotPeak(Blank2022res, mzIdx = 1, centroided = "yes"))
  expect_error(plotPeak(Blank2022res, mzIdx = 1, showPeaks = NA))
})
