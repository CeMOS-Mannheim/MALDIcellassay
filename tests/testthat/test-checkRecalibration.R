test_that("checkRecalibration produces ggplot object", {
  res <- MALDIcellassay::Blank2022res
  p <- checkRecalibration(res)
  expect_true("ggplot" %in% class(p))
  expect_equal(length(p$layers), 5)
  expect_true("GeomLine" %in% class(p$layers[[1]]$geom))
  expect_true("GeomPoint" %in% class(p$layers[[2]]$geom))
  expect_true("GeomVline" %in% class(p$layers[[3]]$geom))
  expect_true("GeomVline" %in% class(p$layers[[4]]$geom))
  expect_true("GeomVline" %in% class(p$layers[[5]]$geom))
})

test_that("checkRecalibration draws centroided data as sticks", {
  res <- MALDIcellassay::Blank2022res
  p <- checkRecalibration(res, centroided = TRUE)
  expect_equal(length(p$layers), 5)
  expect_true("GeomSegment" %in% class(p$layers[[1]]$geom))
  expect_true("GeomPoint" %in% class(p$layers[[2]]$geom))
})

test_that("checkRecalibration can hide the peaks", {
  res <- MALDIcellassay::Blank2022res

  p <- checkRecalibration(res, showPeaks = FALSE)
  expect_equal(length(p$layers), 4)
  expect_true("GeomLine" %in% class(p$layers[[1]]$geom))
  expect_false(any(vapply(p$layers, function(l) "GeomPoint" %in% class(l$geom), logical(1))))

  p <- checkRecalibration(res, centroided = TRUE, showPeaks = FALSE)
  expect_equal(length(p$layers), 4)
  expect_true("GeomSegment" %in% class(p$layers[[1]]$geom))
})

test_that("checkRecalibration stops as intended", {
  res <- MALDIcellassay::Blank2022res
  expect_error(checkRecalibration(1))
  expect_error(checkRecalibration(res, centroided = "yes"))
  expect_error(checkRecalibration(res, centroided = NA))
  expect_error(checkRecalibration(res, showPeaks = "yes"))
  expect_error(checkRecalibration(res, showPeaks = NA))
})
