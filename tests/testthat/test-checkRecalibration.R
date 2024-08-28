test_that("checkRecalibration produces ggplot object", {
  data("Blank2022res")
  p <- checkRecalibration(Blank2022res)
  expect_true("ggplot" %in% class(p))
  expect_equal(length(p$layers), 5)
  expect_true("GeomLine" %in% class(p$layers[[1]]$geom))
  expect_true("GeomLinerange" %in% class(p$layers[[2]]$geom))
  expect_true("GeomVline" %in% class(p$layers[[3]]$geom))
  expect_true("GeomVline" %in% class(p$layers[[4]]$geom))
  expect_true("GeomVline" %in% class(p$layers[[5]]$geom))
})

test_that("checkRecalibration stops as intended", {
  expect_error(checkRecalibration(1))
})
