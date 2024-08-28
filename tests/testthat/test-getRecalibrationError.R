test_that("getRecalibrationError works", {
  data("Blank2022res")
  expect_true(tibble::is_tibble(getRecalibrationError(Blank2022res)))
  expect_true(is.numeric(getRecalibrationError(Blank2022res)[[1]]))
})
