test_that("getAllMz works", {
  expect_equal(length(getAllMz(Blank2022res, excludeNormMz = FALSE)), 23)
  expect_true(all(is.numeric(getAllMz(Blank2022res, FALSE))))
  expect_equal(
    suppressWarnings(
      getAllMz(Blank2022res, excludeNormMz = TRUE)
    ), NA_integer_)
expect_warning(getAllMz(Blank2022res, excludeNormMz = TRUE))
})
