test_that("getSingleSpecIntensity works", {
  expect_equal(length(getSingleSpecIntensity(Blank2022res, 1)), 32)
})
