test_that("filterVariance works", {
  expect_equal(length(filterVariance(1:10, method = "mean", verbose = TRUE)), 5)
  expect_equal(length(filterVariance(c(4, 4, 1:10), method = "median", verbose = TRUE)), 6)
  expect_equal(length(filterVariance(1:10, method = "q25", verbose = TRUE)), 7)
  expect_equal(length(filterVariance(1:10, method = "q75", verbose = TRUE)), 3)
  expect_equal(length(filterVariance(1:10, method = "none", verbose = TRUE)), 10)
})
