test_that("calculateChauvenetCriterion works", {
  set.seed(42)
  #no outlier
  sample <- rnorm(n = 8, mean = 0, sd = 0.01)
  expect_true(all(!calculateChauvenetCriterion(sample)))
  
  # introduce outlier
  sample[1] <- 1
  expect_true(calculateChauvenetCriterion(sample)[1])
  expect_warning(calculateChauvenetCriterion(1:2))
})
