test_that("getVarFilterMethod works", {
  data("Blank2022res")
  expect_equal(getVarFilterMethod(Blank2022res), "mean")
})
