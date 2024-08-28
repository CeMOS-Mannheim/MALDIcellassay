test_that("getSNR works", {
  data("Blank2022res")
  expect_equal(getSNR(Blank2022res), 3)
})
