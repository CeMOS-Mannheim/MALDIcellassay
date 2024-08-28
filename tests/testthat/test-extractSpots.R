test_that("extractSpots extracts spots if present", {
  data("Blank2022spec")
  
  spots <- extractSpots(Blank2022spec)
  expect_true(all(is.character(spots)))
  expect_true(nchar(spots[1])>1) # spots consist of at least two characters (e.g. A1)
})

test_that("extractSpots extracts no spots if not present", {
  dataDir <- system.file("extdata", package="MALDIcellassay")
  s <- loadSpectraMzML(file.path(dataDir, "Koch2024mzML"), verbose = FALSE)
  
  # mzML data has no spots meta data
  spots <- extractSpots(s)
  expect_true(spots[1] == "")
})

test_that("extractSpots stops if no spectra/peak provided", {
  expect_error(extract_spots(1))
})

test_that("extractSpots extracts spots from peaks", {
  data("Blank2022peaks")
  
  spots <- extractSpots(Blank2022peaks)
  expect_true(all(is.character(spots)))
  expect_true(nchar(spots[1])>1) # spots consist of at least two characters (e.g. A1)
})