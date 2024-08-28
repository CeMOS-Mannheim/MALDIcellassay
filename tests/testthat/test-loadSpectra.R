test_that("Spectra loading works bruker format", {
  dataDir <- system.file("extdata", package="MALDIcellassay")
  unzip(file.path(dataDir, "example-raw-spectra.zip"))
  s <- loadSpectra("example-raw-spectra/", verbose = TRUE)
  unlink("example-raw-spectra/", recursive = TRUE)
  
  nm <- as.numeric(names(s))
  
  expect_true(length(nm) == 1)
  expect_true(all(is.numeric(nm)))
  expect_true(MALDIquant::isMassSpectrumList(s))
  expect_true(length(MALDIquant::intensity(s[[1]])) > 1)
})

test_that("Spectra loading works mzML format", {
  dataDir <- system.file("extdata", package="MALDIcellassay")
  s <- loadSpectraMzML(file.path(dataDir, "Koch2024mzML"), verbose = TRUE)
  
  nm <- as.numeric(names(s))
  
  expect_true(length(nm) == 1)
  expect_true(all(is.numeric(nm)))
  expect_true(MALDIquant::isMassSpectrumList(s))
  expect_true(length(MALDIquant::intensity(s[[1]])) > 1)
})