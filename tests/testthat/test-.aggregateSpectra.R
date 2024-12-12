
test_that(".aggregateSpectra stops if spectra not named.", {
  data("Blank2022spec")
  
  spec <- Blank2022spec
  names(spec) <- NULL
  
  expect_error(.aggregateSpectra(spec = spec,
                                 averageMethod = "mean", 
                                 SNR = 3, 
                                 monoisotopicFilter = TRUE,
                                 binTol = 0.1,
                                 normMz = 760.585,
                                 normTol = 0.1,
                                 halfWindowSize = 20,
                                 verbose = FALSE))
})

test_that(".aggregateSpectra stops if spec is not list of spectra.", {
  data("Blank2022spec")
  
  # no spectra
  expect_error(.aggregateSpectra(spec = list(1, 2, 3),
                                 averageMethod = "mean", 
                                 SNR = 3, 
                                 monoisotopicFilter = TRUE,
                                 binTol = 0.1,
                                 normMz = 760.585,
                                 normTol = 0.1,
                                 halfWindowSize = 20,
                                 verbose = FALSE))
  
  # not a list of spectra
  expect_error(.aggregateSpectra(spec = Blank2022spec[[1]],
                                 averageMethod = "mean", 
                                 SNR = 3, 
                                 monoisotopicFilter = TRUE,
                                 binTol = 0.1,
                                 normMz = 760.585,
                                 normTol = 0.1,
                                 halfWindowSize = 20,
                                 verbose = FALSE))
})

test_that(".aggregateSpectra works with spectra.", {
  data("Blank2022spec")
  
  # spectra as input
  s <- .aggregateSpectra(spec = Blank2022spec,
                         averageMethod = "mean", 
                         SNR = 3, 
                         monoisotopicFilter = TRUE,
                         binTol = 0.1,
                         normMz = 761.585,
                         normTol = 0.1,
                         halfWindowSize = 20)
  
  expect_type(s, "list")
  expect_equal(names(s), c("intmat", "avgPeaksBinned", "avgSpec"))
  expect_true(is.matrix(s$intmat))
  expect_true(MALDIquant::isMassSpectrumList(s$avgSpec))
  expect_true(MALDIquant::isMassPeaksList(s$avgPeaksBinned))
  expect_true(length(s$avgSpec) == 8)
  expect_true(length(s$avgPeaksBinned) == 8)
  expect_true(dim(s$intmat)[1] == 8)
  
  # check if filtered (not monoisotopic) normMz was re-added
  mz <- as.numeric(colnames(s$intmat))
  expect_true(!is.na(MALDIquant::match.closest(761.585, mz, tolerance = 0.1)))
  
})

test_that(".aggregateSpectra works with peaks as input", {
  data("Blank2022peaks")
  
  # spectra as input
  s <- .aggregateSpectra(spec = Blank2022peaks,
                         averageMethod = "mean", 
                         SNR = 3, 
                         monoisotopicFilter = TRUE,
                         binTol = 0.1,
                         normMz = 761.585,
                         normTol = 0.1,
                         halfWindowSize = 20)
  
  expect_type(s, "list")
  expect_equal(names(s), c("intmat", "avgPeaksBinned", "avgSpec"))
  expect_true(is.matrix(s$intmat))
  expect_true(MALDIquant::isMassPeaksList(s$avgSpec))
  expect_true(MALDIquant::isMassPeaksList(s$avgPeaksBinned))
  expect_true(length(s$avgSpec) == 8)
  expect_true(length(s$avgPeaksBinned) == 8)
  expect_true(dim(s$intmat)[1] == 8)
  
  # check if filtered (not monoisotopic) normMz was re-added
  mz <- as.numeric(colnames(s$intmat))
  expect_true(!is.na(MALDIquant::match.closest(761.585, mz, tolerance = 0.1)))
  
})
