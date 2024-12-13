test_that("preprocessing works with recal and mz norm", {
  data("Blank2022spec")
  peaks_single <- .detectPeaks(Blank2022spec, SNR = 3)
  
  prc <-.preprocess(peaks_single = peaks_single,
                     spec = Blank2022spec,
                     SinglePointRecal = TRUE,
                     normMz = 760.585,
                     normTol = 0.1,
                     normMeth = "mz",
                     alignTol = 0,
                     allowNoMatches = TRUE, 
                     verbose = FALSE)
  
  expect_equal(names(prc), expected = c("spec", "singlePeaks", "idx", "mzShift", "normFac"))
  expect_true(MALDIquant::isMassSpectrumList(prc$spec))
  expect_true(MALDIquant::isMassPeaksList(prc$singlePeaks))
  expect_true(is.integer(prc$idx))
  expect_true(length(prc$idx)>0)
  expect_true(is.double(prc$mzShift))
  expect_true(length(prc$mzShift)>1)
  expect_true(is.double(prc$normFac))
  expect_named(prc$spec)
  expect_named(prc$singlePeaks)
})

test_that("preprocessing works with peaks as input, recal and mz norm", {
  data("Blank2022peaks")
  
  prc <-.preprocess(peaks_single = Blank2022peaks,
                    spec = Blank2022peaks,
                    SinglePointRecal = TRUE,
                    normMz = 760.585,
                    normTol = 0.1,
                    normMeth = "mz",
                    alignTol = 0,
                    allowNoMatches = TRUE, 
                    verbose = FALSE)
  
  expect_equal(names(prc), expected = c("spec", "singlePeaks", "idx", "mzShift", "normFac"))
  expect_true(MALDIquant::isMassPeaksList(prc$spec))
  expect_true(MALDIquant::isMassPeaksList(prc$singlePeaks))
  expect_true(is.integer(prc$idx))
  expect_true(length(prc$idx)>0)
  expect_true(is.double(prc$mzShift))
  expect_true(length(prc$mzShift)>1)
  expect_true(is.double(prc$normFac))
  expect_named(prc$spec)
  expect_named(prc$singlePeaks)
})

test_that("preprocessing works without recal and TIC norm", {
  data("Blank2022spec")
  peaks_single <- .detectPeaks(Blank2022spec, SNR = 3)
  prc <- .preprocess(peaks_single = peaks_single,
                    spec = Blank2022spec,
                    SinglePointRecal = FALSE,
                    normMz = NULL,
                    normTol = 0.1,
                    normMeth = "TIC",
                    alignTol = 0,
                    allowNoMatches = TRUE, 
                    verbose = FALSE)
  
  expect_equal(names(prc), expected = c("spec", "singlePeaks", "idx", "mzShift", "normFac"))
  expect_true(MALDIquant::isMassSpectrumList(prc$spec))
  expect_true(MALDIquant::isMassPeaksList(prc$singlePeaks))
  expect_true(is.integer(prc$idx))
  expect_true(length(prc$idx)>0)
  expect_true(is.double(prc$mzShift))
  expect_true(length(prc$mzShift)==1)
  expect_true(is.double(prc$normFac))
  expect_named(prc$spec)
  expect_named(prc$singlePeaks)
})

test_that("preprocessing works with alignment", {
  data("Blank2022spec")
  peaks_single <- .detectPeaks(Blank2022spec, SNR = 3)
  prc <- MALDIcellassay:::.preprocess(peaks_single = peaks_single,
                     spec = Blank2022spec,
                     SinglePointRecal = FALSE,
                     normMz = NULL,
                     normTol = 0.1,
                     normMeth = "TIC",
                     alignTol = 0.5,
                     allowNoMatches = TRUE, 
                     verbose = FALSE)
  
  expect_equal(names(prc), expected = c("spec", "singlePeaks", "idx", "mzShift", "normFac"))
  expect_true(MALDIquant::isMassSpectrumList(prc$spec))
  expect_true(MALDIquant::isMassPeaksList(prc$singlePeaks))
  expect_true(is.integer(prc$idx))
  expect_true(length(prc$idx)>0)
  expect_true(is.double(prc$mzShift))
  expect_true(length(prc$mzShift)==1)
  expect_true(is.double(prc$normFac))
  expect_named(prc$spec)
  expect_named(prc$singlePeaks)
})