test_that(".detectPeaks passes through peaks", {
  data("Blank2022peaks")
  peaks <- suppressMessages(.detectPeaks(Blank2022peaks, SNR = 3))
  
  expect_message(.detectPeaks(Blank2022peaks))
  expect_true(MALDIquant::isMassPeaksList(peaks))
  expect_error(.detectPeaks(1))
})

test_that(".detectPeaks turns spectra into peaks", {
  data("Blank2022spec")
  peaks <- .detectPeaks(Blank2022spec, SNR = 3)
  nm <- as.numeric(names(peaks))
  
  expect_true(MALDIquant::isMassPeaksList(peaks))
  expect_true(all(!is.na(nm)))
  expect_true(all(is.numeric(nm)))
}) 
