# helper 
spectraListConstructor <- function(n = 3, metaData = list(test = "test")) {
  l <- lapply(1:n, 
              function(x) 
              {MALDIquant::createMassSpectrum(mass = 1:10, 
                                              intensity = abs(rnorm(10, 0, 1)), 
                                              metaData = metaData)
              })
  return(l)
} 

peakListConstructor <- function(n = 3, metaData = list(test = "test")) {
  l <- lapply(1:n, 
              function(x) 
              {MALDIquant::createMassPeaks(mass = 1:10, 
                                           intensity = abs(rnorm(10, 0, 1)), 
                                           snr = 1:10,
                                           metaData = metaData)
              })
  return(l)
} 

test_that("getMzShift correctly checks input", {
  dummyPeaks <- peakListConstructor()
  expect_error(getMzShift(spectraListConstructor())) # no peaks
  expect_error(getMzShift(dummyPeaks)) # not named
  names(dummyPeaks) <- LETTERS[1:3]
  expect_error(suppressWarnings(getMzShift(dummyPeaks))) # no numeric names
})

test_that("getMzShift stops if a single condition gets filtered completly", {
  dummyPeaks <- peakListConstructor()
  MALDIquant::mass(dummyPeaks[[1]]) <- c(1:4, 6, 7:11)
  names(dummyPeaks) <- 1:3
  
  expect_error(suppressWarnings(getMzShift(dummyPeaks, targetMz = 5.1, tol = 0.15, verbose = TRUE)))
})

test_that("getMzShift finds mass shift correctly", {
  data("Blank2022peaks")
  # successful shift
  dummyPeaks <- peakListConstructor()
  names(dummyPeaks) <- 1:3
  shift <- getMzShift(dummyPeaks, targetMz = 5.1, tol = 0.15, verbose = TRUE)
  expect_equal(names(shift), c("mzshift", "specIdx"))
  expect_equal(length(shift$mzshift), 3)
  expect_true(length(shift$mzshift) == length(shift$specIdx))
  expect_true(is.numeric(shift$mzshift))
  
  shift <- getMzShift(Blank2022peaks, targetMz = 760.585, tol = 0.15, verbose = FALSE)
  expect_equal(length(shift$mzshift), length(Blank2022peaks))
  
  shift <- getMzShift(dummyPeaks, targetMz = 5.001, tol = 1000, verbose = FALSE, tolppm = TRUE)
  expect_equal(length(shift$mzshift), 3)
})

test_that("getMzShift stops if no peak can be found", {
  dummyPeaks <- peakListConstructor()
  names(dummyPeaks) <- 1:3
  expect_error(suppressWarnings(getMzShift(dummyPeaks, targetMz = 5.5, tol = 0.15, verbose = FALSE)))
})
