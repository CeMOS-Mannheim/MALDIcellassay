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

test_that("normalize stops like intended", {
  data("Blank2022peaks")
  data("Blank2022spec")
  
  # spec not named
  spec <- Blank2022spec
  names(spec) <- NULL
  expect_error(normalize(spec = spec, peaks = Blank2022peaks, normMeth = "mz", normMz = 760.585, normTol = 0.1))
  
  # spec not named with numerics/concentrations
  names(spec) <- rep("a", length(spec))
  expect_error(
    suppressWarnings(
      normalize(spec = spec, peaks = Blank2022peaks, normMeth = "mz", normMz = 760.585, normTol = 0.1)
    )
  )
  
  # normMz not in data
  expect_error(
    suppressWarnings(
      normalize(spec = Blank2022spec, peaks = Blank2022peaks, normMeth = "mz", normMz = 0, normTol = 0.1)
    )
  )
  
  # normMz absent in a specific concentration
  spec <- spectraListConstructor()
  peaks <- peakListConstructor()
  
  names(spec) <-1:3
  names(peaks) <- 1:3
  
  MALDIquant::mass(spec[[2]])[10] <- 11
  MALDIquant::mass(peaks[[2]])[10] <- 11
  
  expect_error(
    suppressWarnings(
      normalize(spec = spec, peaks = peaks, normMeth = "mz", normMz = 10, normTol = 0.1)
    )
  )
  
})

test_that("normalize works", {
  data("Blank2022peaks")
  data("Blank2022spec")
  
  n <- normalize(spec = Blank2022spec, peaks = Blank2022peaks, normMeth = "mz", normMz = 760.585, normTol = 0.1)
  expect_type(n, "list")
  expect_equal(length(n), 4)
  expect_equal(names(n), c("spec", "peaks", "factor", "idx"))
  expect_true(MALDIquant::isMassSpectrumList(n$spec))
  expect_true(MALDIquant::isMassPeaksList(n$peaks))
  expect_true(all(is.numeric(n$factor)))
  expect_true(all(is.integer(n$idx)))
  expect_true(length(Blank2022spec) == length(n$factor))
  
  n <- normalize(spec = Blank2022spec, peaks = Blank2022peaks, normMeth = "TIC", normMz = 760.585, normTol = 0.1)
  expect_type(n, "list")
  expect_equal(length(n), 4)
  expect_equal(names(n), c("spec", "peaks", "factor", "idx"))
  expect_true(MALDIquant::isMassSpectrumList(n$spec))
  expect_true(MALDIquant::isMassPeaksList(n$peaks))
  expect_true(all(is.numeric(n$factor)))
  expect_true(all(is.integer(n$idx)))
  expect_true(length(Blank2022spec) == length(n$factor))
  
  n <- normalize(spec = Blank2022spec, peaks = Blank2022peaks, normMeth = "median", normMz = 760.585, normTol = 0.1)
  expect_type(n, "list")
  expect_equal(length(n), 4)
  expect_equal(names(n), c("spec", "peaks", "factor", "idx"))
  expect_true(MALDIquant::isMassSpectrumList(n$spec))
  expect_true(MALDIquant::isMassPeaksList(n$peaks))
  expect_true(all(is.numeric(n$factor)))
  expect_true(all(is.integer(n$idx)))
  expect_true(length(Blank2022spec) == length(n$factor))
})
