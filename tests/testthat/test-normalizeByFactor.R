# helper 
spectraListConstructor <- function(n = 3, metaData = list(test = "test")) {
  l <- lapply(1:n, 
              function(x) 
              {MALDIquant::createMassSpectrum(mass = 1:10, 
                                              intensity = 1:10, 
                                              metaData = metaData)
              })
  return(l)
} 

peakListConstructor <- function(n = 3, metaData = list(test = "test")) {
  l <- lapply(1:n, 
              function(x) 
              {MALDIquant::createMassPeaks(mass = 1:10, 
                                           intensity = 1:10, 
                                           snr = 1:10,
                                           metaData = metaData)
              })
  return(l)
} 

test_that("normalizeByFactor stops as intended", {
  expect_error(normalizeByFactor(spec = spectraListConstructor(), factors = 1:2))
})

test_that("normalizeByFactor works",  {
  normSpec <- normalizeByFactor(spec = spectraListConstructor(), 1:3)
  expect_equal(intensity(normSpec[[1]]), 1:10)
  expect_equal(intensity(normSpec[[2]]), 1:10/2)
  
  normPeaks <- normalizeByFactor(spec = peakListConstructor(), 1:3)
  expect_equal(intensity(normSpec[[1]]), 1:10)
  expect_equal(intensity(normSpec[[2]]), 1:10/2)
})