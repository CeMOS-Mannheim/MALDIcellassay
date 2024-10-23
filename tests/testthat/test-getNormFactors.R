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

test_that("getNormFactors stops if allowNoMatch = FALSE and normMz was not found in one spectrum", {
  peaks <- peakListConstructor()
  MALDIquant::mass(peaks[[1]])[10] <- 12
  
  expect_error(getNormFactors(peaks = peaks, 
                              targetMz = 10, 
                              tol = 100, 
                              tolppm = TRUE, 
                              allowNoMatch = FALSE))
})

test_that("getNormFactors warns if allowNoMatch = FALSE and normMz was not found in one spectrum", {
  peaks <- peakListConstructor()
  MALDIquant::mass(peaks[[1]])[10] <- 12
  
  expect_warning(getNormFactors(peaks = peaks, 
                                targetMz = 10, 
                                tol = 100, 
                                tolppm = TRUE, 
                                allowNoMatch = TRUE))
})

test_that("getNormFactors stops if normMz was not found in any spectrum", {
  peaks <- peakListConstructor()
  
  expect_error(suppressWarnings(
    getNormFactors(peaks = peaks, 
                   targetMz = 12, 
                   tol = 100, 
                   tolppm = TRUE, 
                   allowNoMatch = TRUE)
  )
  )
  
})


test_that("getNormFactors works", {
  peaks <- peakListConstructor()
  n <- getNormFactors(peaks = peaks, targetMz = 10, tol = 100, tolppm = TRUE, allowNoMatch = TRUE)
  expect_type(n, "list")
  expect_equal(names(n), c("norm_factor", "specIdx"))
  expect_equal(length(n), 2)
  expect_equal(length(n$specIdx), 3)
  expect_equal(length(n$norm_factor), 3)
  expect_equal(n$norm_factor, c(10, 10, 10))
})
