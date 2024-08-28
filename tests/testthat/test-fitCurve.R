test_that("fitCurve works", {
  data("Blank2022spec")
  normTol = 0.1
  normMethod = "mz"
  normMz = 760.585
  
  res <- fitCurve(spec = Blank2022spec,
                  SinglePointRecal = TRUE, 
                  normMz = normMz, 
                  alignTol = 0.1, 
                  normTol = normTol,
                  varFilterMethod = "mean",
                  normMeth = normMethod,
                  verbose = TRUE) 
  expect_true(isMALDIassay(res))
  expect_true(MALDIquant::isMassSpectrumList(getAvgSpectra(res)))
  expect_true(MALDIquant::isMassPeaksList(getAvgPeaks(res)))
  expect_true(MALDIquant::isMassPeaksList(getSinglePeaks(res)))
  expect_equal(length(getSpots(res)), length(Blank2022spec))
  expect_true(all(is.character(getSpots(res))))
  expect_equal(length(getCurveFits(res)), 23)
  
  stats <- getPeakStatistics(res, summarise = TRUE)
  expect_equal(colnames(stats), c("mz", "mzIdx", "pEC50", "R2", "log2FC", "SSMD", "V'", "Z'", "CRS"))
  expect_equal(dim(stats)[1], 23)
  expect_equal(length(res@included_specIdx), 32)
  expect_equal(getNormMethod(res), normMethod)
  expect_equal(getNormMzTol(res), normTol)
  expect_equal(length(getAppliedMzShift(res)), 32)
  expect_equal(length(getAppliedNormFactors(res)), 32)
  
  # test monoisotopic peaks filtering
  res <- fitCurve(spec = Blank2022spec,
                  SinglePointRecal = TRUE, 
                  normMz = normMz, 
                  alignTol = 0.1, 
                  normTol = normTol,
                  varFilterMethod = "mean",
                  normMeth = "median",
                  verbose = TRUE, 
                  monoisotopicFilter = TRUE) 
  expect_true(isMALDIassay(res))
  
  # test no normlization and no re-cal
  res <- fitCurve(spec = Blank2022spec,
                  SinglePointRecal = FALSE, 
                  normMz = normMz, 
                  alignTol = 0.1, 
                  normTol = normTol,
                  varFilterMethod = "mean",
                  normMeth = "none",
                  verbose = TRUE, 
                  monoisotopicFilter = TRUE) 
  expect_true(isMALDIassay(res))
})

test_that("fitCurve stops as intended", {
  data("Blank2022spec")
  
  expect_error(fitCurve(spec = Blank2022spec,
                        SinglePointRecal = TRUE, 
                        normMz = NULL, 
                        alignTol = 0.1, 
                        normTol = 0.1,
                        varFilterMethod = "mean",
                        normMeth = "mz",
                        verbose = TRUE))
  # no names given -> error
  spec <- Blank2022spec
  names(spec) <- NULL
  expect_error(fitCurve(spec = spec,
                        SinglePointRecal = TRUE, 
                        normMz = 760.585, 
                        alignTol = 0.1, 
                        normTol = 0.1,
                        varFilterMethod = "mean",
                        normMeth = "mz",
                        verbose = TRUE))
  
  # characters as names -> error
  names(spec) <- rep("test", 32)
  expect_error(
    suppressWarnings(
      fitCurve(spec = spec,
               SinglePointRecal = TRUE, 
               normMz = 760.585, 
               alignTol = 0.1, 
               normTol = 0.1,
               varFilterMethod = "mean",
               normMeth = "mz",
               verbose = TRUE)
    )
  )
  
  # names to short -> error
  names(spec) <- NULL
  names(spec) <- 1:31
  expect_error(fitCurve(spec = spec,
                        SinglePointRecal = TRUE, 
                        normMz = 760.585, 
                        alignTol = 0.1, 
                        normTol = 0.1,
                        varFilterMethod = "mean",
                        normMeth = "mz",
                        verbose = TRUE))
})
