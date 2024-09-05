test_that("MALDIassay class can be initialized", {
  # Create an instance of MALDIassay
  instance <- new("MALDIassay")
  
  # Check if the instance is of class MALDIassay
  expect_s4_class(instance, "MALDIassay")
})

# Test Slot Accessibility
test_that("Slots in MALDIassay are accessible and can be modified", {
  # Create an instance of MALDIassay
  instance <- new("MALDIassay")
  
  # Access and modify slots
  instance@avgSpectra <- list(1)
  instance@avgPeaks <- list(2)
  instance@singlePeaks <- list(3)
  instance@singleSpecSpots <- "spot1"
  instance@normFactors <- c(0.5, 1.0)
  instance@mzShifts <- c(1, 2, 3)
  instance@fits <- list(model1 = "fit1")
  instance@stats <- data.frame(A = 1:3)
  instance@included_specIdx <- c(1, 2)
  instance@settings <- list(option1 = TRUE)
  
  # Check modified values
  expect_equal(instance@avgSpectra, list(1))
  expect_equal(instance@avgPeaks, list(2))
  expect_equal(instance@singlePeaks, list(3))
  expect_equal(instance@singleSpecSpots, "spot1")
  expect_equal(instance@normFactors, c(0.5, 1.0))
  expect_equal(instance@mzShifts, c(1, 2, 3))
  expect_equal(instance@fits, list(model1 = "fit1"))
  expect_equal(instance@stats, data.frame(A = 1:3))
  expect_equal(instance@included_specIdx, c(1, 2))
  expect_equal(instance@settings, list(option1 = TRUE))
})

test_that("MALDIassy_class handles invalid inputs gracefully", {
  avgSpecInvalid <- list("invalidMassSpec")  # Invalid list, not of MassSpectrum class
  validMassSpec <- list(MALDIquant::createMassSpectrum(mass = 1:10, intensity = 1:10))
  
  avgPeaksInvalid <- list("invalidMassPeaks")  # Invalid list, not of MassPeaks class
  validMassPeaks <- list(MALDIquant::createMassPeaks(mass = 1:10, intensity = 1:10, snr = 1:10))
  
  singlePeaksInvalid <- list("invalidSinglePeaks")
  
  # Initialize other valid slots
  spots <- c("spot1", "spot2", "spot3")
  prc <- list(normFac = c(0.5, 1.0), mzShift = c(0.1, 0.2), idx = c(1, 2, 3))
  res_list <- list(model1 = "fit1", model2 = "fit2")
  stat_df <- data.frame(mz = c(1, 1), A = 1:2, B = 4:5)  # Invalid: m/z column must have unique values
  
  # Test invalid avgSpectra
  expect_error(new("MALDIassay",
                   avgSpectra = avgSpecInvalid, avgPeaks = validMassPeaks,
                   singlePeaks = validMassPeaks, singleSpecSpots = spots,
                   normFactors = prc$normFac, mzShifts = prc$mzShift,
                   fits = res_list, stats = stat_df, included_specIdx = prc$idx,
                   settings = list(normMeth = "method")),
               "avgSpectra must be a list of class MALDIquant::MassSpectrum objects.")
  
  # Test invalid avgPeaks
  expect_error(new("MALDIassay",
                   avgSpectra = validMassSpec, avgPeaks = avgPeaksInvalid,
                   singlePeaks = validMassPeaks, singleSpecSpots = spots,
                   normFactors = prc$normFac, mzShifts = prc$mzShift,
                   fits = res_list, stats = stat_df, included_specIdx = prc$idx,
                   settings = list(normMeth = "method")),
               "avgPeaks must be a list of class MALDIquant::MassPeaks objects.")
  
  # Test invalid singlePeaks
  expect_error(new("MALDIassay",
                   avgSpectra = validMassSpec, avgPeaks = validMassPeaks,
                   singlePeaks = singlePeaksInvalid, singleSpecSpots = spots,
                   normFactors = prc$normFac, mzShifts = prc$mzShift,
                   fits = res_list, stats = stat_df, included_specIdx = prc$idx,
                   settings = list(normMeth = "method")),
               "singlePeaks must be a list of class MALDIquant::MassPeaks objects.")
  
  # Test lengths mismatch avgSpectra and avgPeaks
  expect_error(new("MALDIassay",
                   avgSpectra = validMassSpec, avgPeaks = list(),
                   singlePeaks = validMassPeaks, singleSpecSpots = spots,
                   normFactors = prc$normFac, mzShifts = prc$mzShift,
                   fits = res_list, stats = stat_df, included_specIdx = prc$idx,
                   settings = list(normMeth = "method")),
               "Length of avgPeaks can't be 0.")
  
  # Test non-unique m/z in stats
  expect_error(new("MALDIassay",
                   avgSpectra = validMassSpec, avgPeaks = validMassPeaks,
                   singlePeaks = validMassPeaks, singleSpecSpots = spots,
                   normFactors = prc$normFac, mzShifts = prc$mzShift,
                   fits = res_list, stats = stat_df, included_specIdx = prc$idx,
                   settings = list(normMeth = "method")))
})