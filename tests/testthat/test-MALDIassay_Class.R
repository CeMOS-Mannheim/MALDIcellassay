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
  instance@avgSpectra <- Blank2022res@avgSpectra
  instance@avgPeaks <- Blank2022res@avgPeaks
  instance@singlePeaks <- Blank2022res@singlePeaks
  instance@singleSpecSpots <- Blank2022res@singleSpecSpots
  instance@normFactors <- Blank2022res@normFactors
  instance@mzShifts <- Blank2022res@mzShifts
  instance@fits <- Blank2022res@fits
  instance@stats <- Blank2022res@stats
  instance@included_specIdx <- Blank2022res@included_specIdx
  instance@settings <- list(
    Conc = Blank2022res@settings$Conc,
    normMz = 760.585,
    normTol = 0.1,
    varFilterMethod = "none",
    monoisotopicFilter = TRUE,
    alignTol = 0,
    SNR = 3,
    normMeth = "mz",
    binTol = 0.1,
    SinglePointRecal = TRUE
  )
  
  # Test if instance is valid
  expect_true(validObject(instance))
  
  # Test if instance is valid more explicitly (for coverage only)
  expect_true(.validMALDIassay(instance))
  
  # Test show method explicitly
  expect_output(show(instance), "MALDIassay object")
})

test_that("Slots in MALDIassay are accessible and can be modified", {
  # Create an instance of MALDIassay
  instance <- new("MALDIassay")
  
  # Access and modify slots
  instance@avgSpectra <- Blank2022res@avgSpectra
  instance@avgPeaks <- list()
  instance@singlePeaks <- Blank2022res@singlePeaks
  instance@singleSpecSpots <- Blank2022res@singleSpecSpots
  instance@normFactors <- Blank2022res@normFactors
  instance@mzShifts <- Blank2022res@mzShifts
  instance@fits <- Blank2022res@fits
  instance@stats <- Blank2022res@stats
  instance@included_specIdx <- Blank2022res@included_specIdx
  instance@settings <- list(
    Conc = Blank2022res@settings$Conc,
    normMz = 760.585,
    normTol = 0.1,
    varFilterMethod = "none",
    monoisotopicFilter = TRUE,
    alignTol = 0,
    SNR = 3,
    normMeth = "mz",
    binTol = 0.1,
    SinglePointRecal = TRUE
  )
  
  expect_match(.validMALDIassay(instance), 
                "Length of avgPeaks can't be 0.")
  
  instance@avgSpectra <- list()
  instance@avgPeaks <- Blank2022res@avgPeaks
  
  expect_equal(.validMALDIassay(instance), 
                "Length of avgSpectra can't be 0.")
  
  instance@avgSpectra <- list(1:32)
  
  expect_equal(.validMALDIassay(instance), 
                "avgSpectra must be a list of class MALDIquant::MassSpectrum or MALDIquant::MassPeaks objects.")
  
  instance@avgSpectra <- Blank2022res@avgSpectra
  instance@avgPeaks <- list(1:32)
  expect_equal(.validMALDIassay(instance), 
                "avgPeaks must be a list of class MALDIquant::MassPeaks objects.")
  
  instance@avgPeaks <- Blank2022res@avgPeaks
  instance@singlePeaks <- list()
  expect_equal(.validMALDIassay(instance), 
                "Length of singlePeaks can't be 0.")
  
  instance@singlePeaks <- Blank2022res@avgSpectra
  expect_equal(.validMALDIassay(instance), 
                "singlePeaks must be a list of class MALDIquant::MassPeaks objects.")
  
  instance@singlePeaks <- Blank2022res@singlePeaks[-1]
  expect_match(.validMALDIassay(instance), 
              "singlePeaks \\(\\d+\\) and included_specIdx \\(\\d+\\) must have the same length\\.")
  
  instance@singlePeaks <- Blank2022res@singlePeaks
  instance@avgPeaks <- Blank2022res@avgPeaks[-1]
  expect_match(.validMALDIassay(instance), 
               "avgSpectra \\(\\d+\\) and avgPeaks \\(\\d+\\) must have the same length\\.")
  
  instance@avgPeaks <- Blank2022res@avgPeaks
  instance@mzShifts <- Blank2022res@mzShifts[-1]
  expect_match(.validMALDIassay(instance), 
               "singlePeaks \\(\\d+\\) and mzShifts \\(\\d+\\) must have the same length or mzShifts must have length 1\\.")
  
  instance@mzShifts <- Blank2022res@mzShifts
  instance@normFactors <- Blank2022res@normFactors[-1]
  expect_match(.validMALDIassay(instance), 
               "singlePeaks \\(\\d+\\) and normFactors \\(\\d+\\) must have the same length or normFactors must have length 1\\.")
  
  instance@normFactors <- Blank2022res@normFactors
  instance@singleSpecSpots <- Blank2022res@singleSpecSpots[-1]
  expect_match(.validMALDIassay(instance), 
               "singlePeaks \\(\\d+\\) and singleSpecSpots \\(\\d+\\) must have the same length or singleSpecSpots must have length 1\\.")
  
  instance@singleSpecSpots <- Blank2022res@singleSpecSpots
  instance@fits <- Blank2022res@fits[-1]
  expect_match(.validMALDIassay(instance), 
               "fits \\(\\d+\\) and the unique m/z values in stats \\(\\d+\\) must have the same length\\.")
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