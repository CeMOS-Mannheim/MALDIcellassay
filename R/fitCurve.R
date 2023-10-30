#' Fit dose-response curves
#'
#' @param spec                List of MALDIquant::MassSpectrum
#' @param varFilterMethod     Character, function applied for high variance filtering. One of the following options `mean` (default), `median`, `q25`, `q75` or `none` (no filtering).
#' @param conc                Numeric vector, concentration for each spectrum. Length has to be the same as length of spec
#' @param unit                Character, unit of concentration. Used to calculate the concentration in Moles so that pIC50 is correct.
#'                            Set to "M" if you dont want changes in your concentrations.
#' @param monoisotopicFilter  Logical, filter peaks and just use monoisotopic peaks for curve fit.
#' @param averageMethod       Character, aggregation method for average mass spectra ("mean" or "median")
#' @param normMz              Numeric, mz used for normalization AND for single point recalibration.
#' @param normTol             Numeric, tolerance in Dalton to match normMz
#' @param alignTol            Numeric, tolerance for spectral alignment in Dalton.
#' @param binTol              Numeric, tolerance for binning of peaks.
#' @param SNR                 Numeric, signal to noise ratio for peak detection.
#' @param allowNoMatches      Logical, if normMz can not be found in a spectrum, proceed and exclude spectrum or stop
#' @param normMeth            Character, normalization method. Can either be "TIC", "PQM", "median" or "mz". If "mz" then the normMz is used. If none no normalization is done.
#' @param SinglePointRecal    Logical, perform single point recalibration to normMz
#'
#' @return
#' Object of class `MALDIassay`.
#' The most important slot is `fits` which contains the IC50 curve fits.
#' @export
#'
#' @importFrom MALDIquant removeBaseline calibrateIntensity alignSpectra averageMassSpectra detectPeaks binPeaks intensityMatrix match.closest createMassPeaks
#' @importFrom nplr nplr convertToProp getXcurve getYcurve getFitValues getX getY getEstimates getGoodness
#' @importFrom dplyr summarise mutate group_by %>% arrange left_join rename bind_rows filter pull slice_head slice_tail
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot geom_line geom_point scale_x_continuous theme_bw theme element_text labs aes ggsave geom_vline

fitCurve <- function(spec,
                     conc = NA,
                     unit = c("M", "mM", "µM", "nM", "pM", "fM"),
                     varFilterMethod = c("mean", "median", "q25", "q75", "none"),
                     monoisotopicFilter = FALSE,
                     averageMethod = c("mean", "median", "sum"),
                     normMz = NULL,
                     normTol = 0.1,
                     alignTol = 0.01,
                     binTol = 0.0002,
                     SNR = 3,
                     allowNoMatches = TRUE,
                     normMeth = c("mz", "TIC", "PQN", "median", "none"),
                     SinglePointRecal = TRUE) {

  ##### match & evaluate arguments ####
  normMeth <- match.arg(normMeth)
  unit <- match.arg(unit)
  averageMethod <- match.arg(averageMethod)
  varFilterMethod <- match.arg(varFilterMethod)

  unitFactor <- switch (unit,
                        "M" = 1,
                        "mM" = 1e-3,
                        "µM" = 1e-6,
                        "nM" = 1e-9,
                        "pM" = 1e-12,
                        "fM" = 1e-15
  )

  if(normMeth == "mz" & is.null(normMz)) {
    stop("Normalization to m/z is not possible when no m/z was supplied.\n")
  }

  if (!any(is.na(conc))) {
    if(as.numeric(conc)) {
      stop("conc needs to be a numeric vector of concentrations!\n")
    }
    # if conc is given update spectra names with conc
    names(spec) <- as.numeric(conc) * unitFactor
    conc <- as.numeric(conc) * unitFactor
  } else {
    # if conc is not given assume that spectra names are concentrations
    if(length(names(spec)) < 1) {
      stop("No concentrations provided.
           Either name spectra with concentrations or use conc. var
           to set them.\n")
    }
    names(spec) <- as.numeric(names(spec)) * unitFactor
  }
  nm <- names(spec)

  # check spectra for problematic meta data and remove it of needed
  spec <- .repairMetaData(spec)

  # make sure that spectra are in ascending order in regards to concentration
  order <- order(as.numeric(nm))
  nm <- nm[order]
  spec <- spec[order]


  if(!length(nm) == length(spec)) {
    stop("No concentrations provided.
         Either name spectra with concentrations or use 'conc' argument.")
  }

  #### re-calibration ####
  peaks_single <- .detectPeaks(spec, SNR = SNR, method = "SuperSmoother")

  if (SinglePointRecal) {
    # perform single point mass recalibration
    mzShift <- getMzShift(
      peaks = peaks_single,
      tol = normTol,
      targetMz = normMz,
      tolppm = FALSE
    )

    spec <- shiftMassAxis(spec[mzShift$specIdx],
                          mzShift$mzshift)
    peaks_single <- shiftMassAxis(peaks_single[mzShift$specIdx],
                                  mzShift$mzshift)
    included_idx_recal <- mzShift$specIdx

  } else {
    mzShift <- list("mzshift" = 0)
    included_idx_recal <- 1:length(spec)
  }

  #### normalization ####
  cat(MALDIcellassay:::timeNow(), "normalizing... \n")
  norm <- normalize(spec = spec, peaks = peaks_single, normMeth = normMeth)
  spec <- norm$spec
  peaks_single <- norm$peaks
  norm_fac <- norm$factor
  included_specIdx <- norm$idx

  current_names <- names(spec)


  #### alignment ####
  cat(MALDIcellassay:::timeNow(), "aligning spectra... \n")
  wf <- determineWarpingFunctions(l = peaks_single,
                                  tolerance = alignTol,
                                  method = "linear",
                                  allowNoMatches = allowNoMatches)

  spec <- warpMassSpectra(spec,
                          w = wf,
                          emptyNoMatches = allowNoMatches)
  names(spec) <- current_names
  peaks_single <- warpMassPeaks(peaks_single,
                                w = wf,
                                emptyNoMatches = allowNoMatches)
  names(peaks_single) <- current_names


  res_list <- vector("list", length = length(unique(current_names)))
  names(res_list) <- unique(current_names)

  #### average spectra ####
  cat(MALDIcellassay:::timeNow(), "calculating", averageMethod, "spectra... \n")
  spots <- extractSpots(spec)

  avg <- .aggregateSpectra(spec,
                           averageMethod = averageMethod,
                           SNR = SNR,
                           monoisotopicFilter = monoisotopicFilter,
                           binTol = binTol)

  idx <- filterVariance(apply(avg$intmat, 2, var),
                        method = varFilterMethod)

  mzhits <- as.numeric(colnames(avg$intmat))[idx]

  # single spectra data
  allmz <- as.numeric(colnames(avg$intmat))

  singlePeaks <- extractIntensity(mz = allmz,
                                  peaks = peaks_single,
                                  spec = spec,
                                  tol = normTol)

  # fit curves
  cat(MALDIcellassay:::timeNow(), "fitting curves... \n")
  res_list <- calculateCurveFit(intmat = avg$intmat, idx = idx)

  # peak statistics
  stat_df <- calculatePeakStatistics(curveFits = res_list,
                                     singlePeaks = singlePeaks)

  cat(MALDIcellassay:::timeNow(), "Done!", "\n")
  res_class <- new("MALDIassay",
                   avgSpectra = avg$avgSpec,
                   avgPeaks = avg$avgPeaksBinned,
                   singlePeaks = singlePeaks,
                   singleSpecSpots = spots,
                   normFactors = norm_fac$norm_factor,
                   mzShifts = mzShift$mzshift,
                   fits = res_list,
                   stats = stat_df,
                   included_specIdx = included_specIdx,
                   settings = list(
                     Conc = as.numeric(nm),
                     normMz = normMz,
                     normTol = normTol,
                     varFilterMethod = varFilterMethod,
                     monoisotopicFilter = monoisotopicFilter,
                     alignTol = alignTol,
                     SNR = SNR,
                     normMeth = normMeth,
                     binTol = binTol,
                     SinglePointRecal = SinglePointRecal
                   )
  )

  return(res_class)
}
