#' Fit dose-response curves
#'
#' @param spec                List of MALDIquant::MassSpectrum
#' @param dir                 Character, directory for output
#' @param varFilterMethod     Character, function applied for high variance filtering. One of the following options `mean` (default), `median`, `q25`, `q75` or `none` (no filtering).
#' @param conc                Numeric vector, concentration for each spectrum. Length has to be the same as length of spec
#' @param unit                Character, unit of concentration. Used to calculate the concentration in Moles so that pIC50 is correct.
#'                            Set to "M" if you dont want changes in your concentrations.
#' @param normMz              Numeric, mz used for normalization AND for single point recalibration.
#' @param normTol             Numeric, tolerance in Dalton to match normMz
#' @param alignTol            Numeric, tolerance for spectral alignment in Dalton.
#' @param binTol              Numeric, tolerance for binning of peaks.
#' @param SNR                 Numeric, signal to noise ratio for peak detection.
#' @param allowNoMatches      Logical, if normMz can not be found in a spectrum, proceed and exclude spectrum or stop
#' @param normMeth            Character, normalization method. Can either be "TIC", "PQM", "median" or "mz". If "mz" then the normMz is used. If none no normalization is done.
#' @param saveIntensityMatrix Logical, save the intensity matrix as xlsx to the disk.
#' @param SinglePointRecal    Logical, perform single point recalibration to normMz
#' @param fc_thresh           Numeric, threshold for fold change above which curves are plotted. The fold-chage is calculated as max/min for a given m/z.
#' @param markValue           Numeric, value to mark in the resulting plot. Set to NA if no value needs to be marked.
#' @param plot                Logical, should the curves be plotted and written to dir?
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
                     dir,
                     conc = NA,
                     unit = c("M", "mM", "µM", "nM", "pM", "fM"),
                     varFilterMethod = c("mean", "median", "q25", "q75", "none"),
                     monoisotopicFilter = FALSE,
                     normMz = NULL,
                     normTol = 0.1,
                     alignTol = 0.01,
                     binTol = 0.0002,
                     SNR = 3,
                     allowNoMatches = TRUE,
                     normMeth = c("mz", "TIC", "PQN", "median", "none"),
                     saveIntensityMatrix = TRUE,
                     SinglePointRecal = TRUE,
                     fc_thresh = 1,
                     markValue = NA,
                     plot = TRUE) {

  ##### match & evaluate arguments ####
  normMeth <- match.arg(normMeth)
  unit <- match.arg(unit)
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

  if (plot | saveIntensityMatrix) {
    if (missing(dir)) {
      stop("argument `dir` is missing with no default.\n")
    }
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

  #### recal ####

  peaks_single <- detectPeaks(spec, SNR = SNR, method = "SuperSmoother")

  if (SinglePointRecal) {
    # perform single point mass recalibration
    mzShift <- getMzShift(
      peaksdf = peaks2df(peaks_single),
      tol = normTol,
      targetMz = normMz,
      tolppm = FALSE,
      allowNoMatch = TRUE
    )
    cat("found mz", normMz, "in", length(mzShift$specIdx), "/",
        length(spec), "spectra\n")
    cat(MALDIcellassay:::timeNow(), "mzshift was", mean(mzShift$mzshift),
        "in mean and", max(abs(mzShift$mzshift)), " abs. max.\n")

    spec <- shiftMassAxis(spec[mzShift$specIdx],
                          mzShift$mzshift)
    peaks_single <- shiftMassAxis(peaks_single[mzShift$specIdx],
                                  mzShift$mzshift)
    included_idx_recal <- mzShift$specIdx

    if(length(unique(nm)) != length(unique(nm[included_idx_recal]))) {
      # stop if a single condition got filtered completely
      u_nm <- unique(nm)
      u_fil <- unique(nm[included_idx_recal])
      label_removed <- u_nm[which(!(u_nm %in% u_fil))]

      stop("Could not find ", normMz, " in all spectra with label ",
           paste0(label_removed, collapse = ", "),
           ".\n Consider increasing tol.\n")
    }
    # update name string
    nm_new <- nm[included_idx_recal]
    # store included indices to be stored in result object
    included_specIdx <- included_idx_recal

  } else {
    mzShift <- list("mzshift" = 0)
    included_idx_recal <- 1:length(spec)
  }

  #### normalization ####
  cat(MALDIcellassay:::timeNow(), "normalizing... \n")
  switch(normMeth,
         "TIC" = {
           spec <- calibrateIntensity(spec, method = "TIC")
           norm_fac <- list("norm_factor" = 0)
           included_idx_norm <- 1:length(spec)
         },
         "PQN" = {
           spec <- calibrateIntensity(spec, method = "PQN")
           norm_fac <- list("norm_factor" = 0)
           included_idx_norm <- 1:length(spec)
         },
         "median" = {
           spec <- calibrateIntensity(spec, method = "median")
           norm_fac <- list("norm_factor" = 0)
           included_idx_norm <- 1:length(spec)
         },
         "mz" = {
           norm_fac <- getNormFactors(
             peaksdf = peaks2df(peaks_single),
             targetMz = normMz,
             tol = normTol,
             allowNoMatch = TRUE,
             tolppm = TRUE
           )
           spec <- normalizeByFactor(spec[norm_fac$specIdx], norm_fac$norm_factor)
           included_idx_norm <- norm_fac$specIdx

           if(SinglePointRecal) {
             # if recal was done use updated name string
             u_nm <- unique(nm_new)
             u_fil <- unique(nm_new[included_idx_norm])
           } else {
             u_nm <- unique(nm)
             u_fil <- unique(nm[included_idx_norm])
             # if no recal was done store included indices to be written to result object
             included_specIdx <- included_idx_norm
           }

           if(length(u_nm) != length(u_fil)) {
             # stop if a single condition got filtered completely

             label_removed <- u_nm[which(!(u_nm %in% u_fil))]

             stop("Could not find ", normMz, " in all spectra with label ",
                  paste0(label_removed, collapse = ", "),
                  ".\n Consider increasing tol.\n")
           }
         },
         "none" = {
           norm_fac <- list("norm_factor" = 0)
           included_idx_norm <- 1:length(spec)
         }
  )

  current_names <- names(spec)

  #### alignment ####
  cat(MALDIcellassay:::timeNow(), "aligning spectra... \n")
  spec <- alignSpectra(spec,
                       warpingMethod = "linear",
                       tolerance = alignTol,
                       noiseMethod = "SuperSmoother",
                       SNR = SNR,
                       reference = detectPeaks(averageMassSpectra(spec),
                                               method = "SuperSmoother",
                                               SNR = SNR
                       ),
                       allowNoMatches = allowNoMatches,
                       emptyNoMatches = allowNoMatches
  )

  res_list <- vector("list", length = length(unique(current_names)))
  names(res_list) <- unique(current_names)

  #### average spectra ####
  cat(MALDIcellassay:::timeNow(), "calculating average spectra... \n")
  spots <- extractSpots(spec)
  avg_spec <- averageMassSpectra(spec, labels = current_names)
  cat(MALDIcellassay:::timeNow(),
      "building intensity matrix and applying variance filter... \n")
  peaks <- detectPeaks(avg_spec, method = "SuperSmoother", SNR = SNR)

  if(monoisotopicFilter) {
    cat(MALDIcellassay:::timeNow(),
        "Filtering monoisotopic peaks...\n")
    # set it to be less restrictive then default settings
    peaks <- monoisotopicPeaks(peaks,
                               size = 2L:10L,
                               minCor = 0.85,
                               tolerance = 1e-3)
  }

  peaksBinned <- binPeaks(peaks, tolerance = binTol)

  # perform variance filtering
  intmat <- intensityMatrix(peaksBinned, avg_spec)
  cat("      Found", dim(intmat)[2], "peaks in total.\n")
  rownames(intmat) <- names(avg_spec)
  vars <- apply(intmat, 2, var)
  idx <- filterVariance(vars, method = varFilterMethod)

  mzhits <- as.numeric(colnames(intmat))[idx]

  # single spectra data
  allmz <- as.numeric(colnames(intmat))

  singlePeaks <- extractIntensity(
    createMassPeaks(
      mass = allmz,
      intensity = rep(1, length(allmz))
    ),
    spec = spec,
    tol = 0.5
  )

  intmatSingle <- intensityMatrix(singlePeaks, spec)
  rownames(intmatSingle) <- names(spec)

  # fit curves

  cat(MALDIcellassay:::timeNow(), "fitting curves... \n")
  res_list <- calculateCurveFit(intmat = intmat, idx = idx)

  # peak statistics
  stat_df <- calculatePeakStatistics(res_list, intmatSingle)

  if (saveIntensityMatrix) {
    cat(MALDIcellassay:::timeNow(), "writing intensity matrix...", "\n")

    # average spectra
    write.csv(
      x = as_tibble(intmat, rownames = NA),
      file = file.path(
        dir,
        paste0(
          as.character(Sys.Date()),
          "_intensityMatrix_",
          normMeth,
          "norm_avg.csv"
        )
      )
    )

    write.csv(
      x = as_tibble(intmatSingle, rownames = NA),
      file = file.path(
        dir,
        paste0(
          as.character(Sys.Date()), "_",
          normMeth, "_", normMz, "_",
          "_intensityMatrix_",
          normMeth,
          "norm_singleSpec.csv"
        )
      )
    )
    write.csv(
      x = as.data.frame(stat_df),
      file = file.path(
        dir,
        paste0(
          as.character(Sys.Date()), "_",
          normMeth, "_", normMz, "_",
          "_intensityMatrix_",
          normMeth,
          "norm_mzStats.csv"
        )
      )
    )
  }

  cat(MALDIcellassay:::timeNow(), "Done!", "\n")
  res_class <- new("MALDIassay",
                   avgSpectra = avg_spec,
                   avgPeaks = peaksBinned,
                   singlePeaks = singlePeaks,
                   singleSpecSpots = spots,
                   normFactors = norm_fac$norm_factor,
                   mzShifts = mzShift$mzshift,
                   fits = res_list,
                   stats = stat_df,
                   included_specIdx = included_specIdx,
                   settings = list(
                     Conc = as.numeric(nm),
                     dir = dir,
                     normMz = normMz,
                     normTol = normTol,
                     varFilterMethod = varFilterMethod,
                     monoisotopicFilter = monoisotopicFilter,
                     alignTol = alignTol,
                     SNR = SNR,
                     normMeth = normMeth,
                     SinglePointRecal = SinglePointRecal
                   )
  )
  if (plot) {
    cat(MALDIcellassay:::timeNow(), "plotting...", "\n")
    MALDIcellassay:::savePlots(res_class, fc_thresh = fc_thresh)
    cat(MALDIcellassay:::timeNow(), "plotting done!", "\n")
  }

  return(res_class)
}
