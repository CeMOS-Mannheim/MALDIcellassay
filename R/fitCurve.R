#' Fit dose-response curves
#'
#' @param spec                List of MALDIquant::MassSpectrum
#' @param dir                 Character, directory for output
#' @param conc                Numeric vector, concentration for each spectrum. Length has to be the same as length of spec
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
#'
#' @return
#' @export
#'
#' @importFrom MALDIquant removeBaseline calibrateIntensity alignSpectra averageMassSpectra detectPeaks binPeaks intensityMatrix match.closest
#' @importFrom nplr nplr convertToProp getXcurve getYcurve getFitValues getX getY getEstimates getGoodness
#' @importFrom dplyr summarise mutate group_by %>% as_tibble arrange left_join rename bind_rows filter
#' @importFrom tibble tibble
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot geom_line geom_point scale_x_continuous theme_bw theme element_text labs aes ggsave geom_vline

fitCurve <- function(spec,
                     dir,
                     conc = NA,
                     normMz = 760.585,
                     normTol = 0.1,
                     alignTol = 0.05,
                     binTol = 0.25,
                     SNR = 3,
                     allowNoMatches =TRUE,
                     normMeth = c("mz", "TIC", "PQN", "median", "none"),
                     saveIntensityMatrix = TRUE,
                     SinglePointRecal = TRUE,
                     fc_thresh = 1,
                     markValue = NA) {
  normMeth <- match.arg(normMeth)

  if(!any(is.na(conc))) {
    names(spec) <- conc
  }
    nm <- names(spec)

    peaks_single <- detectPeaks(spec, SNR = SNR, method = "SuperSmoother")

  if(SinglePointRecal) {
    # perform single point mass recalibration
    mzShift <- getMzShift(peaksdf = peaks2df(peaks_single),
                          tol = normTol,
                          targetMz = normMz,
                          tolppm = FALSE,
                          allowNoMatch = TRUE)
    cat("found mz", normMz, "in", length(mzShift$specIdx), "/", length(spec), "spectra\n")
    cat(MALDIcellassay:::timeNow(), "mzshift was", mean(mzShift$mzshift), "in mean and", max(abs(mzShift$mzshift)), " abs. max.\n")
    spec <- shiftMassAxis(spec[mzShift$specIdx], mzShift$mzshift)
    peaks_single <- shiftMassAxis(peaks_single[mzShift$specIdx], mzShift$mzshift)
  }

  cat(MALDIcellassay:::timeNow(), "normalizing... \n")
  switch(normMeth,
         "TIC" = {
           spec <- calibrateIntensity(spec, method = "TIC")
         },
         "PQN" = {
           spec <- calibrateIntensity(spec, method = "PQN")
         },
         "median" = {
           spec <- calibrateIntensity(spec, method = "median")
         },
         "mz" = {
           norm_fac <- getNormFactors(peaksdf = peaks2df(peaks_single),
                                      targetMz = normMz,
                                      tol = normTol,
                                      allowNoMatch = TRUE,
                                      tolppm = TRUE)
           spec <- normalizeByFactor(spec[norm_fac$specIdx], norm_fac$norm_factor)
         },
         "none" = {
           # dont normalize
         }
  )

  current_names <- names(spec)

  cat(MALDIcellassay:::timeNow(), "aligning spectra... \n")
  spec <- alignSpectra(spec, warpingMethod = "linear",
                       tolerance = alignTol,
                       allowNoMatches = allowNoMatches,
                       emptyNoMatches = allowNoMatches)

  res_list <- vector("list", length = length(unique(current_names)))
  names(res_list) <- unique(current_names)


  cat(MALDIcellassay:::timeNow(), "calculating average spectra... \n")
  avg_spec <- averageMassSpectra(spec, labels = current_names)
  cat(MALDIcellassay:::timeNow(), "building intensity matrix and applying variance filter... \n")
  peaks <- detectPeaks(avg_spec, method = "SuperSmoother", SNR = SNR)
  peaksBinned <- binPeaks(peaks, tolerance = binTol)

  # perform variance filtering
  intmat <- intensityMatrix(peaksBinned, avg_spec)
  rownames(intmat) <- names(avg_spec)
  vars <- apply(intmat, 2, var)
  idx <- which(vars > mean(vars))
  cat("found", length(idx), "peaks with high variance.\n")

  mzhits <- as.numeric(colnames(intmat))[idx]


  cat(MALDIcellassay:::timeNow(), "fitting curves... \n")
  current_res <- vector("list", length = length(idx))
  names(current_res) <- colnames(intmat[,idx])
  for(j in 1:length(idx)) {
    df <- intmat[,idx[j]] %>%
      as_tibble() %>%
      mutate(conc = rownames(intmat)) %>%
      mutate(conc = as.numeric(conc)) %>%
      arrange(conc)
    concLog <- log10(df$conc)
    if(any(concLog == -Inf)) {
      concLog[which(concLog == -Inf)] <- (min(concLog[which(!concLog == -Inf)])-1)
    }
    df <- df %>%
      mutate(concLog = concLog)
    resp <- convertToProp(y = df$value)
    model <- nplr(x = concLog, y = resp, useLog = FALSE, npars = 4)

    current_res[[j]] <- list(model = model,
                             df = df)

    res_list <- current_res
  }

  cat(MALDIcellassay:::timeNow(), "plotting...", "\n")



  for(mz in as.numeric(names(res_list))) {
    model <- res_list[[as.character(mz)]]$model
    df <- res_list[[as.character(mz)]]$df

    ic50 <- 10^getEstimates(model, targets = 0.5)[,3]
    min <- min(df$value)
    max <- max(df$value)
    fc_window <- max/min

    if(fc_window >= fc_thresh) {
    df_C <- tibble(xC = getXcurve(model), yC = getYcurve(model))
    df_P <- tibble(x = getX(model), y = getY(model))
    getGoodness(model)[[1]] -> R2

    ggplot(data = df_P, aes(x = x, y = y)) +
      geom_line(data = df_C, aes(x = xC, y = yC)) +
      geom_point() +
      scale_x_continuous(labels = c(0, 10^df_P$x[-1]), breaks = df_P$x) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Conc.",
           y = "relative Int. [% of max Int.]",
           title = paste0("mz ", round(mz,2), " Da, R\u00B2=", round(R2,4), "\n",
                           "IC50=", round(ic50,3), " min=", round(min, 4), " max=", round(max, 4), " FC=", round(fc_window, 4))) -> p


    if(!is.na(markValue)) {
      p <- p + geom_vline(aes(xintercept = markValue), linetype = "dashed")
    }
    ggsave(filename = file.path(dir, paste0(as.character(Sys.Date()),"_plotR2_", normMeth, "norm_", round(mz,2),".png")), plot = p)
    }
  }
  cat(MALDIcellassay:::timeNow(), "plotting done!", "\n")
  if(saveIntensityMatrix) {
    cat(MALDIcellassay:::timeNow(), "writing intensity matrix...", "\n")

    # average spectra
    write.csv(x = as_tibble(intmat[,idx], rownames = NA),
               file = file.path(dir,
                                paste0(as.character(Sys.Date()),
                                       "_intensityMatrix_",
                                       normMeth,
                                       "norm_avg.csv")))
    allmz <- as.numeric(colnames(intmat))
    singlePeaks <- extractIntensity(createMassPeaks(mass = allmz,
                                                    intensity = rep(1, length(allmz))),
                                    spec = spec)

    # peak statistics
    fit_df <- lapply(res_list, function(x) {
      model <- x$model
      fc_window <- max(x$df$value)/min(x$df$value)
      res_df <- as_tibble(nplr::getGoodness(model)) %>%
        mutate(fc_window = fc_window)
      return(res_df)
    }) %>%
      bind_rows(.id = "mz") %>%
      rename("R2" = "gof")

    intmatSingle <- intensityMatrix(singlePeaks, spec)
    intmatSingle %>%
      as_tibble() %>%
      mutate(sample = names(spec)) %>%
      gather(mz, int, -sample) %>%
      group_by(sample, mz) %>%
      summarise(min = min(int, na.rm = TRUE),
                mean = mean(int, na.rm = TRUE),
                max = max(int, na.rm = TRUE),
                stdev = sd(int, na.rm = TRUE),
                "cv%" = stdev/mean*100) %>%
      left_join(fit_df, by = "mz") %>%
      filter(!is.na(R2))-> stat_df

    rownames(intmatSingle) <- names(spec)
    write.csv(x = as_tibble(intmatSingle, rownames = NA),
               file = file.path(dir,
                                paste0(as.character(Sys.Date()), "_",
                                       normMeth, "_", normMz, "_",
                                       "_intensityMatrix_",
                                       normMeth,
                                       "norm_singleSpec.csv")))
    write.csv(x = as.data.frame(stat_df),
               file = file.path(dir,
                                paste0(as.character(Sys.Date()), "_",
                                       normMeth, "_", normMz, "_",
                                       "_intensityMatrix_",
                                       normMeth,
                                       "norm_mzStats.csv")))
  }

  cat(MALDIcellassay:::timeNow(), "Done!", "\n")
}
