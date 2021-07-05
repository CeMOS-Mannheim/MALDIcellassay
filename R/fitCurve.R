#' Fit dose-response curves
#'
#' @param spec                List of MALDIquant::MassSpectrum
#' @param dir                 Character, directory for output
#' @param conc                Numeric vector, concentration for each spectrum. Length has to be the same as length of spec
#' @param normMz              Numeric, mz used for normalization AND for single point recalibration.
#' @param normTol             Numeric, tolerance in Dalton to match normMz
#' @param alignTol            Numeric, tolerance for spectral alignment in Dalton.
#' @param binTol              Numeric, tolerance for binning of peaks. D
#' @param SNR                 Numeric, signal to noise ratio for peak detection.
#' @param allowNoMatches      Logical, if normMz can not be found in a spectrum, proceed and exclude spectrum or stop
#' @param normMeth            Character, normalization method. Can either be "TIC" or "mz". If "mz" then the normMz is used.
#' @param saveIntensityMatrix Logical, save the intensity matrix as xlsx to the disk.
#' @param SinglePointRecal    Logical, perform single point recalibration to normMz
#'
#' @return
#' @export
#'
#' @importFrom MALDIquant removeBaseline calibrateIntensity alignSpectra averageMassSpectra detectPeaks binPeaks intensityMatrix match.closest
#' @importFrom nplr nplr convertToProp getXcurve getYcurve getFitValues getX getY
#' @importFrom dplyr summarise mutate group_by %>% as_tibble arrange
#' @importFrom tibble tibble
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot geom_line geom_point scale_x_continuous theme_bw theme element_text labs aes ggsave

fitCurve <- function(spec,
                     dir,
                     conc = NA,
                     normMz = 760.585,
                     normTol = 0.1,
                     alignTol = 0.05,
                     binTol = 0.25,
                     SNR = 3,
                     allowNoMatches =TRUE,
                     normMeth = c("mz", "TIC"),
                     saveIntensityMatrix = TRUE,
                     SinglePointRecal = TRUE) {
  normMeth <- match.arg(normMeth)

  if(!any(is.na(conc))) {
    names(spec) <- conc
  }
    nm <- names(spec)

  if(SinglePointRecal) {
    # perform single point mass recalibration
    peaks <- detectPeaks(spec, SNR = SNR, method = "SuperSmoother")
    mzShift <- getMzShift(peaksdf = peaks2df(peaks),
                          tol = normTol,
                          targetMz = normMz,
                          tolppm = FALSE,
                          allowNoMatch = TRUE)
    cat("found mz", normMz, "in", length(mzShift$specIdx), "/", length(spec), "spectra\n")
    cat(MALDIcellassay:::timeNow(), "mzshift was", mean(mzShift$mzshift), "in mean and", max(abs(mzShift$mzshift)), " abs. max.\n")
    spec <- shiftMassAxis(spec[mzShift$specIdx], mzShift$mzshift)
  }

  cat(MALDIcellassay:::timeNow(), "normalizing... \n")
  switch(normMeth,
         "TIC" = {
           spec <- calibrateIntensity(spec, method = "TIC")
         },
         "mz" = {
           peaks <- detectPeaks(spec, SNR = SNR, method = "SuperSmoother")
           norm_fac <- getNormFactors(peaksdf = peaks2df(peaks),
                                      targetMz = normMz,
                                      tol = normTol,
                                      allowNoMatch = TRUE,
                                      tolppm = TRUE)
           spec <- normalizeByFactor(spec[norm_fac$specIdx], norm_fac$norm_factor)
         }
  )

  cat(MALDIcellassay:::timeNow(), "aligning spectra... \n")
  spec <- alignSpectra(spec, warpingMethod = "linear",
                       tolerance = alignTol,
                       allowNoMatches = allowNoMatches,
                       emptyNoMatches = allowNoMatches)

  res_list <- vector("list", length = length(unique(nm)))
  names(res_list) <- unique(nm)

  current_names <- nm

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
    df_C <- tibble(xC = getXcurve(model), yC = getYcurve(model))
    df_P <- tibble(x = getX(model), y = getY(model))
    df_P %>%
      mutate(yfit = getFitValues(model)) %>%
      mutate(Sres = (y-yfit)^2,
             Stot = (y-mean(y))^2) %>%
      summarise(SStot = sum(Stot),
                SSres =sum(Sres)) %>%
      mutate(R2 = 1-SSres/SStot) %>%
      pull(R2) -> R2

    ggplot(data = df_P, aes(x = x, y = y)) +
      geom_line(data = df_C, aes(x = xC, y = yC)) +
      geom_point() +
      scale_x_continuous(labels = c(0, 10^df_P$x[-1]), breaks = df_P$x) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      labs(x = "Conc.",
           y = "relative Int. [% of max Int.]",
           title = paste0( "mz ", round(mz,2), " Da, R\u00B2=", round(R2,4))) -> p
    ggsave(filename = file.path(dir, paste0(as.character(Sys.Date()),"_plotR2_", normMeth, "norm_", round(mz,2),".png")), plot = p)
  }
  cat(MALDIcellassay:::timeNow(), "plotting done!", "\n")
  if(saveIntensityMatrix) {
    cat(MALDIcellassay:::timeNow(), "writing intensity matrix...", "\n")
    write.csv(x = as_tibble(intmat[,idx], rownames = NA),
               file = file.path(dir,
                                paste0(as.character(Sys.Date()),
                                       "_intensityMatrix_",
                                       normMeth,
                                       "norm_avg.xlsx")))
    singlePeaks <- detectPeaks(spec, method = "SuperSmoother", SNR = SNR)
    singlePeaks <- binPeaks(singlePeaks, tolerance = binTol)
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
                "cv%" = stdev/mean*100) -> stat_df

    rownames(intmatSingle) <- names(spec)
    mzsingle <- as.numeric(colnames(intmatSingle))

    idx_single <- match.closest(mzhits, mzsingle)
    write.csv(x = as_tibble(intmatSingle[,idx_single], rownames = NA),
               file = file.path(dir,
                                paste0(as.character(Sys.Date()),
                                       "_intensityMatrix_",
                                       normMeth,
                                       "norm_singleSpec.xlsx")))
    write.csv(x = as.data.frame(stat_df),
               file = file.path(dir,
                                paste0(as.character(Sys.Date()),
                                       "_intensityMatrix_",
                                       normMeth,
                                       "norm_mzStats.xlsx")))
  }

  cat(MALDIcellassay:::timeNow(), "Done!", "\n")
}
