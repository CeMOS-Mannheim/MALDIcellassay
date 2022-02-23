#' Check the recalibration of spectra from a MALDIassay object
#'
#' Dashed grey lines indicate the mz used for recalibration ± the tolerance.
#' Red dashed line indicate the mz used for recalibration and solid dashed lines indicate peaks.
#'
#' @param object  Object of class MALDIassay
#' @param idx     Numeric, index of spectrum to plot
#'
#' @return
#' ggplot object
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_linerange geom_vline scale_x_continuous scale_y_continuous labs theme_minimal
#' @importFrom scales comma
#' @importFrom MALDIquant mass intensity
#' @importFrom tibble tibble

checkRecalibration <- function(object, idx) {
  
  if(!class(object) == "MALDIassay") {
    stop("object needs to be of class MALDIassay.")
  }
  
  conc <- unique(getConc(object))
  
  
  normMz <- getNormMz(object)
  tol <- getNormMzTol(object)
  spec <- getAvgSpectra(object)
  peaks <- getAvgPeaks(object)
  normMeth <- getNormMethod(object)
  
  if(normMeth == "mz") {
    y_lab <- paste("Intensity normalized to", normMz)
  } else {
    y_lab <- paste0("Intensity (", normMeth, "-normalized)")
  }
  
  df_l <- lapply(X = idx, FUN = function(i) { 
    df <- tibble(mass = mass(spec[[i]]),
                 intensity = intensity(spec[[i]]))
    
  })
  names(df_l) <- conc[idx]
  df <- bind_rows(df_l, .id = "idx")
  
  peakdf_l <- lapply(X = idx, FUN = function(i) {
    peakdf <- tibble(mass = mass(peaks[[i]]),
                     intensity = intensity(peaks[[i]]))
  })
  names(peakdf_l) <- conc[idx]
  peakdf <- bind_rows(peakdf_l, .id = "idx")
  
  ggplot(df, aes(x = mass, y = intensity, col = factor(idx))) +
    geom_line() +
    geom_linerange(data = peakdf, aes(x = mass, ymin = 0, ymax = intensity)) +
    geom_vline(aes(xintercept = normMz - tol), alpha = 0.6, linetype = "dashed") +
    geom_vline(aes(xintercept = normMz + tol), alpha = 0.6, linetype = "dashed") +
    geom_vline(aes(xintercept = normMz), alpha = 0.6, linetype = "dashed", col = "red") +
    scale_x_continuous(limits = c(normMz - 10 * tol, normMz + 10*tol)) +
    scale_y_continuous(limits = c(0, 1.5)) +
    labs(subtitle = paste("Dashed lines:", normMz, "±", tol, "m/z"),
         col = "Conc.",
         x = "m/z",
         y = y_lab) +
    theme_minimal() -> p
  
  if(length(idx) == 1) {
    p <- p + 
      labs(title = paste("Avg. Spectrum of concentration", scales::comma(conc[idx], 0.0000001)))
  }
  return(p)
  
}

#' Calculate remaining calibration error of a MALDIassay object
#'
#' @param object Object of class MALDIassay
#'
#' @return
#' A tibble containing statistics about remaining calibration error
#' @export
#'
#' @importFrom tibble tibble
getRecalibrationError <- function(object) {
  peaks <- peaks2df(getAvgPeaks(object))
  mzdev <- getMzShift(peaksdf = peaks, 
                      tol = getNormMzTol(object), 
                      targetMz = getNormMz(object), 
                      tolppm = FALSE, 
                      allowNoMatch = FALSE)
  
  res_df <- tibble(meanAbs = mean(abs(mzdev$mzshift)),
                   sdAbs = sd(abs(mzdev$mzshift)),
                   maxAbs = max(abs(mzdev$mzshift)),
                   mean = mean(mzdev$mzshift),
                   sd = sd(mzdev$mzshift))
  return(res_df)
}
