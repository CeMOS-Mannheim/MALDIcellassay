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
#' @importFrom ggplot2 ggplot aes geom_line, geom_linerange geom_vline scale_x_continuous scale_y_continuous labs theme_minimal
#' @importFrom scales comma
#' @importFrom MALDIquant mass intensity
#' @importFrom tibble tibble

checkRecalibration <- function(object, idx) {
  
  if(!class(object) == "MALDIassay") {
    stop("object needs to be of class MALDIassay.")
  }
  
  conc <- unique(getConc(object))
  
  normMz <- object@settings$normMz
  tol <- object@settings$normTol
  df <- tibble(mass = mass(object@avgSpectra[[idx]]),
               intensity = intensity(object@avgSpectra[[idx]]))
  peakdf <- tibble(mass = mass(object@avgPeaks[[idx]]),
                   intensity = intensity(object@avgPeaks[[idx]]))
  ggplot(df, aes(x = mass, y = intensity)) +
    geom_line() +
    geom_linerange(data = peakdf, aes(x = mass, ymin = 0, ymax = intensity), col = "red") +
    geom_vline(aes(xintercept = normMz - tol), alpha = 0.6, linetype = "dashed") +
    geom_vline(aes(xintercept = normMz + tol), alpha = 0.6, linetype = "dashed") +
    geom_vline(aes(xintercept = normMz), alpha = 0.6, linetype = "dashed", col = "red") +
    scale_x_continuous(limits = c(normMz - 5 * tol, normMz + 5*tol)) +
    scale_y_continuous(limits = c(0, 1.5)) +
    labs(title = paste("Avg. Spectrum of concentration", scales::comma(conc[idx], 0.0000001)),
         subtitle = paste("Dashed lines:", normMz, "±", tol, "m/z")) +
    theme_minimal() -> p
  return(p)
  
}
