#' Check the recalibration of spectra from a MALDIassay object
#'
#' Dashed gray lines indicate the mz used for recalibration ± the tolerance.
#' Red dashed line indicate the mz used for recalibration and solid dashed lines indicate peaks.
#' The spectrum will show the peak used for recalibration ± 10x the tolerance.
#'
#' @param object  Object of class MALDIassay
#' @param idx     Numeric, index of spectrum to plot
#'
#' @return
#' ggplot object
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_linerange geom_vline scale_x_continuous scale_y_continuous labs
#' @importFrom scales comma
#' @importFrom MALDIquant mass intensity
#' @importFrom tibble tibble
#' @importFrom purrr map
#' @importFrom scales scientific

checkRecalibration <- function(object, idx) {
  if (!class(object) == "MALDIassay") {
    stop("object needs to be of class MALDIassay.")
  }

  normMz <- getNormMz(object)
  tol <- getNormMzTol(object)
  lowerVal <- normMz - 10 * tol
  upperVal <- normMz + 10 * tol

  spec <- getAvgSpectra(object)[idx]
  peaks <- getAvgPeaks(object)[idx]

  names(peaks) <- names(spec)
  normMeth <- getNormMethod(object)

  if (normMeth == "mz") {
    y_lab <- paste("Intensity normalized to", normMz)
  } else {
    y_lab <- paste0("Intensity (", normMeth, "-normalized)")
  }

  df <- map(spec, function(x) {
    df <- tibble(
      mass = mass(x),
      intensity = intensity(x)
    )
  }) %>%
    bind_rows(.id = "conc") %>%
    mutate(conc = scientific(as.numeric(conc)))

  peakdf <- map(peaks, function(x) {
    df <- tibble(
      mass = mass(x),
      intensity = intensity(x)
    )
  }) %>%
    bind_rows(.id = "conc") %>%
    mutate(conc = scientific(as.numeric(conc)))

  p <-
    df %>%
    filter(between(mass, lowerVal, upperVal)) %>%
    ggplot(aes(x = mass, y = intensity, col = factor(conc))) +
    geom_line() +
    geom_linerange(
      data = peakdf %>%
        filter(between(mass, lowerVal, upperVal)),
      aes(x = mass, ymin = 0, ymax = intensity)
    ) +
    geom_vline(aes(xintercept = normMz - tol), alpha = 0.6, linetype = "dashed") +
    geom_vline(aes(xintercept = normMz + tol), alpha = 0.6, linetype = "dashed") +
    geom_vline(aes(xintercept = normMz), alpha = 0.6, linetype = "dashed", col = "red") +
    scale_color_viridis_d(end = 0.75, option = "C") +
    labs(
      subtitle = paste("Dashed lines:", normMz, "±", tol, "m/z"),
      col = "Conc.",
      x = "m/z",
      y = y_lab)

  return(p)
}
