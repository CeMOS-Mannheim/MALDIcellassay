#' Check the recalibration of spectra from a MALDIassay object
#'
#' Dashed gray lines indicate the mz used for re-calibration ± the tolerance.
#' Red dashed line indicate the mz used for re-calibration and solid lines indicate peaks.
#' The spectrum will show the peak used for re-calibration ± 10x the tolerance.
#'
#' Profile data is drawn as lines, centroided data as sticks; the detected peaks
#' that are used for the recalibration are marked with open circles.
#'
#' @param object     Object of class MALDIassay
#' @param idx        Numeric, index of spectrum to plot
#' @param centroided Logical, `TRUE` for centroided data: the spectra are then
#'   drawn as sticks instead of lines
#' @param showPeaks  Logical, mark the detected peaks with open circles
#'
#' @return
#' ggplot object
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_segment geom_vline scale_color_viridis_d labs
#' @importFrom scales comma scientific
#' @importFrom MALDIquant mass intensity
#' @importFrom tibble tibble
#' @importFrom purrr map
#' @importFrom methods is
#' @importFrom dplyr bind_rows filter between mutate .data
#'
#' @examples
#' # see example for `fitCurve()` to see how this data was generated
#' data(Blank2022res)
#' checkRecalibration(Blank2022res, idx = 1:8)
#'
#' # centroided data is drawn as sticks
#' checkRecalibration(Blank2022res, idx = 1:8, centroided = TRUE)
checkRecalibration <- function(object, idx, centroided = FALSE, showPeaks = TRUE) {
  if (!is(object, "MALDIassay")) {
    stop("object needs to be of class MALDIassay.")
  }

  if (!is.logical(centroided) || length(centroided) != 1L || is.na(centroided)) {
    stop("centroided needs to be TRUE or FALSE.")
  }

  if (!is.logical(showPeaks) || length(showPeaks) != 1L || is.na(showPeaks)) {
    stop("showPeaks needs to be TRUE or FALSE.")
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
    tibble(
      mass = mass(x),
      intensity = intensity(x)
    )
  }) %>%
    bind_rows(.id = "conc") %>%
    mutate(conc = scientific(as.numeric(.data$conc)))

  peakdf <- map(peaks, function(x) {
    tibble(
      mass = mass(x),
      intensity = intensity(x)
    )
  }) %>%
    bind_rows(.id = "conc") %>%
    mutate(conc = scientific(as.numeric(.data$conc)))

  in_window <- function(data) {
    filter(data, between(.data$mass, lowerVal, upperVal))
  }

  p <- ggplot()

  if (centroided) {
    # a centroided spectrum has no profile, the centroids are drawn as sticks
    p <- p +
      geom_segment(
        data = in_window(df),
        aes(x = .data$mass, xend = .data$mass,
            y = 0, yend = .data$intensity,
            col = factor(.data$conc)),
        linewidth = 0.4, alpha = 0.85
      )
  } else {
    p <- p +
      geom_line(
        data = in_window(df),
        aes(x = .data$mass, y = .data$intensity, col = factor(.data$conc))
      )
  }

  if (showPeaks) {
    p <- p +
      geom_point(
        data = in_window(peakdf),
        aes(x = .data$mass, y = .data$intensity, col = factor(.data$conc)),
        shape = 1, size = 3, stroke = 1, show.legend = FALSE
      )
  }

  p <- p +
    geom_vline(aes(xintercept = normMz - tol), alpha = 0.6, linetype = "dashed") +
    geom_vline(aes(xintercept = normMz + tol), alpha = 0.6, linetype = "dashed") +
    geom_vline(aes(xintercept = normMz), alpha = 0.6, linetype = "dashed", col = "red") +
    scale_color_viridis_d(end = 0.75, option = "C") +
    labs(
      subtitle = paste("Dashed lines:", normMz, "\u00B1", tol, "m/z"),
      col = "Conc.",
      x = "m/z",
      y = y_lab
    )

  return(p)
}
