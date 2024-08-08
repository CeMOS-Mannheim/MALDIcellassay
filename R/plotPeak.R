#'  Plot a peak of interest from a MALDIassay object
#'
#' @param object object of class MALDIassay
#' @param mzIdx  numeric, index of mass of interest (see \code{getPeakStatistics()})
#' @param tol    numeric, tolerance around peak to plot
#'
#' @return
#' ggplot object
#'
#' @importFrom ggplot2 ggplot aes geom_line scale_y_continuous scale_color_viridis_d labs
#' @importFrom dplyr filter between bind_rows
#' @importFrom tibble tibble
#' @importFrom MALDIquant mass intensity
#' @importFrom forcats fct_reorder
#' @export
plotPeak <- function(object, mzIdx, tol = 0.8) {
  if (missing(mzIdx)) {
    stop("No mzIdx supplied.\n")
  }
  mz <- as.numeric(names(getCurveFits(object))[mzIdx])
  spec <- getAvgSpectra(object)

  df <- map(spec, function(x) {
    tibble(
      mass = mass(x),
      intensity = intensity(x)
    )
  }) %>%
    bind_rows(.id = "conc") %>%
    filter(between(mass, mz - tol, mz + tol)) %>%
    mutate(conc = as.numeric(conc)) %>%
    mutate(conc = fct_reorder(scales::scientific(conc), conc))

  title <- paste0("Profile of m/z ", round(mz, 2), " \u00B1 ", round(tol, 2), "Da")

  p <-ggplot() +
    geom_rect(aes(xmin = mz - getBinTol(object) * mz,
                  xmax = mz + getBinTol(object) * mz,
                  ymin = 0,
                  ymax = max(pull(df, intensity))*1.05),
              alpha=0.2,
              fill="black") +
    geom_line(data = df, aes(x = mass, y = intensity, col = conc)) +
    scale_y_continuous(limits = c(0, NA), expand = c(0,0)) +
    scale_color_viridis_d(end = 0.75, option = "C") +
    labs(
      x = "m/z",
      y = "Intensity",
      col = "Conc.",
      title = title
    )

  return(p)
}

#' Summary plot of a specific m/z with spectrum of peak and dose-response curve
#'
#' @param object    object of class MALDIassay.
#' @param mzIdx     numeric, index of mass of interest (see \code{getPeakStatistics()}).
#' @param tol       numeric, tolerance around peak to plot.
#' @param markValue numeric, value to mark in curve plot with a vertical line.
#' @param errorbars logical, add errorbars to plot representing standard deviation in regards to the measurment replicates.
#'
#' @return
#' ggplot object
#' @export
#'
#' @importFrom ggpubr ggarrange
plotPeakSummary <- function(object, mzIdx, tol = 4, markValue = NA, errorbars = TRUE) {
  if (missing(mzIdx)) {
    stop("No mzIdx supplied.\n")
  }
  peakProfile <- plotPeak(
    object = object,
    mzIdx = mzIdx,
    tol = tol
  )
  curve <- plotCurves(
    object = object,
    mzIdx = mzIdx,
    markValue = markValue,
    errorbars = errorbars
  )
  p <- ggarrange(curve, peakProfile + labs(title = NULL))
  return(p)
}
