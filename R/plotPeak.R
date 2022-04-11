#'  Plot a peak of interest from a MALDIassay object
#'
#' @param object object of class MALDIassay
#' @param mzIdx  numeric, index of mass of interest (see \code{getPeakStatistics()})
#' @param tol    numeric, tolerance around peak to plot
#'
#' @return
#' ggplot object
#'
#' @importFrom ggplot2 ggplot aes geom_line scale_y_continuous scale_color_viridis_d theme_minimal labs
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
  conc <- unique(getConc(object))

  idx <- 1:length(spec)
  df_l <- lapply(X = idx, FUN = function(i) {
    df <- tibble(
      mass = mass(spec[[i]]),
      intensity = intensity(spec[[i]])
    )
  })
  names(df_l) <- conc[idx]
  df <- bind_rows(df_l, .id = "idx")

  title <- paste0("Profile of m/z ", round(mz, 2), " ± ", round(tol, 2), "Da")

  p <- df %>%
    filter(between(mass, mz - tol, mz + tol)) %>%
    mutate(conc = as.numeric(idx)) %>%
    mutate(idx = fct_reorder(idx, conc)) %>%
    ggplot(aes(x = mass, y = intensity, col = idx)) +
    geom_line() +
    scale_y_continuous(limits = c(0, NA)) +
    scale_color_viridis_d(end = 0.75, option = "C") +
    theme_minimal(base_size = 14) +
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
    errorbars = TRUE
  )
  p <- ggarrange(curve, peakProfile + labs(title = NULL))
  return(p)
}
