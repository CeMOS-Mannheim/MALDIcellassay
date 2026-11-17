#'  Plot a peak of interest from a MALDIassay object
#'
#'  Draws the region around the peak of interest using the averaged spectra. For
#'  profile data the spectra are shown as curves, for centroided data the
#'  centroids are shown as sticks.
#'
#'  With `showPeaks = TRUE` the peaks that are used for curve fitting are marked
#'  with open circles: the peaks detected on the averaged spectra for profile
#'  data, and the centroids themselves for centroided data.
#'
#' @param object object of class MALDIassay
#' @param mzIdx  numeric, index of mass of interest (see \code{getPeakStatistics()})
#' @param tol    numeric, tolerance around peak to plot
#' @param centroided logical, whether the data are centroided. If \code{TRUE}
#'   the centroids are drawn as sticks instead of profile curves.
#' @param showPeaks logical, whether the peaks used for curve fitting are marked
#'   with open circles (default \code{TRUE}).
#'
#' @return
#' ggplot object
#'
#' @importFrom ggplot2 ggplot aes geom_line scale_y_continuous scale_color_viridis_d labs geom_rect geom_point geom_segment
#' @importFrom dplyr filter between bind_rows .data
#' @importFrom tibble tibble
#' @importFrom MALDIquant mass intensity detectPeaks
#' @importFrom forcats fct_reorder
#' @export
#' @examples
#'
#' data(Blank2022res)
#' plotPeak(Blank2022res, mzIdx = 2)
#' 
#' # centroided data: centroids are drawn as sticks
#' plotPeak(Blank2022res, mzIdx = 2, centroided = TRUE)
plotPeak <- function(object, mzIdx, tol = 0.8, centroided = FALSE,
                     showPeaks = TRUE) {
  if (missing(mzIdx) | is.null(mzIdx) |is.na(mzIdx)) {
    stop("No mzIdx supplied.\n")
  }

  if(tol <= 0 | is.null(tol) | is.na(tol)) {
    stop("tol needs to be a positive number.\n")
  }

  if (!is.logical(centroided) || length(centroided) != 1L || is.na(centroided)) {
    stop("centroided needs to be TRUE or FALSE.\n")
  }

  if (!is.logical(showPeaks) || length(showPeaks) != 1L || is.na(showPeaks)) {
    stop("showPeaks needs to be TRUE or FALSE.\n")
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
    filter(between(.data$mass, mz - tol, mz + tol)) %>%
    mutate(conc = as.numeric(.data$conc)) %>%
    mutate(conc = fct_reorder(scales::scientific(.data$conc), .data$conc))

  # peaks that end up in the curve fitting
  peaks_df <- NULL
  if (showPeaks) {
    if (centroided) {
      # the centroids are the data
      peaks_df <- df
    } else {
      detected <- detectPeaks(spec, SNR = getSNR(object),
                              method = "SuperSmoother")
      # detectPeaks() drops the names of the spectra, keep the concentrations
      names(detected) <- names(spec)

      peaks_df <- map(detected, function(x) {
        tibble(
          mass = mass(x),
          intensity = intensity(x)
        )
      }) %>%
        bind_rows(.id = "conc") %>%
        filter(between(.data$mass, mz - tol, mz + tol)) %>%
        mutate(conc = as.numeric(.data$conc)) %>%
        mutate(conc = fct_reorder(scales::scientific(.data$conc), .data$conc))
    }
  }

  ymax <- if (nrow(df) > 0) {
    max(pull(df, .data[["intensity"]])) * 1.05
  } else {
    1
  }

  kind <- if (centroided) "Centroids" else "Profile"
  title <- paste0(kind, " of m/z ", round(mz, 2), " \u00B1 ", round(tol, 2), "Da")

  p <-ggplot() +
    geom_rect(aes(xmin = mz - getBinTol(object) * mz,
                  xmax = mz + getBinTol(object) * mz,
                  ymin = 0,
                  ymax = ymax),
              alpha=0.2,
              fill="black")

  if (centroided) {
    p <- p +
      geom_segment(data = df,
                   aes(x = .data[["mass"]], xend = .data[["mass"]],
                       y = 0, yend = .data[["intensity"]],
                       col = .data[["conc"]]),
                   linewidth = 0.4, alpha = 0.85)
  } else {
    p <- p +
      geom_line(data = df,
                aes(x = .data[["mass"]], y = .data[["intensity"]],
                    col = .data[["conc"]]))
  }

  if (!is.null(peaks_df) && nrow(peaks_df) > 0) {
    p <- p +
      geom_point(data = peaks_df,
                 aes(x = .data[["mass"]], y = .data[["intensity"]],
                     col = .data[["conc"]]),
                 shape = 1, size = 3, stroke = 1, show.legend = FALSE)
  }

  p <- p +
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
