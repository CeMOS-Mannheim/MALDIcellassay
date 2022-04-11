#' generate ggplot objects for each of the curve fits in a MALDIassay object
#'
#' @param object    object of class MALDIassay
#' @param fc_thresh numeric, max/min fold-change above which plots should be generated
#' @param markValue numeric, a value to add as a reference to the plots
#' @param R2_thresh numeric, min. R-squared (goodness of curve fit) to plot curve
#' @param mzIdx     numeric, indicies of mz values to plot (see \code{getPeakStatistics()}). Note, fc_thresh and R2_thresh filters do not apply if mzIdx is set!
#' @param errorbars logical, add errorbars to plot representing standard deviation in regards to the measurment replicates.
#'
#' @return
#' list of ggplot objects
#'
#' @importFrom ggplot2 geom_errorbar geom_point ggplot aes geom_line scale_x_continuous theme_bw theme element_text labs
#' @importFrom dplyr group_by summarise
#' @importFrom tibble tibble
#' @importFrom nplr getGoodness getEstimates getXcurve getYcurve getX getY convertToProp
#' @export
plotCurves <- function(object, fc_thresh = 1, R2_tresh = 0, markValue = NA, mzIdx = NULL, errorbars = FALSE) {
  stopIfNotIsMALDIassay(object)
  if (is.null(mzIdx)) {
    res_list <- getCurveFits(object)
  } else {
    res_list <- getCurveFits(object)[mzIdx]
  }
  len <- length(names(res_list))
  mz_vals <- as.numeric(names(res_list))

  p_list <- vector("list", length = len)
  for (i in 1:len) {
    mz <- mz_vals[i]
    model <- res_list[[as.character(mz)]]$model
    df <- res_list[[as.character(mz)]]$df

    ic50 <- 10^getEstimates(model, targets = 0.5)[, 3]
    min <- min(df$value)
    max <- max(df$value)

    fc_window <- MALDIcellassay:::calculateFC(df)
    R2 <- getGoodness(model)[[1]]

    if ((abs(fc_window) >= fc_thresh & R2 >= R2_tresh) | !is.null(mzIdx)) {
      df_C <- tibble(xC = getXcurve(model), yC = getYcurve(model))
      df_P <- tibble(x = getX(model), y = getY(model))

      int <- vapply(getSinglePeaks(object), function(x) {
        targetmass <- mz
        mass <- mass(x)
        idx <- match.closest(targetmass, mass, tolerance = 0.01)
        int <- intensity(x)
        return(int[idx])
      }, numeric(1))

      df_singlePeaks <- tibble(
        conc = getConc(object),
        int_raw = int,
        int = convertToProp(
          y = int,
          T0 = min,
          Ctrl = max
        )
      ) %>%
        group_by(conc) %>%
        summarise(
          y = mean(int),
          ymin = y - sd(int),
          ymax = y + sd(int)
        )

      p <- ggplot(data = df_P, aes(x = x, y = y)) +
        geom_line(data = df_C, aes(x = xC, y = yC)) +
        geom_point() +
        scale_x_continuous(labels = c(0, 10^df_P$x[-1]), breaks = df_P$x) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
          x = "Conc.",
          y = "relative Int. [% of max Int.]",
          title = paste0(
            "mz ", round(mz, 2), " Da, R\u00B2=", round(R2, 3), "\n",
            "pIC50=", round(-log10(ic50), 3),
            " min=", round(min, 3),
            " max=", round(max, 3),
            " FC=", round(fc_window, 2)
          )
        )

      if (errorbars) {
        p <- p +
          geom_errorbar(
            data = df_singlePeaks, aes(
              x = log10(conc),
              y = y,
              ymin = ymin,
              ymax = ymax
            ),
            alpha = 0.5
          )
      }

      if (!is.na(markValue)) {
        p <- p + geom_vline(aes(xintercept = markValue), linetype = "dashed")
      }
      p_list[[i]] <- p
      names(p_list) <- mz_vals
    }
  }
  # check for empty entries (result of filtering for FC or R2) and remove them
  idx <- vapply(p_list, function(x) {
    length(x) > 0
  }, FUN.VALUE = TRUE)
  if (sum(!idx) == len) {
    stop("Nothing to plot. Condsider decreasing fc_thresh or R2_tresh.\n")
  }

  if (length(p_list) == 1) {
    return(p_list[[1]])
  }

  return(p_list[idx])
}
