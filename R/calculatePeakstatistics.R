#' Convert numeric vector into proportions
#'
#' @param x Numeric vector
#'
#' @return
#' Numeric vector in proportions (0 - 1)
.absToProp <- function(x) {
  prop <- (x - min(x)) / (max(x) - min(x))

  return(prop)
}

#' Get value of curve for a give target effect size (e.g. 50%)
#'
#' @param model   `nplr::nplr` object
#' @param target  Numeric, typically 0.5 for the EC50 or IC50
#'
#' @return
#' Numeric, value on the fit curve.
#' Used to provide the correct target for `nplr::getEstimates()`
.getFitValue <- function(model, target) {
  x <- getYcurve(model)
  prop <- .absToProp(x)
  idx <- which.min(abs(prop - target))

  return(x[idx])
}


#' Calculate peak statistics
#'
#' @param curveFits       list of curve fits as returned by `MALDIcellassay::calculateCurveFit()`.
#' @param singlePeaks     list of MALDIquant::MassPeaks.
#' @param spec            list of MALDIquant::MassSpectrum.
#'
#' @return
#' A tibble with peak statistics.
#' @export
#' @importFrom dplyr join_by
#' @importFrom MALDIquant intensityMatrix isMassPeaksList

calculatePeakStatistics <- function(curveFits, singlePeaks, spec) {
  if (!is.list(curveFits)) {
    stop("curveFits must be a list of curve fits. See MALDIcellassay::calculateCurveFit().\n")
  }
  if (!isMassPeaksList(singlePeaks)) {
    stop("singlePeaks must be a list of MALDIquant::MassPeaks.\n")
  }
  if (is.null(names(singlePeaks))) {
    stop("singlePeaks must have concentrations as names!\n")
  }

  intensityMatrix <- intensityMatrix(singlePeaks, spec)
  rownames(intensityMatrix) <- names(singlePeaks)

  fit_df <- lapply(curveFits, function(x) {
    model <- x$model
    pIC50 <-
      -suppressMessages(
      suppressWarnings(
        getEstimates(model,
                     targets = .getFitValue(model, 0.5)))[, 3]
      )

    fc_window <- MALDIcellassay:::calculateFC(x$df)

    res_df <-
      suppressMessages(
      suppressWarnings(
        as_tibble(
          nplr::getGoodness(model)
        )
      )
    ) %>%
      mutate(
        fc_window = fc_window,
        pIC50 = pIC50
      )
    return(res_df)
  }) %>%
    bind_rows(.id = "mz") %>%
    rename("R2" = "gof")

  stat_df <-
    intensityMatrix %>%
    as_tibble() %>%
    mutate(sample = rownames(intensityMatrix)) %>%
    gather(mz, int, -sample) %>%
    arrange(as.numeric(mz)) %>%
    group_by(sample, mz) %>%
    summarise(
      min = min(int, na.rm = TRUE),
      mean = mean(int, na.rm = TRUE),
      max = max(int, na.rm = TRUE),
      stdev = sd(int, na.rm = TRUE),
      "cv%" = stdev / mean * 100
    ) %>%
    ungroup() %>%
    left_join(fit_df, by = join_by(mz)) %>%
    filter(!is.na(R2)) %>%
    mutate(mzIdx = as.numeric(as.factor(as.numeric(mz))))
  return(stat_df)
}
