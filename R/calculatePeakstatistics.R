#' Calculate peak statistics
#'
#' @param curveFits       list of curve fits as returned by ```MALDIcellassay::calculateCurveFit()```
#' @param intensityMatrix matrix of intensities with rownames as concentrations of the respective spectra
#'
#' @return
#' A tibble with peak statistics
#' @export
#' @importFrom dplyr join_by

calculatePeakStatistics <- function(curveFits, intensityMatrix) {
  if (!is.list(curveFits)) {
    stop("curveFits must be a list of curve fits. See MALDIcellassay::calculateCurveFit().\n")
  }
  if (!is.matrix(intensityMatrix)) {
    stop("intensityMatrix must be of class matrix. See MALDIquant::intensityMatrix().\n")
  }
  if (is.null(rownames(intensityMatrix))) {
    stop("intensityMatrix must have concentrations as rownames!\n")
  }

  fit_df <- lapply(curveFits, function(x) {
    model <- x$model
    pIC50 <- -suppressMessages(
      suppressWarnings(
        getEstimates(model, targets = 0.5)[, 3]
      )
    )
    fc_window <- MALDIcellassay:::calculateFC(x$df)
    res_df <- suppressMessages(
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

  stat_df <- intensityMatrix %>%
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
