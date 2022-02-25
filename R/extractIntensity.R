#' Extract intensity using peaks as template
#'
#' @param peaks MALDIquant::MassPeaks list
#' @param spec  MALDIquant::MassSpectrum list
#' @param tol   numeric, tolerance in Da
#'
#' @return MALDIquant::MassPeaks list with extracted intensities from spec at m/z of peaks = pseudo peaks.
#' Useful in combination with sdMassSpectrum to get standard deviation of peaks as intensity matrix.
#' @export
extractIntensity <- function(peaks, spec, tol = 0.1) {
  if(!MALDIquant:::isMassPeaksList(peaks)) {

    if(!MALDIquant:::isMassPeaks(peaks)) {
      stop(sQuote("peaks"), " has to be MALDIquant::MassPeaks or list thereof!")
    }
  }

  if (length(peaks) != length(spec)) {
    if(MALDIquant:::isMassPeaks(peaks)) {
      res <- lapply(1:length(spec), FUN = function(i) {
        idx <- match.closest(x = mass(peaks),
                             table = mass(spec[[i]]),
                             tolerance = tol)
        createMassPeaks(mass = mass(peaks),
                        intensity = intensity(spec[[i]])[idx],
                        snr = rep(NA_integer_, length(idx)))
      })
      return(unlist(res))
    }
    stop("For each item in ", sQuote("peaks"), " there must be a spectrum in ",
         sQuote("spec"), "!")
  }

  res <- lapply(1:length(spec), FUN = function(i) {
    idx <- match.closest(x = mass(peaks[[i]]),
                         table = mass(spec[[i]]),
                         tolerance = tol)
    createMassPeaks(mass = mass(peaks[[i]]),
                    intensity = intensity(spec[[i]])[idx],
                    snr = rep(NA_integer_, length(idx)))
  })
  return(unlist(res))
}
