#' Detect peaks
#'
#' @param spec    MALDIquant::MassSpectrum or list thereof
#' @param SNR     Numeric, signal to noise value
#' @param method  Character, method see MALDIquant::detectPeaks
#'
#' @details
#' Just a wrapper around MALDIquant::detectPeaks to ensure that the returned peak list is named.
#'
#' @return
#' List of MALDIquant::MassPeaks with the same names as `spec`
.detectPeaks <- function(spec, SNR, method = "SuperSmoother") {
  stopifnot(isMassSpectrum(spec)| isMassSpectrumList(spec))

  peaks <- detectPeaks(spec, SNR = SNR, method = "SuperSmoother")
  names(peaks) <- names(spec)

  return(peaks)
}
