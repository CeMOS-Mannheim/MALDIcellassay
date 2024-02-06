#' Detect peaks
#'
#' @param spec    MALDIquant::MassSpectrum or list thereof
#' @param SNR     Numeric, signal to noise value
#' @param method  Character, method see MALDIquant::detectPeaks
#' @param halfWindowSize Numeric, defines width of window for peak detection. See `MALDIquant::detectPeaks()`.
#'
#' @details
#' Just a wrapper around MALDIquant::detectPeaks to ensure that the returned peak list is named.
#'
#' @return
#' List of MALDIquant::MassPeaks with the same names as `spec`
.detectPeaks <- function(spec, SNR, method = "SuperSmoother", halfWindowSize = 3) {
  stopifnot(isMassSpectrum(spec)| isMassSpectrumList(spec))

  peaks <- detectPeaks(spec,
                       SNR = SNR,
                       method = "SuperSmoother",
                       halfWindowSize = halfWindowSize)
  names(peaks) <- names(spec)

  return(peaks)
}
