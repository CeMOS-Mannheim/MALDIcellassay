#' Shift mass axis
#'
#' @param spec    List of MALDIquant::MassSpectrum or MALDIquant::MassPeaks
#' @param mzdiff  Numeric vector, see getMzShift()
#'
#' @return List of MALDIquant::MassSpectrum or MALDIquant::MassPeaks with shifted mass axis.
#'
#' @importFrom MALDIquant isMassSpectrum isMassPeaks isMassSpectrumList isMassPeaksList
#' @export
shiftMassAxis <- function(spec, mzdiff) {
  if (isMassSpectrum(spec) || isMassPeaks(spec)) {
    spec@mass <- spec@mass + mzdiff
    return(spec)
  }

  if (isMassSpectrumList(spec) || isMassPeaksList(spec)) {
    if (!(length(spec) == length(mzdiff))) {
      stop("length(spec) != length(mzdiff) !\n")
    }
    for (i in 1:length(spec)) {
      spec[[i]]@mass <- spec[[i]]@mass + mzdiff[i]
    }
    return(spec)
  }
  stop("spec needs to be a MALDIquant::MassSpectrum or MALDIquant::MassPeaks or a list of these. \n")
}
