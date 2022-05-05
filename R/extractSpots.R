#' Extract the spot coordinates
#'
#' @param spec list of MALDIquant::MassSpectrum or MALDIquant::MassPeaks objects
#'
#' @return
#' Character vector of spot names. If multiple spots are used (e.g. for average spectra) they will be concatenate.
extractSpots <- function(spec) {
  if(!isMassSpectrumList(spec)) {
    if(!isMassPeaksListt(spec)) {
      stop("spec needs to be either a list of MassSpectra or MassPeaks.\n")
    }
  }
  spot <- vapply(spec,
                 function(x) {
                   paste(metaData(x)$patch, collapse = " ")
                   }, FUN.VALUE = character(1))
  return(spot)
}
