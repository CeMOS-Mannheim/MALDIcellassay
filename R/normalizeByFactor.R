#' Apply normalization factors to spectra
#'
#' @param spec         List of MALDIquant::MassSpectrum or MALDIquant::MassPeaks
#' @param factors      Numeric vector of normalization factors. See getNormFactors().
#'
#' @return             List of normalized Spectra or Peaks
#'
#' @importFrom MALDIquant intensity
#' @export
normalizeByFactor <- function(spec, factors) {
  if(!length(spec) == length(factors)) {
    stop("Number of spectra and normalization factors not equal!\n")
  }

  spec_res <- spec

  for(i in 1:length(factors)) {
    spec_res[[i]]@intensity <- intensity(spec_res[[i]])/factors[i]
  }
  return(spec_res)
}
