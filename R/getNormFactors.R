#' Get normalization factors from peak data.frame
#'
#' @param peaks         List of Object of class MALDIquant::MassPeaks 
#' @param targetMz      Numeric, target mass
#' @param tol           Numeric, tolerance around targetMz
#' @param tolppm        Logical, is the tolerance provided in ppm (TRUE) or Daltion (FALSE)
#' @param allowNoMatch  Logical, stop if targetMz is not fround in single spectrum?
#'                      If TRUE spectra without targetMz match will be excluded.
#'
#' @return              List with two entries:
#'                                       norm_factor The normalization factor for each spectrum
#'                                       specIdx     The index of the spectra with a match for targetMz
#'
#' @importFrom dplyr pull %>% filter arrange
#' @export
#' 
#' @examples
#' data(Blank2022peaks)
#' getNormFactors(Blank2022peaks, targetMz = 760.585, tol = 0.1, tolppm = FALSE)
getNormFactors <- function(peaks, targetMz, tol, tolppm = TRUE, allowNoMatch = TRUE) {
  allIdx <- 1:length(peaks)
  
  if (tolppm) {
    tol <- tol * targetMz / 1e6
  } 
  
  norm_factor <- map_vec(peaks, 
                         function(x) {
                           mz <- mass(x)
                           
                           idx <- match.closest(targetMz, 
                                                table = mz, 
                                                tolerance = tol, 
                                                nomatch = NA_integer_)
                           
                           return(intensity(x)[idx])
                         })
  specIdx <- which(!is.na(norm_factor))
  norm_factor <- na.omit(norm_factor)
  
  if (!all(allIdx %in% specIdx)) {
    if (!allowNoMatch) {
      stop("Could not find ", targetMz, " for all spectra! Consider adjusting tol.\n")
    }
    warning("Could not find ", targetMz, " in spectrum ", paste(which(!(allIdx %in% specIdx)), collapse = ", "), ".\n")
    specIdx <- which(allIdx %in% specIdx)
  } else {
    specIdx <- allIdx
  }
  if (length(specIdx) < 1) {
    stop("Could not find targetMz in any spectrum! Consider adjusting tol.\n")
  }
  return(list(
    norm_factor = norm_factor,
    specIdx = specIdx
  ))
}
