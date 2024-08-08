#' Column wise standard deviation
#'
#'
#' also a function of sgibb
#' see https://stackoverflow.com/questions/17549762/is-there-such-colsd-in-r
#' @param x     matrix
#' @param na.rm logical
#'
#' @return sd
#' @export
colSdColMeans <- function(x, na.rm = TRUE) {
  if (na.rm) {
    n <- colSums(!is.na(x)) # thanks @flodel
  } else {
    n <- nrow(x)
  }
  colVar <- colMeans(x * x, na.rm = na.rm) - (colMeans(x, na.rm = na.rm))^2
  return(sqrt(colVar * n / (n - 1)))
}

#' Still fork from sgibb MALDIquant
#'
#' @param l       see MALDIquant
#' @param mergeMetaData see MALDIquant
#'
#' @return see MALDIquant
sdMassSpectraFun <- function(l, mergeMetaData = TRUE) {

  ## merge metaData
  if (mergeMetaData) {
    metaData <- MALDIquant:::.mergeMetaData(lapply(l, function(x) x@metaData))
  } else {
    metaData <- list()
  }

  ## use the first non empty spectrum as reference
  i <- which(!vapply(l, isEmpty, logical(1L)))[1L]
  if (!is.na(i)) {
    mass <- l[[i]]@mass
  } else {
    mass <- NA_real_
  }

  ## interpolate not existing masses
  approxSpectra <- lapply(l, MALDIquant:::approxfun)

  ## get interpolated intensities
  intensityList <- lapply(approxSpectra, function(x) x(mass))

  ## create a matrix which could merged
  m <- do.call(rbind, intensityList)

  ## merge intensities
  intensity <- colSdColMeans(m, na.rm = TRUE)

  ## create an empty spectrum if all intensities are NaN
  if (is.nan(intensity[1L])) {
    intensity <- double()
    mass <- double()
  }

  createMassPeaks(mass = mass, intensity = intensity, snr = rep(NA_integer_, length(intensity)),metaData = metaData)
}





#' Compute standard-deviation spectra
#'
#' This is a fork from sgibb's MALDIquant::averageMassSpectra() function.
#' It is now able to compute "standard-deviation spectra".
#'
#' @param l      list, list of MassSpectrum objects.
#' @param labels list, list of factors (one for each MassSpectrum object) to do groupwise averaging.
#' @param ...    arguments to be passed to underlying functions (currently only mc.cores is supported).
#'
#' @return
#' Returns a single (no labels given) or a list (labels given) of standard-deviation spectra as MassSpectrum objects.
#' @export

sdMassSpectrum <- function(l, labels, ...) {

  #MALDIquant:::.stopIfNotIsMassSpectrumList(l)
  MALDIquant:::.doByLabels(
    l = l, labels = labels, FUN = sdMassSpectraFun,
  )
}
