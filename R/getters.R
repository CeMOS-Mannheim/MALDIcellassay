#' Extract the concentrations used in a MALDIassay
#'
#' @param object Object of class MALDIassay
#'
#' @return
#' Numeric vector, concentrations used in a MALDIassay
#' @export

getConc <- function(object) {
  if(!class(object) == "MALDIassay") {
    stop("object needs to be of class MALDIassay.")
  }
  return(object@settings$Conc)
}

#' Extract m/z used for normalization
#'
#' @param object Object of class MALDIassay
#'
#' @return
#' Numeric, m/z used for normalization 
#' @export

getNormMz <- function(object) {
  if(!class(object) == "MALDIassay") {
    stop("object needs to be of class MALDIassay.")
  }
  return(object@settings$normMz)
}

#' Extract tolerance used for normalization
#'
#' @param object Object of class MALDIassay
#'
#' @return
#' Numeric, tolerance used for normalization 
#' @export
getNormMzTol <- function(object) {
  if(!class(object) == "MALDIassay") {
    stop("object needs to be of class MALDIassay.")
  }
  return(object@settings$normTol)
}


#' Extract average spectra
#'
#' @param object Object of class MALDIassay
#'
#' @return
#' List of MALDIquantMassSpectrum
#' @export
getAvgSpectra <- function(object) {
  if(!class(object) == "MALDIassay") {
    stop("object needs to be of class MALDIassay.")
  }
  return(object@avgSpectra)
}

#' Extract peaks of average spectra
#'
#' @param object Object of class MALDIassay
#'
#' @return
#' List of MALDIquantMassPeaks
#' @export
getAvgPeaks <- function(object) {
  if(!class(object) == "MALDIassay") {
    stop("object needs to be of class MALDIassay.")
  }
  return(object@avgPeaks)
}

#' Extract normalization method
#'
#' @param object Object of class MALDIassay
#'
#' @return
#' Character, normalization method used.
#' @export
getNormMethod <- function(object) {
  if(!class(object) == "MALDIassay") {
    stop("object needs to be of class MALDIassay.")
  }
  return(object@settings$normMeth)
}

