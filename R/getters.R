#' Extract the concentrations used in a MALDIassay
#'
#' @param object Object of class MALDIassay
#'
#' @return
#' Numeric vector, concentrations used in a MALDIassay
#' @export

getConc <- function(object) {
  stopIfNotIsMALDIassay(object)
  return(object@settings$Conc[object@included_specIdx])
}

#' Extract m/z used for normalization
#'
#' @param object Object of class MALDIassay
#'
#' @return
#' Numeric, m/z used for normalization
#' @export

getNormMz <- function(object) {
  stopIfNotIsMALDIassay(object)
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
  stopIfNotIsMALDIassay(object)
  return(object@settings$normTol)
}

#' Extract SNR used for peak detection
#'
#' @param object Object of class MALDIassay
#'
#' @return
#' Numeric, SNR used for peak detection
#' @export
getSNR <- function(object) {
  stopIfNotIsMALDIassay(object)
  return(object@settings$SNR)
}

#' Extract average spectra
#'
#' @param object Object of class MALDIassay
#'
#' @return
#' List of MALDIquantMassSpectrum
#' @export
getAvgSpectra <- function(object) {
  stopIfNotIsMALDIassay(object)
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
  stopIfNotIsMALDIassay(object)
  return(object@avgPeaks)
}

#' Extract peaks of single spectra spectra (before average calculation)
#'
#' @param object Object of class MALDIassay
#'
#' @return
#' List of MALDIquantMassPeaks
#' @export
getSinglePeaks <- function(object) {
  stopIfNotIsMALDIassay(object)
  return(object@singlePeaks)
}


#' Extract normalization method
#'
#' @param object Object of class MALDIassay
#'
#' @return
#' Character, normalization method used.
#' @export
getNormMethod <- function(object) {
  stopIfNotIsMALDIassay(object)
  return(object@settings$normMeth)
}

#' Extract variance filtering method
#'
#' @param object Object of class MALDIassay
#'
#' @return
#' Character of variance filtering method used
#' @export
getVarFilterMethod <- function(object) {
  stopIfNotIsMALDIassay(object)
  return(object@settings$varFilterMethod)
}

#' Extract applied mz-shift
#'
#' @param object Object of class MALDIassay
#'
#' @return
#' Numeric vector of mz-shits applied to spectra
#' @export
getAppliedMzShift <- function(object) {
  stopIfNotIsMALDIassay(object)
  return(object@mzShifts)
}

#' Extract curve fits
#'
#' @param object Object of class MALDIassay
#'
#' @return
#' List, containing the data used to do the fits as well as the nlpr curve fit .
#' @export
getCurveFits <- function(object) {
  stopIfNotIsMALDIassay(object)
  return(object@fits)
}

#' Extract directory path
#'
#' @param object Object of class MALDIassay
#'
#' @return
#' List, containing the data used to do the fits as well as the nlpr curve fit .
#' @export
getDirectory <- function(object) {
  stopIfNotIsMALDIassay(object)
  return(object@settings$dir)
}

#' Extract peak statistics
#'
#' @param object Object of class MALDIassay
#' @param summarise Logical, return summarised results (one result per mz and not per mz and spectra)
#'
#' @return
#' A tibble with peak statistics (R², fold-change, CV%, etc.)
#' @export
getPeakStatistics <- function(object, summarise = FALSE) {
  stopIfNotIsMALDIassay(object)
  stats <- object@stats

  if (summarise) {
    stats <- stats %>%
      mutate(mz = round(as.numeric(mz), 3)) %>%
      group_by(mz, mzIdx) %>%
      summarise(
        pIC50 = first(pIC50),
        R2 = first(R2),
        wgof = first(wgof),
        FC = first(fc_window)
      )
  }

  return(stats)
}

#' Calculate remaining calibration error of a MALDIassay object
#'
#' @param object Object of class MALDIassay
#'
#' @return
#' A tibble containing statistics about remaining calibration error
#' @export
#'
#' @importFrom tibble tibble
getRecalibrationError <- function(object) {
  peaks <- peaks2df(getAvgPeaks(object))
  mzdev <- getMzShift(
    peaksdf = peaks,
    tol = getNormMzTol(object),
    targetMz = getNormMz(object),
    tolppm = FALSE,
    allowNoMatch = FALSE
  )

  res_df <- tibble(
    meanAbs = mean(abs(mzdev$mzshift)),
    sdAbs = sd(abs(mzdev$mzshift)),
    maxAbs = max(abs(mzdev$mzshift)),
    mean = mean(mzdev$mzshift),
    sd = sd(mzdev$mzshift)
  )
  return(res_df)
}

#' Get the mz value associated with a mzIdx
#'
#' @param object Object of class MALDIassay
#' @param mzIdx  numeric, index of mass of interest (see \code{getPeakStatistics()})
#'
#' @return
#' numeric, mz value
#' @export
getMzFromMzIdx <- function(object, mzIdx) {
  stopIfNotIsMALDIassay(object)
  mz <- getPeakStatistics(object, TRUE) %>%
    filter(mzIdx == {{mzIdx}}) %>%
    pull(mz)

  if(length(mz)>1) {
    warning("Something is wrong. There are multiple mz-values with this index!\n")
  }
  return(mz)
}
