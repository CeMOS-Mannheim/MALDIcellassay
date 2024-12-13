#' Class MALDIassay
#'
#' A class for holding MALDI assay related information.
#'
#' @param object MALDIassay.
#'
#' @importFrom dplyr arrange group_by first summarise desc
#' @importFrom MALDIquant mass
#' @importFrom utils head
setClass("MALDIassay",
  representation = representation(
    avgSpectra = "list",
    avgPeaks = "list",
    singlePeaks = "list",
    singleSpecSpots = "character",
    normFactors = "numeric",
    mzShifts = "numeric",
    fits = "list",
    stats = "data.frame",
    included_specIdx = "numeric",
    settings = "list"
  )
)

.validMALDIassay <- function(object) {
  if (length(object@avgSpectra) < 1) {
    return("Length of avgSpectra can't be 0.")
  }
  
  if (!MALDIquant:::.isMassObjectList(object@avgSpectra)) {
    return("avgSpectra must be a list of class MALDIquant::MassSpectrum or MALDIquant::MassPeaks objects.")
  }
  
  if (length(object@avgPeaks) < 1) {
    return("Length of avgPeaks can't be 0.")
  }
  
  if (!isMassPeaksList(object@avgPeaks)) {
    return("avgPeaks must be a list of class MALDIquant::MassPeaks objects.")
  }
  
  if (length(object@singlePeaks) < 1) {
    return("Length of singlePeaks can't be 0.")
  }
  
  if (!isMassPeaksList(object@singlePeaks)) {
    return("singlePeaks must be a list of class MALDIquant::MassPeaks objects.")
  }
  
  if (length(object@avgSpectra) != length(object@avgPeaks)) {
    return(paste0("avgSpectra (", length(object@avgSpectra),
                  ") and avgPeaks (", length(object@avgPeaks),
                  ") must have the same length."))
  }
  
  if (length(object@singlePeaks) != length(object@included_specIdx)) {
    return(paste0("singlePeaks (", length(object@singlePeaks),
                  ") and included_specIdx (", length(object@included_specIdx),
                  ") must have the same length."))
  }
  
  if (!(length(object@singlePeaks) == length(object@mzShifts) || length(object@mzShifts) == 1)) {
    return(paste0("singlePeaks (", length(object@singlePeaks),
                  ") and mzShifts (", length(object@mzShifts),
                  ") must have the same length or mzShifts must have length 1."))
  }
  
  if (!(length(object@singlePeaks) == length(object@normFactors) || length(object@normFactors) == 1)) {
    return(paste0("singlePeaks (", length(object@singlePeaks),
                  ") and normFactors (", length(object@normFactors),
                  ") must have the same length or normFactors must have length 1."))
  }
  
  if (!(length(object@singlePeaks) == length(object@singleSpecSpots) || length(object@singleSpecSpots) == 1)) {
    return(paste0("singlePeaks (", length(object@singlePeaks),
                  ") and singleSpecSpots (", length(object@singleSpecSpots),
                  ") must have the same length or singleSpecSpots must have length 1."))
  }
  
  if (length(object@fits) != length(unique(object@stats[["mz"]]))) {
    return(paste0("fits (", length(object@fits),
                  ") and the unique m/z values in stats (", length(unique(object@stats[["mz"]])),
                  ") must have the same length."))
  }
  
  TRUE
}

setValidity("MALDIassay", method=.validMALDIassay)

#' Print summary of MALDIassay object
#'
#' @param object object of class MALDIassay
#'
#' @return nothing, prints to console
#'
#' @noRd
show_MALDIassay <- function(object) {
  mz <- round(getNormMz(object), digits = 2)
  varFilterMethod <- getVarFilterMethod(object)
  tol <- getNormMzTol(object)
  numPeaksTotal <- length(mass(getSinglePeaks(object)[[1]]))
  hiVarPeaks <- length(unique(getPeakStatistics(object)$mz))
  mzdev <- getRecalibrationError(object)
  meanMzShift <- round(mean(getAppliedMzShift(object)), 4)
  absMaxMzShift <- round(max(abs(getAppliedMzShift(object))), 4)
  conc <- getConc(object)

  # Compose normalization information
  if (getNormMethod(object) == "mz") {
    normStr <- paste("Normalization on m/z", mz, "\u00B1", tol, "Da.\n")
  } else {
    normStr <- paste("Normalization using", getNormMethod(object), "method.\n")
  }


  cat("------MALDIassay object------\n")
  cat("\n")
  cat("Including", length(unique(conc)), "concentrations,\n")
  cat("ranging from", min(conc), "to", max(conc), ".")
  cat("\n")
  cat(normStr)
  cat("\n")
  if (object@settings$SinglePointRecal) {
    cat("Single point recalibation on", mz, "with", tol, "Da tolerance.\n")
    cat("Avg. mass shift before recal.:", meanMzShift, "Da. Max abs. shift:", absMaxMzShift, "Da.\n")
    cat("Avg. mass shift after recal. :", round(mzdev$mean, 4), "Da. Max abs. shift:", round(mzdev$meanAbs, 4), "Da.\n")
    cat("\n")
  }

  cat(paste0("Found ", numPeaksTotal, " peaks (SNR ", getSNR(object), ") and ", hiVarPeaks, " high variance peaks\n"))
  cat("using variance filtering method:", paste0(varFilterMethod, ".\n"))
  cat("\n")

  cat("Top5-features based on Fold-Change and R\u00B2:\n")
  print(getPeakStatistics(object, summarise = TRUE) %>%
    arrange(desc(.data$CRS), desc(.data$log2FC)) %>%
    as.data.frame() %>%
    head(n = 5))
}

setMethod(
  "show", signature(object = "MALDIassay"),
  show_MALDIassay
)
