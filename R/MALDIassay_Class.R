#' Class MALDIassay
#'
#' A class for holding MALDI assay related information.
#'
#' @param object MALDIassay.
#'
#' @importFrom dplyr arrange group_by first summarise desc
#' @importFrom MALDIquant mass


setClass("MALDIassay",
         representation = representation(avgSpectra = "list",
                                         avgPeaks = "list",
                                         singlePeaks = "list",
                                         normFactors = "numeric",
                                         mzShifts = "numeric",
                                         fits = "list",
                                         stats = "data.frame",
                                         settings = "list")
)



show_MALDIassay <- function(object) {
  
  mz <- round(getNormMz(object), digits = 2)
  varFilterMethod <- object@settings$varFilterMethod
  tol <-  getNormMzTol(object)
  numPeaksTotal <- length(mass(object@singlePeaks[[1]]))
  hiVarPeaks <- length(unique(object@stats$mz))
  mzdev <- getRecalibrationError(object)
  # Compose normalization information
  if(getNormMethod(object) == "mz") {
    normStr <- paste("Normalization on m/z", mz, "±", tol, "Da.\n")
  } else {
    normStr <- paste("Normalization using", getNormMethod(object), "method.\n")
  }
  
  
  cat("------MALDIassay object------\n")
  cat("\n")
  cat("Including", length(getAvgPeaks(object)), "different concentrations.\n")
  cat(normStr)
  cat("\n")
  if(object@settings$SinglePointRecal) {
    cat("Single point recalibation on", mz, "with", tol, "Da tolerance.\n")
    cat("Avg. mass shift before recal.:", round(mean(object@mzShifts), 4), "Da. Max abs. shift:", round(max(abs(object@mzShifts)), 4), "Da.\n")
    cat("Avg. mass shift after recal. :", round(mzdev$mean, 4), "Da. Max abs. shift:", round(mzdev$meanAbs, 4), "Da.\n")
    cat("\n")
  }
  
  cat("Found", numPeaksTotal, "peaks above SNR", object@settings$SNR,  "and", hiVarPeaks, "high variance peaks\n")
  cat("using variance filtering method:", paste0(varFilterMethod, ".\n"))
  cat("\n")
  
  cat("Top5-features based on Fold-chage and R²:\n")
  print(getPeakStatistics(object) %>%
          mutate(mz = round(as.numeric(mz), 3)) %>%
          group_by(mz) %>%
          summarise(mzIdx = first(mzIdx),
                    R2 = dplyr::first(round(R2,4)),
                    wgof = dplyr::first(round(wgof,4)),
                    FC = dplyr::first(round(fc_window, 4))) %>%
          arrange(desc(FC), desc(R2)) %>%
          as.data.frame() %>%
          head(n = 5))
}

setMethod("show", signature(object = "MALDIassay"),
          MALDIcellassay:::show_MALDIassay)


