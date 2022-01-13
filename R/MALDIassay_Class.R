#' Class MALDIassay
#'
#' A class for holding MALDI assay related information.
#'
#' @param object MALDIassay.
#'
#' @importFrom dplyr arrange group_by first summarise
#' @importFrom MALDIquant mass


setClass("MALDIassay",
         representation = representation(avgSpectra = "list",
                                         singlePeaks = "list",
                                         normFactors = "numeric",
                                         mzShifts = "numeric",
                                         fits = "list",
                                         stats = "data.frame",
                                         settings = "list")
)

setMethod("show", signature(object = "MALDIassay"),
          function(object) {

            mz <- round(object@settings$normMz, digits = 2)
            tol <-  object@settings$normTol
            numPeaksTotal <- length(mass(object@singlePeaks[[1]]))
            hiVarPeaks <- length(unique(object@stats$mz))
            # Compose normalization information
            if(normMeth == "mz") {
              normStr <- paste("Normalization on m/z", mz, "±", tol, "Da.\n")
            } else {
              normStr <- paste("Normalization using", normMeth, "method.\n")
            }


            cat("------MALDIassay object------\n")
            cat("\n")
            cat(normStr)
            cat("\n")
            if(object@settings$SinglePointRecal) {
              cat("Single point recalibation on", mz, "with", tol, "Da tolerance.\n")
              cat("Avg. mass shift:", round(mean(object@mzShifts), 4), "Da. Max abs. shift:", round(max(abs(object@mzShifts)), 4), "Da.\n")
              cat("\n")
            }

            cat("Found", numPeaksTotal, "peaks above SNR", object@settings$SNR,  "and", hiVarPeaks, "high variance peaks.\n")
            cat("\n")

            cat("Top-features based on R² and max/min-Fold-chage:\n")
            print(object@stats %>%
              group_by(mz) %>%
              summarise(R2 = first(R2),
                        wgof = first(wgof),
                        FC = first(fc_window)) %>%
              arrange(desc(R2), desc(FC)))
          })
