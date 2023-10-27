#' Preprocessing function to average spectra, detect and bin peaks before turning them into an intensity matrix
#'
#' @param spec                List of MALDIquant::MassSpectrum
#' @param averageMethod       Character, method for aggregation: "mean", "median" or "sum"
#' @param SNR                 Numeric, Signal noise value for peak detection
#' @param monoisotopicFilter  Logical, filter monoisotopic peaks
#' @param binTol              Numeric, tolerance for binning
#'
#' @return
#' List of lists with intensity matrix, average spectra and average peaks
.aggregateSpectra <- function(spec, averageMethod, SNR, monoisotopicFilter, binTol) {
  nm <- names(spec)
  stopifnot(!is.null(nm))
  stopifnot(isMassSpectrumList(spec))

  avg_spec <- averageMassSpectra(spec,
                                 labels = nm,
                                 method = averageMethod)

  cat(MALDIcellassay:::timeNow(),
      "building intensity matrix and applying variance filter... \n")
  peaks <- detectPeaks(avg_spec, method = "SuperSmoother", SNR = SNR)

  if(monoisotopicFilter) {
    cat(MALDIcellassay:::timeNow(),
        "Filtering monoisotopic peaks...\n")
    # set it to be less restrictive then default settings
    peaks <- monoisotopicPeaks(peaks,
                               size = 2L:10L,
                               minCor = 0.85,
                               tolerance = 1e-3)
  }

  peaksBinned <- binPeaks(peaks, tolerance = binTol)

  # perform variance filtering
  intmat <- intensityMatrix(peaksBinned, avg_spec)

  cat("      Found", dim(intmat)[2], "peaks in total.\n")
  rownames(intmat) <- names(avg_spec)

  return(list(intmat = intmat,
              avgPeaksBinned = peaksBinned,
              avgSpec = avg_spec))
}
