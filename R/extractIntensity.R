#' Extract intensity using peaks as template
#'
#' @param peaks MALDIquant::MassPeaks list
#' @param spec  MALDIquant::MassSpectrum list
#' @param tol   numeric, tolerance in Da
#'
#' @return MALDIquant::MassPeaks list with extracted intensities from spec at m/z of peaks = pseudo peaks.
#' Useful in combination with sdMassSpectrum to get standard deviation of peaks as intensity matrix.
#' @export
extractIntensity <- function(mz, peaks, spec, tol) {
  if(!(length(peaks) == length(spec))) {
    stop("length of peaks and spec must match.\n")
  }

  res_peaks <-purrr::map(seq_along(peaks),
                         function(i) {
                           #browser()
                           mz_peaks <- mass(peaks[[i]])
                           mz_spec <- mass(spec[[i]])

                           peak_mz_idx <- match.closest(x = mz,
                                                        table = mz_peaks,
                                                        tol = tol)
                           spec_mz_idx <- match.closest(x = mz,
                                                        table = mz_spec,
                                                        tol = Inf)

                           # extract intensity from peaks is possible
                           # otherwise use noise from spectra
                           # resulting in a dense intensity matrix
                           new_int <- ifelse(!is.na(peak_mz_idx),
                                             yes = intensity(peaks[[i]])[peak_mz_idx],
                                             no = intensity(spec[[i]])[spec_mz_idx])

                           new_snr <- ifelse(!is.na(peak_mz_idx),
                                             yes = snr(peaks[[i]])[peak_mz_idx],
                                             no = NA_integer_)

                           createMassPeaks(mass = mz,
                                           intensity = new_int,
                                           snr = new_snr,
                                           metaData = metaData(peaks[[i]]))
                         })
  names(res_peaks) <- names(peaks)
  return(res_peaks)
}


# extractIntensity <- function(peaks, spec, tol = 0.1) {
#   if (!MALDIquant:::isMassPeaksList(peaks)) {
#     if (!MALDIquant:::isMassPeaks(peaks)) {
#       stop(sQuote("peaks"), " has to be MALDIquant::MassPeaks or list thereof!")
#     }
#   }
#
#   if (length(peaks) != length(spec)) {
#     if (MALDIquant:::isMassPeaks(peaks)) {
#       res <- lapply(1:length(spec), FUN = function(i) {
#         idx <- match.closest(
#           x = mass(peaks),
#           table = mass(spec[[i]]),
#           tolerance = tol
#         )
#         int <- intensity(spec[[i]])[idx]
#         if(any(is.na(int))) {
#           na_idx <- which(is.na(int))
#           warning("Could not match intensity for spectrum idx",
#                   i, "and m/z", round(mass(peaks)[na_idx], 2),
#                   ".\n Replacing by zero.\n")
#           int[na_idx] <- 0
#         }
#         createMassPeaks(
#           mass = mass(peaks),
#           intensity = int,
#           snr = rep(NA_integer_, length(idx))
#         )
#       })
#       return(unlist(res))
#     }
#     stop(
#       "For each item in ", sQuote("peaks"), " there must be a spectrum in ",
#       sQuote("spec"), "!"
#     )
#   }
#
#   res <- lapply(1:length(spec), FUN = function(i) {
#     idx <- match.closest(
#       x = mass(peaks[[i]]),
#       table = mass(spec[[i]]),
#       tolerance = tol
#     )
#     createMassPeaks(
#       mass = mass(peaks[[i]]),
#       intensity = intensity(spec[[i]])[idx],
#       snr = rep(NA_integer_, length(idx)),
#       metaData = metaData(spec)
#     )
#   })
#   return(unlist(res))
# }
