.getPQN <- function(l) {
  ## copied and modified from MALDIquant https://github.com/sgibb/MALDIquant/blob/master/R/calibrateIntensity-functions.R

  ## 1. normalization
  l <- calibrateIntensity(l, method="TIC")
  ## 2. median spectrum
  reference <- .averageMassSpectra(l, fun=.colMedians, mergeMetaData=FALSE)

  factor <- vapply(l,
                   function(x) {
                     ## 3. quotient
                     q <- approxfun(x)(reference@mass) / reference@intensity
                     ## 4. median of quotient
                     m <- median(q, na.rm=TRUE)
                   }, numeric(1))
  return(factor)
}


#' Normalize spectra and peaks
#'
#' @param spec     List of MALDIquant::MassSpectrum
#' @param peaks    List of MALDIquant::MassPeaks
#' @param normMeth Character, normalization method. Options are "TIC", "PQN", "median" and "mz".
#' @param normMz   Numeric, mz used to normalize.
#' @param normTol  Numeric, tolerance around normMz.
#'
#' @return
#' List of lists of normalized MALDIquant::MassSpectrum and normalized MALDIquant::MassPeaks
#'
#' @export
normalize <- function(spec, peaks, normMeth, normMz, normTol) {
  nm <- names(spec)
  stopifnot(!is.null(nm))
  stopifnot(is.numeric(as.numeric(nm)))

  switch(normMeth,
         "TIC" = {
           tic <- purrr::map_dbl(spec, totalIonCurrent)

           spec <- normalizeByFactor(spec, tic)
           peaks <- normalizeByFactor(peaks, tic)

           norm_fac <- tic
           included_specIdx <- 1:length(spec)
         },
         "PQN" = {
           pqn <- .getPQN(spec)

           spec <- normalizeByFactor(spec, pqn)
           peaks <- normalizeByFactor(peaks, pqn)

           norm_fac <- pqn
           included_specIdx <- 1:length(spec)
         },
         "median" = {
           median <- vapply(spec,
                            FUN = function(x)
                              MALDIquant:::.scalingFactor(object = x,
                                                          method = "median"),
                            numeric(1))

           spec <- normalizeByFactor(spec, median)
           peaks <- normalizeByFactor(peaks, median)

           norm_fac <- median
           included_specIdx <- 1:length(spec)
         },
         "mz" = {
           mzNorm <- getNormFactors(
             peaksdf = peaks2df(peaks),
             targetMz = normMz,
             tol = normTol,
             allowNoMatch = TRUE,
             tolppm = FALSE
           )

           norm_fac <- mzNorm$norm_factor
           spec <- normalizeByFactor(spec[mzNorm$specIdx], mzNorm$norm_factor)
           peaks <- normalizeByFactor(peaks[mzNorm$specIdx], mzNorm$norm_factor)
           included_specIdx <- mzNorm$specIdx

           u_nm <- unique(nm)
           u_fil <- unique(nm[included_specIdx])
           if(length(u_nm) != length(u_fil)) {
             # stop if a single condition got filtered completely

             label_removed <- u_nm[which(!(u_nm %in% u_fil))]

             stop("Could not find ", normMz, " in all spectra with label ",
                  paste0(label_removed, collapse = ", "),
                  ".\n Consider increasing tol.\n")
           }
         },
         "none" = {
           norm_fac <- list("norm_factor" = 0)
           included_specIdx <- 1:length(spec)
         }
  )

  return(list(spec = spec,
              peaks = peaks,
              factor = norm_fac,
              idx = included_specIdx))
}
