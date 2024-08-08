#' Filter for high variance signals
#'
#' @param vars    Numeric vector, variances of signals
#' @param method  Character, filtering method. One of "mean" (default), "median", "q25", "q75" (25 and 75% quantile) or "none".
#'
#' @return Indices of spectra with a high variance
#'
#' @export
#' 
#' @importFrom stats median quantile
filterVariance <- function(vars, method = c("mean", "median", "q25", "q75", "none")) {
  method <- match.arg(method)
  switch(method,
    "mean" = {
      idx <- which(vars > mean(vars))
      cat("      Found", length(idx), "peaks with high variance using `mean` method.\n")
    },
    "median" = {
      idx <- which(vars > median(vars))
      cat("      Found", length(idx), "peaks with high variance using `median` method..\n")
    },
    "q25" = {
      idx <- which(vars > quantile(vars, 0.25))
      cat("      Found", length(idx), "peaks with high variance using 25%-quantile method..\n")
    },
    "q75" = {
      idx <- which(vars > quantile(vars, 0.75))
      cat("      Found", length(idx), "peaks with high variance using 75%-quantile method..\n")
    },
    "none" = {
      # get all indicies, no filtering applied
      idx <- 1:length(vars)
      cat("      No variance filtering applied. Using all peaks.\n")
    }
  )
  return(idx)
}
