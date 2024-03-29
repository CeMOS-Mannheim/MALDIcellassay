.getSigmaFit <- function(res) {
  fits <- getCurveFits(res)

  purrr::map_dbl(fits, function(x) {
    yfit =  getFitValues(x$model)
    y =  getY(x$model)
    n <- length(y)

    sigmaFit <- sqrt(1/n*sum((y - yfit)^2))
    return(sigmaFit)
  })
}

.calculateZPrime <- function(pos, neg) {
  if(length(pos) < 2 | length(neg) < 2) {
    return(NA_integer_)
  }

  meanPos <- mean(pos, na.rm = TRUE)
  meanNeg <- mean(neg, na.rm = TRUE)
  sdPos <- sd(pos, na.rm = TRUE)
  sdNeg <- sd(neg, na.rm = TRUE)

  z <- 1 - (3 * (sdPos + sdNeg) / (abs(meanPos - meanNeg)))

  return(z)
}

.calculateSSMD <- function(pos, neg) {
  # strictly standardized mean difference
  if(length(pos) < 2 | length(neg) < 2) {
    return(NA_integer_)
  }

  meanPos <- mean(pos, na.rm = TRUE)
  meanNeg <- mean(neg, na.rm = TRUE)
  sdPos <- sd(pos, na.rm = TRUE)
  sdNeg <- sd(neg, na.rm = TRUE)

  ssmd <- abs(meanPos - meanNeg) /  sqrt(sdNeg^2 + sdPos^2)

  return(ssmd)
}

.getTopAndBottomIntensityMatrix <- function(res, nConc) {
  concs <- unique(getConc(res))

  topIdx <- which(getConc(res) %in% concs[1:nConc])
  botIdx <- which(getConc(res) %in% concs[(length(concs)-nConc+1):length(concs)])

  intmatPos <- getIntensityMatrix(res,
                                  avg = FALSE,
                                  excludeNormMz = FALSE)[botIdx,]
  intmatNeg <- getIntensityMatrix(res,
                                  avg = FALSE,
                                  excludeNormMz = FALSE)[topIdx,]

  return(list(pos = intmatPos,
              neg = intmatNeg))
}

.calculateMetric <- function(res, nConc = 2, fun) {

  int <- .getTopAndBottomIntensityMatrix(res = res,
                                         nConc = nConc)


  metric <- vapply(X = seq_along(getAllMz(res,
                                          excludeNormMz = FALSE)),
                   FUN = function(i) {
                     fun(pos = int$pos[, i],
                         neg = int$neg[, i])
                   },
                   FUN.VALUE = numeric(1))

  return(metric)
}

#' Calculate Z'-factor of assay quality
#'
#' @param res      Object of class MALDIassay
#' @param internal Logical, currently only the internal implementation,
#'                 using `nConc` top and bottom concentrations, is implemented.
#' @param nConc    Numeric, number of top and bottom concentrations to be used
#'                 to calculate the pseudo positive and negative control.
#'                 Only used if `internal` is TRUE
#'
#' @details
#' The most common way to measure the quality of an assay is the so-called Z'-factor,
#' which describes the separation of the positive and negative control in terms of their standard deviations \eqn{\sigma_p} and \eqn{\sigma_n}.
#' The Z'-factor is defined as [Ji-Hu Zhang et al., A simple statistical parameter for use in evaluation and validation of high throughput screening assays](https://journals.sagepub.com/doi/abs/10.1177/108705719900400206).
#' \deqn{Z' = 1 - (3 * (\sigma_p+\sigma_n))/|\mu_p-\mu_n|}
#'
#' where \eqn{\mu_p} and \eqn{\mu_p} is the mean value of the positive (response expected) and negative (no response expected) control, respectively.
#' Therefore, the assay quality is **independent of the shape of the concentration response curve** and solely depend on two control values.
#'
#' Note, if `internal` is set to TRUE, the `nConc` highest concentrations are assumed as positive control,
#' whereas the `nConc` lowest concentrations are used as negative.
#'
#' |**Value** |**Interpretation** |
#' |--|--|
#' |Z' ~ 1 | perfect assay |
#' |1 > Z' > 0.5 | excellent assay |
#' |0.5 > Z' > 0 | moderate assay |
#' |Z' = 0 | good only for yes/no response |
#' |Z' < 0 | unacceptable  |#'
#'
#' @return
#' Numeric vector of Z'-factors.
#' @importFrom purrr map_dbl
#' @export
calculateZPrime <- function(res, internal = TRUE, nConc = 2) {
  if(!internal) {
    cat("Currently only the internal implementation,
        using nConc top and bottom concentrations, is implemented.\n")
    return()
  }

  .calculateMetric(res = res,
                   nConc = nConc,
                   fun = .calculateZPrime)
}


#' Calculate strictly standardized mean difference (SSMD)
#'
#' @param res      Object of class MALDIassay
#' @param internal Logical, currently only the internal implementation,
#'                 using `nConc` top and bottom concentrations, is implemented.
#' @param nConc    Numeric, number of top and bottom concentrations to be used
#'                 to calculate the pseudo positive and negative control.
#'                 Only used if `internal` is TRUE
#'
#' @details
#' The strictly standardized mean difference (SSMD) is a measure of effect size.
#' It is the mean divided by the standard deviation of a difference between the positve and negative control.
#'
#' \deqn{γ=|\mu_n - \muμ_p|/\sqrt{\sigma_n^2 + \sigma_p^2}}
#'
#' The SSMD can be easily be interpreted as it denotes the difference between positve and negative controls in units of standard deviation.
#'
#'
#'
#' @return
#' Numeric vector of strictly standardized mean differences (SSMD)
#' @export
calculateSSMD <- function(res, internal = TRUE, nConc = 2) {
  if(!internal) {
    cat("Currently only the internal implementation,
        using nConc top and bottom concentrations, is implemented.\n")
    return()
  }

  .calculateMetric(res = res,
                   nConc = nConc,
                   fun = .calculateSSMD)
}

#'  Calculate V'-Factor
#'
#' @param res      Object of class MALDIassay
#' @param internal Logical, currently only the internal implementation,
#'                 using `nConc` top and bottom concentrations, is implemented.
#'
#' @details
#' The V'-factor is a generalization of the Z'-factor to a dose-response curve.
#' See [M.-A. Bray and A. Carpenter, Advanced assay development guidelines for image-based high content screening and analysis](https://www.ncbi.nlm.nih.gov/books/NBK126174/pdf/Bookshelf_NBK126174.pdf) for details.
#' It is defined as:
#' \deqn{V' = 1 - 6 * \sigma_f/|\mu_p - \mu_n|}
#'
#' with
#'
#' \deqn{\sigma_f = \sqrt{1/N * \sum{y_fit - y_measured}^2}}
#'
#' In other words, \eqn{\sigma_f} is the standard deviation of residuals.
#'
#' Note, we do not need to estimate the variance for the mean of the positive and negative value.
#' So, this function uses the top and bottom asymptote directly instead of taking the top and bottom concentrations in consideration.
#'
#' @return
#' Numeric vector of V'-factors
#' @export
calculateVPrime <- function(res, internal = TRUE) {
  if(!internal) {
    cat("Currently only the internal implementation,
        using nConc top and bottom concentrations, is implemented.\n")
    return()
  }

  param <- getFittingParameters(res)

  meanPos <- as.numeric(param$top)
  meanNeg <- as.numeric(param$bottom)

  sigmaFit <- .getSigmaFit(res)

  v <- 1 - 6 * (sigmaFit/(abs(meanNeg - meanPos)))

  return(v)
}

.integrateMSWD <- function(chi2, nu) {
  result <- integrate(
    function(x, nu) {
      (x^(nu/2 - 1) * exp(-x/2)) / (2^(nu/2) * gamma(nu/2))
    },
    lower = chi2,
    upper = Inf,
    nu = nu,
    stop.on.error = FALSE)

  return(result$value)
}

#' Calculate reduced Chi-squared statistic / mean squared weighted deviation (MSWD)
#'
#' @param res      Object of class MALDIassay
#'
#' @details
#' The reduced chi-square statistic is used extensively in goodness of fit testing.
#' It is also known as mean squared weighted deviation (**MSWD**) in isotopic dating
#' and variance of unit weight in the context of weighted least squares.
#'
#' \deqn{\chi^2_\nu =\frac{\chi^2}{\nu}}
#'
#' \deqn{\chi^2=\sum_{N}{\frac{1}{\sigma^2_i}*[y_i - f(x_i)]^2}}
#'
#' \deqn{\nu = N - m - 1}
#' @return
#' Numeric vector of MSWD
#' @export
#'
calculateMSWD <- function(res) {
  fits <- getCurveFits(res)

  intmat <- intensityMatrix(getSinglePeaks(res))

  pval <- purrr::map_dbl(seq_along(fits),
                 function(i) {
                   fit <- fits[[i]]$model
                   int <- intmat[,i]
                   conc <- getConc(res)

                   vars <- vapply(1:length(unique(conc)),
                                  function(j) {
                                    var(int[conc==unique(conc)[j]])
                                  },
                                  FUN.VALUE = numeric(1))

                   if(all(vars==0)) {
                     return(NA)
                   }

                   yfit <-  getFitValues(fit)
                   yobs <- getY(fit)

                   chisq <- sum(1 / vars * (yobs - yfit)^2)

                   redChisq <- chisq/(length(yfit) - getPar(fit)$npar - 1)

                   q <- .integrateMSWD(chi2 = chisq, nu = (length(yfit) - getPar(fit)$npar - 1))

                   return(q)
                 })

  return(q)
}
