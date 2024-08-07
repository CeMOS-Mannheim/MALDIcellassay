% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qualitiyMetrics.R
\name{calculateSSMD}
\alias{calculateSSMD}
\title{Calculate strictly standardized mean difference (SSMD)}
\usage{
calculateSSMD(res, internal = TRUE, nConc = 2)
}
\arguments{
\item{res}{Object of class MALDIassay}

\item{internal}{Logical, currently only the internal implementation,
using \code{nConc} top and bottom concentrations, is implemented.}

\item{nConc}{Numeric, number of top and bottom concentrations to be used
to calculate the pseudo positive and negative control.
Only used if \code{internal} is TRUE}
}
\value{
Numeric vector of strictly standardized mean differences (SSMD)
}
\description{
Calculate strictly standardized mean difference (SSMD)
}
\details{
The strictly standardized mean difference (SSMD) is a measure of effect size.
It is the mean divided by the standard deviation of a difference between the positive and negative control.

\deqn{\gamma=\frac{\mid\mu_n - \mu_p\mid}{\sqrt{\sigma_n^2 + \sigma_p^2}}}

The SSMD can be easily be interpreted as it denotes the difference between positve and negative controls in units of standard deviation.
}
