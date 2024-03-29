% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getters.R
\name{getSpots}
\alias{getSpots}
\title{Get the spot coordinates of spectra}
\usage{
getSpots(object, singleSpec = TRUE)
}
\arguments{
\item{object}{Object of class MALDIassay}

\item{singleSpec}{Logical, extract the spot coordinates of single spectra (default) or from average spectra.}
}
\value{
character vector of spot coordinates. In case of average spectra multiple spots are concatenated.
}
\description{
Get the spot coordinates of spectra
}
