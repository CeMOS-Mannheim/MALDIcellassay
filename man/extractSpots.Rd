% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractSpots.R
\name{extractSpots}
\alias{extractSpots}
\title{Extract the spot coordinates}
\usage{
extractSpots(spec)
}
\arguments{
\item{spec}{list of MALDIquant::MassSpectrum or MALDIquant::MassPeaks objects}
}
\value{
Character vector of spot names. If multiple spots are used (e.g. for average spectra) they will be concatenate.
}
\description{
Extract the spot coordinates
}
