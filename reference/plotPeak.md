# Plot a peak of interest from a MALDIassay object

Draws the region around the peak of interest using the averaged spectra.
For profile data the spectra are shown as curves, for centroided data
the centroids are shown as sticks.

## Usage

``` r
plotPeak(object, mzIdx, tol = 0.8, centroided = FALSE, showPeaks = TRUE)
```

## Arguments

- object:

  object of class MALDIassay

- mzIdx:

  numeric, index of mass of interest (see
  [`getPeakStatistics()`](https://cemos-mannheim.github.io/MALDIcellassay/reference/getPeakStatistics.md))

- tol:

  numeric, tolerance around peak to plot

- centroided:

  logical, whether the data are centroided. If `TRUE` the centroids are
  drawn as sticks instead of profile curves.

- showPeaks:

  logical, whether the peaks used for curve fitting are marked with open
  circles (default `TRUE`).

## Value

ggplot object

## Details

With `showPeaks = TRUE` the peaks that are used for curve fitting are
marked with open circles: the peaks detected on the averaged spectra for
profile data, and the centroids themselves for centroided data.

## Examples

``` r

data(Blank2022res)
plotPeak(Blank2022res, mzIdx = 2)


# centroided data: centroids are drawn as sticks
plotPeak(Blank2022res, mzIdx = 2, centroided = TRUE)
```
