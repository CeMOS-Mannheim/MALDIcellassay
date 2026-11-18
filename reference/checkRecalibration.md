# Check the recalibration of spectra from a MALDIassay object

Dashed gray lines indicate the mz used for re-calibration ± the
tolerance. Red dashed line indicate the mz used for re-calibration and
solid lines indicate peaks. The spectrum will show the peak used for
re-calibration ± 10x the tolerance.

## Usage

``` r
checkRecalibration(object, idx, centroided = FALSE, showPeaks = TRUE)
```

## Arguments

- object:

  Object of class MALDIassay

- idx:

  Numeric, index of spectrum to plot

- centroided:

  Logical, `TRUE` for centroided data: the spectra are then drawn as
  sticks instead of lines

- showPeaks:

  Logical, mark the detected peaks with open circles

## Value

ggplot object

## Details

Profile data is drawn as lines, centroided data as sticks; the detected
peaks that are used for the recalibration are marked with open circles.

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
checkRecalibration(Blank2022res, idx = 1:8)


# centroided data is drawn as sticks
checkRecalibration(Blank2022res, idx = 1:8, centroided = TRUE)
```
