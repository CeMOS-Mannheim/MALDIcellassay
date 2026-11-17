# Get the intensity matrix of single spectra for all fitted curves

Get the intensity matrix of single spectra for all fitted curves

## Usage

``` r
getIntensityMatrix(object, avg = FALSE, excludeNormMz = FALSE)
```

## Arguments

- object:

  Object of class MALDIassay

- avg:

  Logical, return single spectra intensity matrix (default) or average
  spectra intensity matrix

- excludeNormMz:

  Logical, exclude normMz from intensity matrix.

## Value

A matrix with columns as *m/z* values and rows as concentrations/spectra

## Details

Note that the returned matrix only contains *m/z* values that were
actually fitted. If a variance filtering step was applied this will not
include **all** *m/z* values. If you wish to get a matrix of **all**
*m/z* values use `MALDIquant::intensityMatrix(getSinglePeaks(object))`.
For average spectra intensity matrix with **all** *m/z* values use
`MALDIquant::intensityMatrix(getAvgPeaks(object), getAvgSpectra(object))`.

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
head(getIntensityMatrix(Blank2022res, avg = TRUE, excludeNormMz = TRUE) )
#> Warning: Could not find normMz. No exlusion possible.
#> [1] NA NA NA NA NA NA
```
