# Extract applied mz-shift

Extract applied mz-shift

## Usage

``` r
getAppliedMzShift(object)
```

## Arguments

- object:

  Object of class MALDIassay

## Value

Numeric vector of mz-shits applied to spectra

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
head(getAppliedMzShift(Blank2022res))
#> [1] -0.0828 -0.0828 -0.0828 -0.0828 -0.0828 -0.0828
```
