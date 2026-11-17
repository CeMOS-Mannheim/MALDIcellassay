# Extract applied normalization factors

Extract applied normalization factors

## Usage

``` r
getAppliedNormFactors(object)
```

## Arguments

- object:

  Object of class MALDIassay

## Value

Numeric vector of normalization factors applied to spectra

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
head(getAppliedNormFactors(Blank2022res))
#> [1] 16141 17827 16463 15979 11687 14565
```
