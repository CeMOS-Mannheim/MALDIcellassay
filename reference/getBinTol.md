# Get binning tolerance

Get binning tolerance

## Usage

``` r
getBinTol(object)
```

## Arguments

- object:

  Object of class MALDIassay

## Value

Numeric, tolerance used for binning in Dalton.

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
getBinTol(Blank2022res)
#> [1] 2e-04
```
