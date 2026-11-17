# Extract tolerance used for normalization

Extract tolerance used for normalization

## Usage

``` r
getNormMzTol(object)
```

## Arguments

- object:

  Object of class MALDIassay

## Value

Numeric, tolerance used for normalization

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
getNormMzTol(Blank2022res)
#> [1] 0.1
```
