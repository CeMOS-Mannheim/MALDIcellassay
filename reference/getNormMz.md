# Extract m/z used for normalization

Extract m/z used for normalization

## Usage

``` r
getNormMz(object)
```

## Arguments

- object:

  Object of class MALDIassay

## Value

Numeric, m/z used for normalization

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
getNormMz(Blank2022res) 
#> [1] 760.585
```
