# Extract SNR used for peak detection

Extract SNR used for peak detection

## Usage

``` r
getSNR(object)
```

## Arguments

- object:

  Object of class MALDIassay

## Value

Numeric, SNR used for peak detection

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
getSNR(Blank2022res)
#> [1] 3
```
