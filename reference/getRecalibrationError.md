# Calculate remaining calibration error of a MALDIassay object

Calculate remaining calibration error of a MALDIassay object

## Usage

``` r
getRecalibrationError(object)
```

## Arguments

- object:

  Object of class MALDIassay

## Value

A tibble containing statistics about remaining calibration error

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
getRecalibrationError(Blank2022res)
#> # A tibble: 1 × 5
#>   meanAbs sdAbs maxAbs    mean    sd
#>     <dbl> <dbl>  <dbl>   <dbl> <dbl>
#> 1  0.0126     0 0.0126 -0.0126     0
```
