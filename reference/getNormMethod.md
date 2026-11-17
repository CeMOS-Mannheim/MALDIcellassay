# Extract normalization method

Extract normalization method

## Usage

``` r
getNormMethod(object)
```

## Arguments

- object:

  Object of class MALDIassay

## Value

Character, normalization method used.

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
getNormMethod(Blank2022res)
#> [1] "mz"
```
