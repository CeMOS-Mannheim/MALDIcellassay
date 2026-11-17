# Extract variance filtering method

Extract variance filtering method

## Usage

``` r
getVarFilterMethod(object)
```

## Arguments

- object:

  Object of class MALDIassay

## Value

Character of variance filtering method used

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
getVarFilterMethod(Blank2022res)
#> [1] "mean"
```
