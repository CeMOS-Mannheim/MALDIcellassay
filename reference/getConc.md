# Extract the concentrations used in a MALDIassay

Extract the concentrations used in a MALDIassay

## Usage

``` r
getConc(object)
```

## Arguments

- object:

  Object of class MALDIassay

## Value

Numeric vector, concentrations used in a MALDIassay

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
head(getConc(Blank2022res))
#> [1] 0.00 0.00 0.00 0.00 0.04 0.04
```
