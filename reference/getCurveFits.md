# Extract curve fits

Extract curve fits

## Usage

``` r
getCurveFits(object)
```

## Arguments

- object:

  Object of class MALDIassay

## Value

List, containing the data used to do the fits as well as the nlpr curve
fit .

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
fits <- getCurveFits(Blank2022res)
```
