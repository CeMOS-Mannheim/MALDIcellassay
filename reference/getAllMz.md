# Get all mz value of an MALDIassay-object

Get all mz value of an MALDIassay-object

## Usage

``` r
getAllMz(object, excludeNormMz = FALSE)
```

## Arguments

- object:

  Object of class MALDIassay

- excludeNormMz:

  Logical, remove normMz from list of mz values.

## Value

numeric vector of mz values

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
head(getAllMz(Blank2022res))
#> [1] 406.9142 428.0041 438.8908 439.0423 566.9615 676.4854
```
