# Get the mz value associated with a mzIdx

Get the mz value associated with a mzIdx

## Usage

``` r
getMzFromMzIdx(object, mzIdx)
```

## Arguments

- object:

  Object of class MALDIassay

- mzIdx:

  numeric, index of mass of interest (see
  [`getPeakStatistics()`](https://cemos-mannheim.github.io/MALDIcellassay/reference/getPeakStatistics.md))

## Value

numeric, mz value

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
getMzFromMzIdx(Blank2022res, mzIdx = 2)
#> [1] 428.004
```
