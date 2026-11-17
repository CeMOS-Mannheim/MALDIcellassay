# Extract peak statistics

Extract peak statistics

## Usage

``` r
getPeakStatistics(object, summarise = FALSE)
```

## Arguments

- object:

  Object of class MALDIassay

- summarise:

  Logical, return summarized results (one result per mz and not per mz
  and spectra)

## Value

A tibble with peak statistics (R², fold-change, CV%, etc.)

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
head(getPeakStatistics(Blank2022res, summarise = TRUE))
#> # A tibble: 6 × 9
#>      mz mzIdx pEC50    R2 log2FC  SSMD  `V'`  `Z'`   CRS
#>   <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1  407.     1  0.28  0.76  -1.17  0.99 -0.17 -3.24     0
#> 2  428.     2  0.22  0.98  -1.84  1.83  0.6  -1.25     0
#> 3  439.     3 -1.07  1      8.81  0.46  1    -5.54     0
#> 4  439.     4 -0.64  0.84   0.73  0.57  0.1  -4.97     0
#> 5  567.     5 -0.7   0.92   1.53  1.29  0.32 -1.63     0
#> 6  676.     6 -0.99  0.93   0.92  1.79  0.53 -0.96     0
```
