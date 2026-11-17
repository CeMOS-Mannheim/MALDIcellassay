# Get fitting parameters

Get fitting parameters

## Usage

``` r
getFittingParameters(object)
```

## Arguments

- object:

  Object of class MALDIassay

## Value

tibble of fitting parameters for each fitted m/z-value

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
head(getFittingParameters(Blank2022res))
#> # A tibble: 6 × 8
#>   mzIdx mz                npar  bottom   top   xmid   scal     s
#>   <dbl> <chr>            <dbl>   <dbl> <dbl>  <dbl>  <dbl> <dbl>
#> 1     1 406.914163450922     4 0.112   0.252 -0.251 -0.772     1
#> 2     2 428.004132260612     4 0.0595  0.212 -0.217 -5.57      1
#> 3     3 438.890824557592     4 0.00130 0.584  1.08   6.22      1
#> 4     4 439.042346897966     4 0.652   1.08   0.686  1.39      1
#> 5     5 566.961527775518     4 0.0706  0.203  0.694  2.37      1
#> 6     6 676.485447847528     4 0.551   1.04   1.04   1.95      1
```
