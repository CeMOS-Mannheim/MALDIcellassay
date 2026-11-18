# Get mass shift for target mz

Get mass shift for target mz

## Usage

``` r
getMzShift(
  peaks,
  targetMz,
  tol,
  tolppm = FALSE,
  verbose = TRUE,
  allowNoMatch = TRUE
)
```

## Arguments

- peaks:

  List of MALDIquant::MassPeak

- targetMz:

  Numeric, target mass

- tol:

  Numeric, tolerance around targetMz

- tolppm:

  Logical, tolerance supplied in ppm

- verbose:

  Logical, print logs to the console.

- allowNoMatch:

  Logical, allow no matches.

## Value

List with two entries: `MzShift` The mass shift for each spectrum
`specIdx` The index of the spectra with a match for targetMz

## Examples

``` r
data(Blank2022peaks)
getMzShift(Blank2022peaks, targetMz = 760.585, tol = 0.1, tolppm = FALSE)
#> found mz 760.585 in 32 / 32 spectra
#> 06:12 mzshift was -0.07964664 in mean and 0.08278404  abs. max.
#> $mzshift
#>           0           0           0           0        0.04        0.04 
#> -0.08278404 -0.08278404 -0.08278404 -0.08278404 -0.08278404 -0.08278404 
#>        0.04        0.04        0.12        0.12        0.12        0.12 
#> -0.08278404 -0.08278404 -0.08278404 -0.08278404 -0.08278404 -0.08278404 
#>        0.37        0.37        0.37        0.37        1.11        1.11 
#> -0.08278404 -0.08278404 -0.08278404 -0.08278404 -0.08278404 -0.08278404 
#>        1.11        1.11          10          10          10          10 
#> -0.08278404 -0.08278404 -0.08278404 -0.08278404 -0.08278404 -0.08278404 
#>        3.33        3.33        3.33        3.33          30          30 
#> -0.08278404 -0.08278404 -0.08278404 -0.08278404 -0.05768481 -0.05768481 
#>          30          30 
#> -0.05768481 -0.05768481 
#> 
#> $specIdx
#>  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
#> [26] 26 27 28 29 30 31 32
#> 
```
