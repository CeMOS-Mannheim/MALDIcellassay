# Get normalization factors from peak data.frame

Get normalization factors from peak data.frame

## Usage

``` r
getNormFactors(peaks, targetMz, tol, tolppm = TRUE, allowNoMatch = TRUE)
```

## Arguments

- peaks:

  List of Object of class MALDIquant::MassPeaks

- targetMz:

  Numeric, target mass

- tol:

  Numeric, tolerance around targetMz

- tolppm:

  Logical, is the tolerance provided in ppm (TRUE) or Daltion (FALSE)

- allowNoMatch:

  Logical, stop if targetMz is not fround in single spectrum? If TRUE
  spectra without targetMz match will be excluded.

## Value

             List with two entries:
                                      norm_factor The normalization factor for each spectrum
                                      specIdx     The index of the spectra with a match for targetMz

## Examples

``` r
data(Blank2022peaks)
getNormFactors(Blank2022peaks, targetMz = 760.585, tol = 0.1, tolppm = FALSE)
#> $norm_factor
#>     0     0     0     0  0.04  0.04  0.04  0.04  0.12  0.12  0.12  0.12  0.37 
#> 16141 17827 16463 15979 11687 14565 15224 11857  7949 13138 19570 13271 10036 
#>  0.37  0.37  0.37  1.11  1.11  1.11  1.11    10    10    10    10  3.33  3.33 
#> 12771 16056 16157 17366 14420 15080 18108 10366 17430 12866 14851 12895 14920 
#>  3.33  3.33    30    30    30    30 
#> 16715 13496  8735 14091 13168  7955 
#> 
#> $specIdx
#>  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
#> [26] 26 27 28 29 30 31 32
#> 
```
