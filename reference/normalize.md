# Normalize spectra and peaks

Normalize spectra and peaks

## Usage

``` r
normalize(spec, peaks, normMeth, normMz, normTol)
```

## Arguments

- spec:

  List of MALDIquant::MassSpectrum or MALDIquant::MassPeaks objects

- peaks:

  List of MALDIquant::MassPeaks

- normMeth:

  Character, normalization method. Options are "TIC", "median" and "mz".

- normMz:

  Numeric, mz used to normalize.

- normTol:

  Numeric, tolerance around normMz.

## Value

List of lists of normalized MALDIquant::MassSpectrum, normalized
MALDIquant::MassPeaks, normalization factors as well as indicies of
spectra containing the `normMz` in case of `normMeth = "mz"`,

## Details

Typically MassSpectrum lists should be provided as `spec`, still its
possible to provide peaks for both arguments `spec` and `peaks` in the
case that no raw data (continuous spectra) is available. In this case
the operations performed on `spec` and `peaks` are equal.

## Examples

``` r
data(Blank2022spec)
data(Blank2022peaks)
norm <- normalize(Blank2022spec, Blank2022peaks, normMeth = "mz", normMz = 760.585, normTol = 0.1)

# normalization factors
norm$factor 
#>     0     0     0     0  0.04  0.04  0.04  0.04  0.12  0.12  0.12  0.12  0.37 
#> 16141 17827 16463 15979 11687 14565 15224 11857  7949 13138 19570 13271 10036 
#>  0.37  0.37  0.37  1.11  1.11  1.11  1.11    10    10    10    10  3.33  3.33 
#> 12771 16056 16157 17366 14420 15080 18108 10366 17430 12866 14851 12895 14920 
#>  3.33  3.33    30    30    30    30 
#> 16715 13496  8735 14091 13168  7955 
```
