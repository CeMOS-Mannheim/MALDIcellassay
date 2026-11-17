# Extract the spot coordinates

Extract the spot coordinates

## Usage

``` r
extractSpots(spec)
```

## Arguments

- spec:

  list of MALDIquant::MassSpectrum or MALDIquant::MassPeaks objects

## Value

Character vector of spot names. If multiple spots are used (e.g. for
average spectra) they will be concatenate.

## Examples

``` r
data(Blank2022spec)
head(extractSpots(Blank2022spec))
#>     0     0     0     0  0.04  0.04 
#> "O13" "O14" "P13" "P14" "M13" "M14" 
```
