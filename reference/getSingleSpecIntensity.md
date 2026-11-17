# Extract the intensities of single spectra for a given mzIdx

Extract the intensities of single spectra for a given mzIdx

## Usage

``` r
getSingleSpecIntensity(object, mz_idx)
```

## Arguments

- object:

  Object of class MALDIassay

- mz_idx:

  Integer, index of mz

## Value

Numeric vector, intensities of mzIdx

## Examples

``` r

# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
head(getSingleSpecIntensity(Blank2022res, 2))
#>         0         0         0         0      0.04      0.04 
#> 0.2338765 0.1602625 0.3010994 0.1878716 0.3570634 0.1138345 
```
