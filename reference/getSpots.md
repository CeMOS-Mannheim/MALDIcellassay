# Get the spot coordinates of spectra

Get the spot coordinates of spectra

## Usage

``` r
getSpots(object, singleSpec = TRUE)
```

## Arguments

- object:

  Object of class MALDIassay

- singleSpec:

  Logical, extract the spot coordinates of single spectra (default) or
  from average spectra.

## Value

character vector of spot coordinates. In case of average spectra
multiple spots are concatenated.

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
# spots per spectrum
getSpots(Blank2022res, singleSpec = TRUE)
#>     0     0     0     0  0.04  0.04  0.04  0.04  0.12  0.12  0.12  0.12  0.37 
#> "O13" "O14" "P13" "P14" "M13" "M14" "N13" "N14" "K13" "K14" "L13" "L14" "I13" 
#>  0.37  0.37  0.37  1.11  1.11  1.11  1.11  3.33  3.33  3.33  3.33    10    10 
#> "I14" "J13" "J14" "G13" "G14" "H13" "H14" "E13" "E14" "F13" "F14" "C13" "C14" 
#>    10    10    30    30    30    30 
#> "D13" "D14" "A13" "A14" "B13" "B14" 

#spots per concentration
getSpots(Blank2022res, singleSpec = FALSE)
#>                 0              0.04              0.12              0.37 
#> "O13 O14 P13 P14" "M13 M14 N13 N14" "K13 K14 L13 L14" "I13 I14 J13 J14" 
#>              1.11              3.33                10                30 
#> "G13 G14 H13 H14" "E13 E14 F13 F14" "C13 C14 D13 D14" "A13 A14 B13 B14" 
```
