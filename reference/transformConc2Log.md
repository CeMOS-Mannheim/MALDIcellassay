# Convert concentration to log10 and replace zero's

Convert concentration to log10 and replace zero's

## Usage

``` r
transformConc2Log(conc)
```

## Arguments

- conc:

  numeric, concentrations.

## Value

numeric, log10 transformed concentrations

## Examples

``` r
transformConc2Log(c(0.1, 0.01,0.001))
#> [1] -1 -2 -3
```
