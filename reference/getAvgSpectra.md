# Extract average spectra

Extract average spectra

## Usage

``` r
getAvgSpectra(object)
```

## Arguments

- object:

  Object of class MALDIassay

## Value

List of MALDIquantMassSpectrum

## Examples

``` r

# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
getAvgSpectra(Blank2022res)[[1]]
#> S4 class type            : MassSpectrum         
#> Number of m/z values     : 21985                
#> Range of m/z values      : 399.861 - 899.94     
#> Range of intensity values: 6.717e-05 - 1.805e+00
#> Memory usage             : 356.375 KiB          
#> Name1                    : 0.O13                
#> Name2                    : 0.O14                
#> Name3                    : 0.P13                
#> Name4                    : 0.P14                
#> File1                    : Y:\9-Martina Blank\R\Example_data_Rscript\20191209_EOC_plate24hR1\SAHA\0\0_O13\1\1SRef\fid
#> File2                    : Y:\9-Martina Blank\R\Example_data_Rscript\20191209_EOC_plate24hR1\SAHA\0\0_O14\1\1SRef\fid
#> File3                    : Y:\9-Martina Blank\R\Example_data_Rscript\20191209_EOC_plate24hR1\SAHA\0\0_P13\1\1SRef\fid
#> File4                    : Y:\9-Martina Blank\R\Example_data_Rscript\20191209_EOC_plate24hR1\SAHA\0\0_P14\1\1SRef\fid
```
