# Extract peaks of average spectra

Extract peaks of average spectra

## Usage

``` r
getAvgPeaks(object)
```

## Arguments

- object:

  Object of class MALDIassay

## Value

List of MALDIquantMassPeaks

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
getAvgPeaks(Blank2022res)[[1]]
#> S4 class type            : MassPeaks            
#> Number of m/z values     : 169                  
#> Range of m/z values      : 406.914 - 896.58     
#> Range of intensity values: 3.853e-03 - 1.805e+00
#> Range of snr values      : 3 - 268.316          
#> Memory usage             : 16.977 KiB           
#> Name1                    : 0.O13                
#> Name2                    : 0.O14                
#> Name3                    : 0.P13                
#> Name4                    : 0.P14                
#> File1                    : Y:\9-Martina Blank\R\Example_data_Rscript\20191209_EOC_plate24hR1\SAHA\0\0_O13\1\1SRef\fid
#> File2                    : Y:\9-Martina Blank\R\Example_data_Rscript\20191209_EOC_plate24hR1\SAHA\0\0_O14\1\1SRef\fid
#> File3                    : Y:\9-Martina Blank\R\Example_data_Rscript\20191209_EOC_plate24hR1\SAHA\0\0_P13\1\1SRef\fid
#> File4                    : Y:\9-Martina Blank\R\Example_data_Rscript\20191209_EOC_plate24hR1\SAHA\0\0_P14\1\1SRef\fid
```
