# Extract peaks of single spectra spectra (before average calculation)

Extract peaks of single spectra spectra (before average calculation)

## Usage

``` r
getSinglePeaks(object)
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
getSinglePeaks(Blank2022res)[[1]]
#> S4 class type            : MassPeaks            
#> Number of m/z values     : 216                  
#> Range of m/z values      : 406.775 - 898.561    
#> Range of intensity values: 4.956e-04 - 1.909e+00
#> Range of snr values      : NA - NA              
#> Memory usage             : 15.688 KiB           
#> Name                     : 0.O13                
#> File                     : Y:\9-Martina Blank\R\Example_data_Rscript\20191209_EOC_plate24hR1\SAHA\0\0_O13\1\1SRef\fid
```
