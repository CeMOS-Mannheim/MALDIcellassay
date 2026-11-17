# Extract intensity using peaks as template

Extract intensity using peaks as template

## Usage

``` r
extractIntensity(mz, peaks, spec, tol)
```

## Arguments

- mz:

  numeric, mz values to be extracted from the peaks/spectra

- peaks:

  MALDIquant::MassPeaks list

- spec:

  MALDIquant::MassSpectrum list

- tol:

  numeric, tolerance in Da

## Value

MALDIquant::MassPeaks list with extracted intensities from spec at m/z
of peaks = pseudo peaks. Useful in combination with sdMassSpectrum to
get standard deviation of peaks as intensity matrix.

## Examples

``` r
data(Blank2022peaks)
data(Blank2022spec)

int <- extractIntensity(mz = c(409, 423, 440), 
                        peaks = Blank2022peaks, 
                        spec = Blank2022spec, 
                        tol = 0.2)
head(int)
#> $`0`
#> S4 class type            : MassPeaks           
#> Number of m/z values     : 3                   
#> Range of m/z values      : 409 - 440           
#> Range of intensity values: 4.91e+02 - 1.497e+03
#> Range of snr values      : 10.777 - 31.482     
#> Memory usage             : 10.719 KiB          
#> Name                     : 0.O13               
#> File                     : Y:\9-Martina Blank\R\Example_data_Rscript\20191209_EOC_plate24hR1\SAHA\0\0_O13\1\1SRef\fid
#> 
#> $`0`
#> S4 class type            : MassPeaks           
#> Number of m/z values     : 3                   
#> Range of m/z values      : 409 - 440           
#> Range of intensity values: 3.76e+02 - 1.128e+03
#> Range of snr values      : 10.359 - 29.223     
#> Memory usage             : 10.719 KiB          
#> Name                     : 0.O14               
#> File                     : Y:\9-Martina Blank\R\Example_data_Rscript\20191209_EOC_plate24hR1\SAHA\0\0_O14\1\1SRef\fid
#> 
#> $`0`
#> S4 class type            : MassPeaks           
#> Number of m/z values     : 3                   
#> Range of m/z values      : 409 - 440           
#> Range of intensity values: 5.16e+02 - 1.276e+03
#> Range of snr values      : 10.839 - 25.391     
#> Memory usage             : 10.719 KiB          
#> Name                     : 0.P13               
#> File                     : Y:\9-Martina Blank\R\Example_data_Rscript\20191209_EOC_plate24hR1\SAHA\0\0_P13\1\1SRef\fid
#> 
#> $`0`
#> S4 class type            : MassPeaks           
#> Number of m/z values     : 3                   
#> Range of m/z values      : 409 - 440           
#> Range of intensity values: 3.39e+02 - 1.025e+03
#> Range of snr values      : 9.707 - 27.484      
#> Memory usage             : 10.719 KiB          
#> Name                     : 0.P14               
#> File                     : Y:\9-Martina Blank\R\Example_data_Rscript\20191209_EOC_plate24hR1\SAHA\0\0_P14\1\1SRef\fid
#> 
#> $`0.04`
#> S4 class type            : MassPeaks           
#> Number of m/z values     : 3                   
#> Range of m/z values      : 409 - 440           
#> Range of intensity values: 5.24e+02 - 1.214e+03
#> Range of snr values      : 11.025 - 25.344     
#> Memory usage             : 10.734 KiB          
#> Name                     : 0_04.M13            
#> File                     : Y:\9-Martina Blank\R\Example_data_Rscript\20191209_EOC_plate24hR1\SAHA\0.04\0_M13\1\1SRef\fid
#> 
#> $`0.04`
#> S4 class type            : MassPeaks          
#> Number of m/z values     : 3                  
#> Range of m/z values      : 409 - 440          
#> Range of intensity values: 3.13e+02 - 7.93e+02
#> Range of snr values      : 10.578 - 24.118    
#> Memory usage             : 10.734 KiB         
#> Name                     : 0_04.M14           
#> File                     : Y:\9-Martina Blank\R\Example_data_Rscript\20191209_EOC_plate24hR1\SAHA\0.04\0_M14\1\1SRef\fid
#> 
```
