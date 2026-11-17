# load bruker MALDI target plate spectra

load bruker MALDI target plate spectra

## Usage

``` r
loadSpectra(Dir, filter = NA, nameSpectra = TRUE, verbose = TRUE)
```

## Arguments

- Dir:

  Character, parent directory of spectra.

- filter:

  Character vector, filter out spectra which match the given vector.

- nameSpectra:

  Logical, if TRUE the spectra in the resulting list will be named
  according to the dirname.

- verbose:

  Logical, print logs to the console.

## Value

List of MALDIquant::MassSpectra

## Examples

``` r
dataDir <- system.file("extdata", package="MALDIcellassay")
unzip(file.path(dataDir, "example-raw-spectra.zip"))

loadSpectra("example-raw-spectra/")
#> 10:48 Loading spectra...
#> 
#> Progress: 100%  
#> $`0.0001`
#> S4 class type            : MassSpectrum         
#> Number of m/z values     : 86835                
#> Range of m/z values      : 1802.533 - 24896.647 
#> Range of intensity values: 1.225e+05 - 8.714e+05
#> Memory usage             : 1.336 MiB            
#> Name                     : 0_0001.J1            
#> File                     : /home/runner/work/MALDIcellassay/MALDIcellassay/docs/reference/example-raw-spectra/0.0001/0_J1/1/1SLin/fid
#> 

unlink("example-raw-spectra/", recursive = TRUE)
```
