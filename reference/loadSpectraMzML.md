# load mzML spectra

load mzML spectra

## Usage

``` r
loadSpectraMzML(
  Dir,
  filter = NA,
  nameSpectra = TRUE,
  verbose = TRUE,
  centroided = FALSE
)
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

  Logical, print logs to console

- centroided:

  Logical, if FALSE (default) a list of MALDIquant::MassSpectrum is
  returned, otherwise a list of MALDIquant::MassPeaks are returned.

## Value

           List of MALDIquant::MassSpectra

## Examples

``` r
dataDir <- system.file("extdata", package="MALDIcellassay")

loadSpectraMzML(file.path(dataDir, "Koch2024mzML"))
#> 06:12 Loading spectra...
#> 
#> Progress: 100%  
#> $`0.0001`
#> S4 class type            : MassSpectrum         
#> Number of m/z values     : 86835                
#> Range of m/z values      : 1802.533 - 24896.647 
#> Range of intensity values: 1.225e+05 - 8.714e+05
#> Memory usage             : 1.327 MiB            
#> File                     : /home/runner/work/_temp/Library/MALDIcellassay/extdata/Koch2024mzML/0.0001.mzML
#> 
```
