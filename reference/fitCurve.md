# Fit dose-response curves

Fit dose-response curves

## Usage

``` r
fitCurve(
  spec,
  unit = c("M", "mM", "uM", "nM", "pM", "fM"),
  varFilterMethod = c("mean", "median", "q25", "q75", "none"),
  monoisotopicFilter = FALSE,
  averageMethod = c("mean", "median", "sum"),
  normMz = NULL,
  normTol = 0.1,
  alignTol = 0.01,
  binTol = 2e-04,
  SNR = 3,
  halfWindowSize = 3,
  allowNoMatches = TRUE,
  normMeth = c("mz", "TIC", "median", "none"),
  SinglePointRecal = TRUE,
  verbose = TRUE
)
```

## Arguments

- spec:

  List of MALDIquant::MassSpectrum

- unit:

  Character, unit of concentration. Used to calculate the concentration
  in Moles so that pIC50 is correct. Set to "M" if you dont want changes
  in your concentrations.

- varFilterMethod:

  Character, function applied for high variance filtering. One of the
  following options `mean` (default), `median`, `q25`, `q75` or `none`
  (no filtering).

- monoisotopicFilter:

  Logical, filter peaks and just use monoisotopic peaks for curve fit.

- averageMethod:

  Character, aggregation method for average mass spectra ("mean" or
  "median")

- normMz:

  Numeric, mz used for normalization AND for single point recalibration.

- normTol:

  Numeric, tolerance in Dalton to match normMz

- alignTol:

  Numeric, tolerance for spectral alignment in Dalton.

- binTol:

  Numeric, tolerance for binning of peaks.

- SNR:

  Numeric, signal to noise ratio for peak detection.

- halfWindowSize:

  2ction. See
  [`MALDIquant::detectPeaks()`](https://rdrr.io/pkg/MALDIquant/man/detectPeaks-methods.html).

- allowNoMatches:

  Logical, if normMz can not be found in a spectrum, proceed and exclude
  spectrum or stop

- normMeth:

  Character, normalization method. Can either be "TIC", "median" or
  "mz". If "mz" then the normMz is used. If none no normalization is
  done.

- SinglePointRecal:

  Logical, perform single point recalibration to normMz

- verbose:

  Logical, print logs to console.

## Value

Object of class `MALDIassay`. The most important slot is `fits` which
contains the IC50 curve fits.

## Examples

``` r
data(Blank2022spec)

fitCurve(spec = Blank2022spec,
         SinglePointRecal = TRUE, 
         normMz = 760.585, 
         alignTol = 0.1, 
         normTol = 0.1,
         varFilterMethod = "mean") 
#> found mz 760.585 in 32 / 32 spectra
#> 06:12 mzshift was -0.07964664 in mean and 0.08278404  abs. max.
#> 06:12 normalizing... 
#> 06:12 aligning spectra... 
#> 06:12 calculating mean spectra... 
#> 06:12 Building intensity matrix and applying variance filter... 
#>       Found 216 peaks in total.
#> 06:12 fitting curves... 
#>       Found 23 peaks with high variance using `mean` method.
#> fitting 406.775 (1/23)
#> fitting 406.914 (2/23)
#> fitting 407.933 (3/23)
#> fitting 408.897 (4/23)
#> fitting 409.006 (5/23)
#> fitting 409.96 (6/23)
#> fitting 412.753 (7/23)
#> fitting 412.906 (8/23)
#> fitting 422.889 (9/23)
#> fitting 423.933 (10/23)
#> fitting 424.746 (11/23)
#> fitting 424.903 (12/23)
#> fitting 425.784 (13/23)
#> fitting 425.836 (14/23)
#> fitting 425.913 (15/23)
#> fitting 427.873 (16/23)
#> fitting 428.004 (17/23)
#> fitting 428.848 (18/23)
#> fitting 428.971 (19/23)
#> fitting 437.903 (20/23)
#> fitting 438.04 (21/23)
#> fitting 438.891 (22/23)
#> fitting 439.042 (23/23)
#> 06:12 Done! 
#> ------MALDIassay object------
#> 
#> Including 8 concentrations,
#> ranging from 0 to 30 .
#> Normalization on m/z 760.58 ± 0.1 Da.
#> 
#> Single point recalibation on 760.58 with 0.1 Da tolerance.
#> Avg. mass shift before recal.: -0.0796 Da. Max abs. shift: 0.0828 Da.
#> Avg. mass shift after recal. : -0.0127 Da. Max abs. shift: 0.0127 Da.
#> 
#> Found 216 peaks (SNR 3) and 23 high variance peaks
#> using variance filtering method: mean.
#> 
#> Top5-features based on Fold-Change and R²:
#>        mz mzIdx pEC50   R2 log2FC SSMD   V'    Z'  CRS
#> 1 832.601    22 -0.66 0.99   4.17 4.33 0.98  0.14 75.3
#> 2 810.623    19 -0.54 1.00   1.49 5.04 0.91  0.27 67.3
#> 3 811.624    20 -0.58 1.00   1.59 4.67 0.90  0.22 64.9
#> 4 804.556    18 -0.72 0.99   1.81 3.97 0.83  0.08 56.2
#> 5 835.677    23 -0.95 0.99   1.24 3.66 0.79 -0.02 40.9
```
