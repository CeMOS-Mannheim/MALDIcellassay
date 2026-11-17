# Example spectra from Blank 2022

MALDI mass spectrometry data of EOC cells treated with different
concentrations of SAHA. It is used to demonstrate the usage of
`MALDIcellassay`.

## Usage

``` r
data(Blank2022spec)
```

## Format

A list of MALDIquant::MassSpectrum-objects named with the respective
concentration.

## Details

The concentrations include: 0, 0.04, 0.12, 0.37, 1.11, 3.33, 10 and 30
uM of SAHA at 4 replicates each. The original spectra were trimmed to
400-900 Da mass-range to keep the file size small.

## References

Blank, M., Enzlein, T. & Hopf, C. LPS-induced lipid alterations in
microglia revealed by MALDI mass spectrometry-based cell fingerprinting
in neuroinflammation studies. Sci Rep 12, 2908 (2022).
https://doi.org/10.1038/s41598-022-06894-1
