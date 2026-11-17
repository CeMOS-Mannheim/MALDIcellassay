# Example MALDIassay object from Blank 2022

Object of class MALDIcellassay from MALDI mass spectrometry data of EOC
cells treated with different concentrations of SAHA. It is used to
demonstrate the usage of `MALDIcellassay`.

## Usage

``` r
data(Blank2022res)
```

## Format

Matrix with concentrations of original spectra as rownames and
m/z-values as colnames.

## Details

The concentrations include: 0, 0.04, 0.12, 0.37, 1.11, 3.33, 10 and 30
uM of SAHA at 4 replicates each. The original spectra were trimmed to
400-900 Da mass-range to keep the file size small. The peaks are the
result of
`fitCurve(spec = Blank2022spec, SinglePointRecal = TRUE, normMz = 760.585, alignTol = 0.1, normTol = 0.1)`

## References

Blank, M., Enzlein, T. & Hopf, C. LPS-induced lipid alterations in
microglia revealed by MALDI mass spectrometry-based cell fingerprinting
in neuroinflammation studies. Sci Rep 12, 2908 (2022).
https://doi.org/10.1038/s41598-022-06894-1
