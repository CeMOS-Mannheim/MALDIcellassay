# Shift mass axis

Shift mass axis

## Usage

``` r
shiftMassAxis(spec, mzdiff)
```

## Arguments

- spec:

  List of MALDIquant::MassSpectrum or MALDIquant::MassPeaks

- mzdiff:

  Numeric vector, see getMzShift()

## Value

List of MALDIquant::MassSpectrum or MALDIquant::MassPeaks with shifted
mass axis.

## Examples

``` r
data(Blank2022spec)
# raw mz
head(Blank2022spec[[1]]@mass)
#> [1] 400.0029 400.0210 400.0392 400.0574 400.0756 400.0938

# shifted mz
shifted <-shiftMassAxis(Blank2022spec[1:2], c(0.5, 0.5))
head(shifted[[1]]@mass)
#> [1] 400.5029 400.5210 400.5392 400.5574 400.5756 400.5938
```
