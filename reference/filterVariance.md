# Filter for high variance signals

Filter for high variance signals

## Usage

``` r
filterVariance(
  vars,
  method = c("mean", "median", "q25", "q75", "none"),
  verbose = TRUE
)
```

## Arguments

- vars:

  Numeric vector, variances of signals

- method:

  Character, filtering method. One of "mean" (default), "median", "q25",
  "q75" (25 and 75% quantile) or "none".

- verbose:

  Logical, print logs to console.

## Value

Indices of spectra with a high variance

## Examples

``` r
data(Blank2022intmat)
# get variance of each peak
vars <- apply(Blank2022intmat, 2, var)
highVarIndicies <- filterVariance(vars, method = "mean", verbose = TRUE)
#>       Found 57 peaks with high variance using `mean` method.
```
