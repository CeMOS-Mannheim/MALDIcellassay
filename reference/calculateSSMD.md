# Calculate strictly standardized mean difference (SSMD)

Calculate strictly standardized mean difference (SSMD)

## Usage

``` r
calculateSSMD(res, nConc = 2)
```

## Arguments

- res:

  Object of class MALDIassay

- nConc:

  Numeric, number of top and bottom concentrations to be used to
  calculate the pseudo positive and negative control.

## Value

Numeric vector of strictly standardized mean differences (SSMD)

## Details

The strictly standardized mean difference (SSMD) is a measure of effect
size. It is the mean divided by the standard deviation of a difference
between the positive and negative control.

\$\$\gamma=\frac{\mid\mu_n - \mu_p\mid}{\sqrt{\sigma_n^2 +
\sigma_p^2}}\$\$

The SSMD can be easily be interpreted as it denotes the difference
between positve and negative controls in units of standard deviation.

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)

calculateSSMD(Blank2022res, nConc = 2)       
#>  [1] 0.9929269 1.8295971 0.4589847 0.5734228 1.2892462 1.7888130 3.1918283
#>  [8] 1.8648946 1.6817642 1.6596281 1.6481736 1.2971227 4.0481543 2.3039603
#> [15] 2.4210874 3.8748702 3.8690381 3.9669764 5.0425580 4.6732764 4.3425562
#> [22] 4.3308703 3.6642894
```
