# Calculate Z'-factor of assay quality

Calculate Z'-factor of assay quality

## Usage

``` r
calculateZPrime(res, nConc = 2)
```

## Arguments

- res:

  Object of class MALDIassay

- nConc:

  Numeric, number of top and bottom concentrations to be used to
  calculate the pseudo positive and negative control.

## Value

Numeric vector of Z'-factors.

## Details

The most common way to measure the quality of an assay is the so-called
Z'-factor, which describes the separation of the positive and negative
control in terms of their standard deviations \\\sigma_p\\ and
\\\sigma_n\\. The Z'-factor is defined as [Ji-Hu Zhang et al., A simple
statistical parameter for use in evaluation and validation of high
throughput screening assays](https://pubmed.ncbi.nlm.nih.gov/10838414/).
\$\$Z' = 1 - (3 \* (\sigma_p+\sigma_n))/\|\mu_p-\mu_n\|\$\$

where \\\mu_p\\ and \\\mu_p\\ is the mean value of the positive
(response expected) and negative (no response expected) control,
respectively. Therefore, the assay quality is **independent of the shape
of the concentration response curve** and solely depend on two control
values.

Note, the `nConc` highest concentrations are assumed as positive
control, whereas the `nConc` lowest concentrations are used as negative.

|                |                               |
|----------------|-------------------------------|
| **Value**      | **Interpretation**            |
| Z' ~ 1         | perfect assay                 |
| 1 \> Z' \> 0.5 | excellent assay               |
| 0.5 \> Z' \> 0 | moderate assay                |
| Z' = 0         | good only for yes/no response |
| Z' \< 0        | unacceptable                  |

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)
calculateZPrime(Blank2022res, nConc = 2)       
#>  [1] -3.237607219 -1.252868716 -5.541275313 -4.970064620 -1.631465147
#>  [6] -0.957522437 -0.223832609 -0.902466878 -1.017670816 -1.058958794
#> [11] -1.562843832 -2.115946846 -0.024041256 -0.812729462 -0.691440743
#> [16] -0.009526589  0.030512012  0.077456813  0.266325098  0.216154352
#> [21]  0.151846342  0.140734203 -0.020104856
 
```
