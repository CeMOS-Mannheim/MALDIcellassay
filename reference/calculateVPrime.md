# Calculate V'-Factor

Calculate V'-Factor

## Usage

``` r
calculateVPrime(res)
```

## Arguments

- res:

  Object of class MALDIassay

## Value

Numeric vector of V'-factors

## Details

The V'-factor is a generalization of the Z'-factor to a dose-response
curve. See [M.-A. Bray and A. Carpenter, Advanced assay development
guidelines for image-based high content screening and
analysis](https://www.ncbi.nlm.nih.gov/books/NBK126174/pdf/Bookshelf_NBK126174.pdf)
for details. It is defined as: \$\$V' = 1 - 6 \* \sigma_f/\|\mu_p -
\mu_n\|\$\$

with

\$\$\sigma_f = \sqrt{1/N \* \sum{y_fit - y_measured}^2}\$\$

In other words, \\\sigma_f\\ is the standard deviation of residuals.

Note, we do not need to estimate the variance for the mean of the
positive and negative value. So, this function uses the top and bottom
asymptote directly instead of taking the top and bottom concentrations
in consideration.

## Examples

``` r
# see example for `fitCurve()` to see how this data was generated
data(Blank2022res)

calculateVPrime(Blank2022res) 
#> 406.914163450922 428.004132260612 438.890824557592 439.042346897966 
#>      -0.16705820       0.59861762       0.99681015       0.09886635 
#> 566.961527775518 676.485447847528 703.586286847746 725.549662904297 
#>       0.32258091       0.52759605       0.54818121       0.42672045 
#> 726.555208487752 727.552245980755 741.536315662242 756.548750613787 
#>       0.51068609       0.76295188       0.32371110      -0.11915923 
#> 772.540528346963  782.57644400297 783.582654272217 798.563919325163 
#>       0.25670312       0.90600765       0.90835642       0.59311636 
#> 799.567517217833 804.556260133438 810.622915627018 811.624361313474 
#>       0.51790518       0.83299718       0.91141686       0.90021277 
#> 824.581961130221 832.600992026117  835.67738341154 
#>       0.33031185       0.97768687       0.78933321 
```
