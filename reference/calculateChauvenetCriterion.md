# Calculate Chauvenet's criterion for outlier detection

Calculate Chauvenet's criterion for outlier detection

## Usage

``` r
calculateChauvenetCriterion(x)
```

## Arguments

- x:

  numeric, values (e.g. intensities) to test for outliers

## Value

logical vector, TRUE for detected outliers.

## Details

Note that, as for all outlier detection criteria: Excluding data points
from your measurement should only be conducted with extreme care. Even
if this (or any other) function tells you that a data point is an
outlier, you might still want to have it in your sample population
especially if you are not sure if your data is normal distributed. See
[Wikipedia](https://en.wikipedia.org/wiki/Chauvenet%27s_criterion) for
details of the algorithm.

## Examples

``` r
set.seed(42) 

#no outlier
sample <- rnorm(n = 8, mean = 0, sd = 0.01)
calculateChauvenetCriterion(sample)
#> [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE

# introduce outlier
sample[1] <- 1
calculateChauvenetCriterion(sample)
#> [1]  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
```
