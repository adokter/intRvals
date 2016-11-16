# intRval
This package can be used to analyse interval data with missed arrival observations. I developed it for estimating goose faecal output rates (dropping rates).

When measuring intervals of dropping events in the field, observers regularly fail to see a dropping, leading to longer intervals (at integer multiples of the true dropping interal). Inter-arrival times are fitted to a distribution that accounts for these missed observations.

The package corrects mean and variance of the rate for the effects of missed observations, and provides simple summary statistics and tests for comparing means and variances (analogous to R's native ``t.test`` and ``var.test`` functions)

###Installation in R
```
library(devtools)
install_github("adokter/intRval")
```
