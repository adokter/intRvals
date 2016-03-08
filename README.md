# droprate
This package can be used to analyse interval data with missed event observations. I developed it for estimating goose faecal output rates (dropping rates). When measuring intervals of dropping events in the field, it often happens that observers fail to see a dropping, leading to longer intervals (at integer multiples of the true dropping interal). Event intervals are fitted to a distribution that accounts for these missed observations. The package corrects mean and variance of the rate for the effects of missed observations, and provides simple summary statistics and tests for comparing means and variances (analogous to R's ``t.test`` and ``var.test`` default functions)

###Installation in R
```
library(devtools)
install_github("adokter/droprate")
```
