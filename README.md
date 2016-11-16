# intRval
IntRval calculates means and variances of arrival intervals corrected for missed arrival observations, and compares means and variances of groups of interval data.

The package was designed originally for analysing dropping intervals of grazing geese to estimate their faecal output, but can be used to analyse general interval data where intervals are derived from distinct arrival observations.

Intervals are defined as the time between observed arrival events (e.g. the time between one excreted droppings to the next) The package provides a way of taking into account missed observations (excreted droppings), which lead to occasional observed intervals at integer multiples of the true arrival interval.

When observing intervals of arrival events, observers may fail to observe arrival events (e.g. a dropping excretion in the case of geese). With this package inter-arrival times can be fitted to a distribution that accounts for these missed observations.

The package corrects mean and variance of the rate for the effects of missed observations, and provides simple summary statistics and tests for comparing means and variances (analogous to R's native ``t.test`` and ``var.test`` functions)

###Installation in R
```
library(devtools)
install_github("adokter/intRval")
```
