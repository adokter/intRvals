# intRvals
Package intRvals calculates means and variances of arrival intervals (and arrival rates) corrected for missed arrival observations, and compares means and variances of groups of interval data.

### Installation in R
```
library(devtools)
install_github("adokter/intRvals")
```

### General
The central function of package `intRvals` is `estinterval`, which is used to estimate the
mean arrival interval (and its standard deviation) from interval data with missed arrivals. This is
achieved by fitting the theoretical probability density `intervalpdf` to the interval data

The package can be used to analyse general interval data where
intervals are derived from distinct arrival observations.
For example, the authors have used it to analyze dropping intervals
of grazing geese for estimating their faecal output.

Intervals are defined as the time between observed arrival events (e.g. the time between one excreted droppings to the next)
The package provides a way of taking into account missed observations
(e.g. defecations), which lead to occasional observed intervals at integer multiples of the
true arrival interval.

### Typical workflow
* Fit interval model `m` to an interval dataset `d` using `estinterval`, as in `m=estinterval(d)`.
* Visually inspect model fits using `plot.intRvals`, as in `plot(m)`.
* Use `anova.intRvals` to check whether the missed event probability was signficantly different from zero, as in `anova(m)`
* Also use `anova.intRvals` to perform model selection between competing models `m1`,`m2` for the same interval dataset `d`, as in `anova(m1,m2)`
* Compare means and variances between different interval datasets `d1`,`d2` using `ttest` and `vartest`

### Other useful functionality
* `fold` provides functionality to fold observed intervals back to their fundamental interval
* `fundamental` tests which intervals are fundamental, i.e. intervals not containing a missed arrival observation
* `interval2rate` converts interval estimates to rates
* `partition` estimates and tests for the presence of within-subject variation
* `intervalsim` simulates a set of observed intervals

The package comes with a example interval dataset `goosedrop`

### References
* Dokter, A.M., et al. 2017. Analysing time-ordered event data with missed observations, Ecology and Evolution, 2017, in press.
* Bédard, J. & Gauthier, G. 1986. Assessment of faecal output in geese. Journal of Applied Ecology, 23, 77-90.
* Owen, M. 1971. The Selection of Feeding Site by White-Fronted Geese in Winter. Journal of Applied Ecology 8: 905-917.
