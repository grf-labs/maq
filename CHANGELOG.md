# Changelog
All notable changes to `maq` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2023-06-27
First CRAN beta release. The R package currently supports

- Fitting Qini curves for an arbitrary number of arms using `maq(cate.hat, cost.hat, max.budget, Y.eval, ...)`, where `cate.hat` and `cost.hat` are CATE and cost estimates obtained via some function learned on a training set, `Y.eval` are evaluation scores on a test set (for example inverse-propensity weighted outcomes) and `max.budget` is the maximum mean budget/unit to fit the curve on. Setting the option `target.with.covariates` to `FALSE` yields a baseline Qini curve that can be used to assess the value of treatment targeting based on covariates.
- Computing point estimates and standard errors for any point on the curve with `average_gain()` (supporting clustered standard errors if fit with clusters).
- Comparing arbitrary points on different curves using `difference_gain()`, yielding standard errors that account for the correlation arising from curves fit on the same evaluation data.
- Retrieving the underlying "induced" policy at arbitrary spend levels with `predict()`.
- Basic plotting functionality with `plot()`.
- Retrieving the full gain/spend/allocation path with `summary(maq.object)`.

The Python package (source install only) currently supports

- Fitting the Qini curve, and retrieving standard errors for any point on the curve, as well as the implied policy \pi (remaining functionality can be added on the wrapper side - the heavy lifting is done in core/C++, and the input/output is language agnostic).
