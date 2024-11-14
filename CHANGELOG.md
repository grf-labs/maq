# Changelog
All notable changes to `maq` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.5.0] - 2024-11-??

### Added
- Invisibly return data when calling plot to allow for optional customization with other plotting libraries. [#94](https://github.com/grf-labs/maq/pull/94)

### Fixed
- Omit empty CI plots lines when R = 0. [#76](https://github.com/grf-labs/maq/pull/76)

## [0.4.0] - 2024-04-18

### Added
- Add `get_aipw_scores` for calculating AIPW scores given user-supplied estimates, and add some minor polish to `get_ipw_scores`. [#72](https://github.com/grf-labs/maq/pull/72), [#71](https://github.com/grf-labs/maq/pull/71)

### Fixed
- Draw a horizontal line when constructing a new plot depending on user-supplied xlim. [#73](https://github.com/grf-labs/maq/pull/73)
- Fix a minor discrepancy in `integrated_difference`. [#48](https://github.com/grf-labs/maq/pull/48)

## [0.3.1] - 2023-10-08

### Fixed
- Fix a bug in `integrated_difference` where the AUC measure is wrong if \bar B exceeds the point at which the curve plateaus. [#44](https://github.com/grf-labs/maq/pull/44)

## [0.3.0] - 2023-09-10

### Changed (breaking)
- Change the `maq` function signature to make `budget` an optional argument. The default behavior (`budget = NULL`) is to fit the Qini curve up to a maximum spend/unit where each unit that is expected to benefit, is treated. [#41](https://github.com/grf-labs/maq/pull/41)

### Added
- Add `integrated_difference(object.lhs, object.rhs, spend)` for estimating the area between two Qini curves up to some maximum budget `spend`. [#42](https://github.com/grf-labs/maq/pull/42)

## [0.2.0] - 2023-09-03

### Added
- Add a `type = c("matrix", "vector")` option to `predict.maq`, optionally returning predictions in the set {0, 1, ..., K} if `type = "vector"`. [#29](https://github.com/grf-labs/maq/pull/29)
- Add a convenience function `get_ipw_scores` to construct evaluation scores via IPW. [#28](https://github.com/grf-labs/maq/pull/28)
- Add more documentation on statistical details from the paper. [#27](https://github.com/grf-labs/maq/pull/27)

### Changed
- Have `predict.maq` return a standard dense matrix, and remove dependence on the sparse Matrix package. [#30](https://github.com/grf-labs/maq/pull/30)

### Fixed
- Fix `horizontal.line` in plot, the option for whether a curve added to the plot should extend all the way to the right of a main plot, if the added Qini curve's spend path stops before the spend path of the main plot.

## [0.1.0] - 2023-06-27
First CRAN beta release (this changelog tracks the R package). The R package currently supports

- Fitting Qini curves for an arbitrary number of arms using `maq(cate.hat, cost.hat, max.budget, Y.eval, ...)`, where `cate.hat` and `cost.hat` are CATE and cost estimates obtained via some function learned on a training set, `Y.eval` are evaluation scores on a test set (for example inverse-propensity weighted outcomes) and `max.budget` is the maximum mean budget/unit to fit the curve on. Setting the option `target.with.covariates` to `FALSE` yields a baseline Qini curve that can be used to assess the value of treatment targeting based on covariates.
- Computing point estimates and standard errors for any point on the curve with `average_gain()` (supporting clustered standard errors if fit with clusters).
- Comparing arbitrary points on different curves using `difference_gain()`, yielding standard errors that account for the correlation arising from curves fit on the same evaluation data.
- Retrieving the underlying "induced" policy at arbitrary spend levels with `predict()`.
- Basic plotting functionality with `plot()`.
- Retrieving the full gain/spend/allocation path with `summary(maq.object)`.

The Python bindings (source install only) currently supports

- Fitting the Qini curve, and retrieving standard errors for any point on the curve, as well as the implied policy \pi (remaining functionality can be added on the wrapper side - the heavy lifting is done in core/C++, and the input/output is language agnostic).
