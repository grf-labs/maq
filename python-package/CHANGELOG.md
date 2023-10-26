# Changelog
All notable changes to `maq/python-package` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2023-10-25

### Added
- Added functionality to sync up features with R package.

## [0.1.0] - 2023-06-27
The basic Python bindings (source install only) currently supports

- Fitting the Qini curve, and retrieving standard errors for any point on the curve, as well as the implied policy \pi (remaining functionality can be added on the wrapper side - the heavy lifting is done in core/C++, and the input/output is language agnostic).
