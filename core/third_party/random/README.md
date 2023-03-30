# third_party/random

These are the llvm `random` and `algorithm` headers bundled to avoid compiler dependency on random number generation. This means that with the same random seed, MAQ will deliver exactly the same bootstrapped standard errors across platforms, which is convenient for reproducibility. See https://github.com/grf-labs/grf/tree/master/core/third_party/random for details.
