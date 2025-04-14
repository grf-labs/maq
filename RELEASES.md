# maq Release Process

Open a “Prepare the v.x. release PR”

**1) Prepare release notes**

Go over PRs since last PR named “Prepare the x.x.x release.”, add all changes to CHANGELOG.
Make different notes for R and Python.

**Bump version in DESCRIPTION**

Bump patch release version (1.0.0 → 1.0.1), or minor version if it’s a bigger release (1.0.2 → 1.1.0).

**Python**

Bump the version in `setup.py` and in Python's `README.md` install link. Update the Python CHANGELOG.
If breaking Python changes, make sure to reflect these in the causalML notebook too. Reflect Python versions with a pyx.x.x tag.

**2) Make sure there are no performance regressions**

See below for example scripts to run.

**3) Create release tarball**

Add
```
^tests/testthat/test_((?!cran).).*
```
to `.Rbuildignore` (do not commit the change)

Run

`R CMD build .`

**4) Run CRAN checks**

Run `R CMD check --as-cran --run-donttest <release tarball>`.

**Check package using win-builder**

Upload the release tarball to https://win-builder.r-project.org/upload.aspx to check it against the latest R-devel and also R-release on Windows. To test a larger release extra carefully on even more systems you can check it here as well: https://builder.r-hub.io/.  First upload a tarball that includes all the R tests not intended to be run on CRAN.

**5) Submit to CRAN**

Upload package here: https://cran.r-project.org/submit.html

Add the tarball to local release folder (not committed)

Tag the release

```
git tag -a v1.1.0 5342d1de47f7604d86c4e59825bb7e4928a97c16
git push --tags
```

Where the hash is for the commit “Prepare the x.x.x release”.

## Make sure there are no performance regressions

**1) Check timing**

File: https://github.com/grf-labs/maq/blob/master/r-package/maq/tests/benchmarks/benchmark.R


**2) Check memory usage**

Run the massif tool from valgrind on a simple script.

```
> R -d "valgrind --tool=massif" -f valgrind.R
> ms_print massif.out.34681 &> ms_print.out
```

The script `valgrind.R` is defined as follows:

```R
library(maq)
budget <- 1000
n <- 15000
K <- 5
reward <- matrix(0.1 + rnorm(n * K), n, K)
cost <- 0.05 + matrix(runif(n * K), n, K)
mq <- maq(reward, cost, budget, reward, R = 150)
```

For a big release, can also run `R CMD check --as-cran --run-donttest --use-valgrind <development tarball>` on a Linux machine with Valgrind for a thorough stress test (will take very long).

## Previous performance test results

**0.6.0**

Only a minor release with no perf/C++ touches.

**0.5.0**

Only a minor release with no perf/C++ touches.

**0.4.0**

Only a minor release with no perf/C++ touches.

**0.3.1**

Only a patch release touching integrated_difference

**0.3.0**

Only a minor release with no perf/C++ touches.

**0.2.0**

Only a minor release with no perf/C++ touches.

**0.1.0**

(Machine: 2015 Macbook Pro 8 cores, R version 4.2.2)

perf

```
# > print(list(b1,b2,b3,b41,b42,b5), digits = 3)
# [[1]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, reward, R = 200, paired.inference = FALSE) 1.93 1.97    2   1.99 2.02 2.12    10
#
# [[2]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, reward, R = 200, paired.inference = FALSE) 3.84 3.85 3.87   3.87 3.89 3.92    10
#
# [[3]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, reward, R = 200, paired.inference = FALSE) 15.3 15.3 15.4   15.3 15.3 15.8     5
#
# [[4]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, reward, R = 0, paired.inference = FALSE) 1.91 1.92 1.94   1.93 1.95 2.01    10
#
# [[5]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, reward, R = 200, paired.inference = FALSE) 31.7 31.8 32.9   32.1 32.4 40.7    10
#
# [[6]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, reward, R = 0, paired.inference = FALSE) 10.2 10.2 10.2   10.2 10.2 10.2     5
```

memory usage

```
    MB
201.6^                                                                       #
     |                                                                       #
     |                                                                      :#
     |                                                      @:@:@:@@:@::::@::#
     |                                                :@   :@:@:@:@@:@::::@::#
     |                                               @@@::::@:@:@:@@:@::::@::#
     |                                             @@@@@::::@:@:@:@@:@::::@::#
     |                                            @@@@@@::::@:@:@:@@:@::::@::#
     |                                           @@@@@@@::::@:@:@:@@:@::::@::#
     |                          :    :         :@@@@@@@@::::@:@:@:@@:@::::@::#
     |                         :: ::::::::   :::@@@@@@@@::::@:@:@:@@:@::::@::#
     |                       @::::::::: ::::::::@@@@@@@@::::@:@:@:@@:@::::@::#
     |                      :@::::::::: ::: ::::@@@@@@@@::::@:@:@:@@:@::::@::#
     |                     @:@::::::::: ::: ::::@@@@@@@@::::@:@:@:@@:@::::@::#
     |            @:   ::@:@:@::::::::: ::: ::::@@@@@@@@::::@:@:@:@@:@::::@::#
     |         :@:@:::@::@:@:@::::::::: ::: ::::@@@@@@@@::::@:@:@:@@:@::::@::#
     |        ::@:@:::@::@:@:@::::::::: ::: ::::@@@@@@@@::::@:@:@:@@:@::::@::#
     |    ::::::@:@:::@::@:@:@::::::::: ::: ::::@@@@@@@@::::@:@:@:@@:@::::@::#
     |   ::: :::@:@:::@::@:@:@::::::::: ::: ::::@@@@@@@@::::@:@:@:@@:@::::@::#
     | @@::: :::@:@:::@::@:@:@::::::::: ::: ::::@@@@@@@@::::@:@:@:@@:@::::@::#
   0 +----------------------------------------------------------------------->Gi
     0                                                                   5.330
```
