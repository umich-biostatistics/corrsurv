Get started with *corrsurv*
================

## About

Collection of two-sample tests for treatment effects with paired
censored survival data and recurrent events survival data. The methods
are implementations of three papers by Susan Murray and Nabihah Tayob.

## Install

First, install and load the `corrsurv` package from Github using the
following code:

If you have not installed devtools, run the following:

``` r
install.packages("devtools")
```

If devtools is sucessfully installed, run:

``` r
if(!require(corrsurv)) {
  library(devtools)
  install_github('umich-biostatistics/corrsurv') 
}
```

Load the package:

``` r
library(corrsurv)
```

## Summary of methods

This package implements methods from the following three papers by Tayob
and Murray:

### Method 1: Murray, 2000

##### Description

Perform two-sample tests for treatment effects with paired censored
survival data.

##### Reference

Murray, Susan. Nonparametric Rank-Based Methods for Group Sequential
Monitoring of Paired Censored Survival Data. 2000. Biometrics, 56,
pp. 984-990.

Jump to the `pairtest()` section for the methods from this paper.

### Method 2: Tayob and Murray, 2014

##### Description

Perform the Tayob and Murray two-sample recurrent events test.

##### Reference

Tayob, N. and Murray, S., 2014. Nonparametric tests of treatment effect
based on combined endpoints for mortality and recurrent events.
Biostatistics, 16(1), pp.73-83.

Jump to the `TM()` section for the methods from this paper.

### Method 3: Tayob and Murray, 2016

##### Description

Estimate the tau-restricted mean survival across multiple follow-up
intervals.

##### Reference

Tayob, N. and Murray, S., 2016. Nonparametric restricted mean analysis
across multiple follow-up intervals. Statistics & probability letters,
109, pp.152-158.

Jump to the `TM2()` section for the methods from this paper.

## Results

``` r
vignette("learn-corrsurv")
```
