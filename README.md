Get started with pairsurv
================

## About

This package contains tools for analyzing paired survival data. Three
methods are contained for the assessment of treatment effects in paired
survival data are approached using the methods of Murray (2000), Tayob,
Murray (2015), and Tayob, Murray (2016). The details of these methods
are introduced below. See the package documentation for additional
details.

## Install

The package is currently hosted on Github. If you have not already, you
will need to install the devtools package before using the package:

``` r
install.packages("devtools")
```

Once installed, simply run the following in R:

``` r
devtools::install_github("umich-biostatistics/pairsurv")
```

## Usage

#### Method of Murray, 2000 (pairtest)

For the method of Murray (2000), use the `pairtest()` function.

Here is an example R script using the `pairdata.csv` data set.

``` r
library(pairsurv)

data(pairdata)

eyeresults=matched.survtest.9.10.99(pairdata$x1, pairdata$delta1, pairdata$x2, pairdata$delta2, 3711)

(integration_upper_limit= eyeresults$upperlim)

(logrank_statistic=eyeresults$lgrk.stat.p)
(logrank_statistic_p=2*(pnorm(logrank_statistic, lower.tail=T)) )#since negative statistic

(gehan_statistic=eyeresults$gehan.stat.p)
(gehan_statistic_p=2*(pnorm(gehan_statistic, lower.tail=T)) )#since negative statistic

(yls_statistic=eyeresults$yls.stat.p)
(yls_statistic_p=2*(pnorm(yls_statistic, lower.tail=F))) #since positive statistic

(pepe_flem_statistic=eyeresults$pf.stat.p)
(pepe_flem_statistic_p=2*(pnorm(pepe_flem_statistic, lower.tail=F))) #since positive statistic


(logrank_assuming_indep=eyeresults$ lgrk.nopair.stat.p)
(logrank_assuming_indep_p=2*(pnorm(logrank_assuming_indep, lower.tail=T)) )#since negative statistic

(gehan_assuming_indep= eyeresults$gehan.nopair.stat.p)
(gehan_assuming_indep_p=2*(pnorm(gehan_assuming_indep, lower.tail=T))) #since negative statistic

(yls_assuming_indep= eyeresults$ yls.nopair.stat.p)
(yls_assuming_indep_p=2*(pnorm(yls_assuming_indep, lower.tail=F))) #since positive statistic

(pf_assuming_indep= eyeresults$ pf.nopair.stat.p)
(pf_assuming_indep_p=2*(pnorm(pf_assuming_indep, lower.tail=F))) #since positive statistic


nstuff=((3711*3711)/(3711+3711))^(.5)


#estimate area between survival curves and 95% conf. int.
(yls.diff=eyeresults$yls.num/nstuff)
(yls.upper=yls.diff+1.96*sqrt(eyeresults$yls.type.var)/nstuff)
(yls.lower=yls.diff-1.96*sqrt(eyeresults$yls.type.var)/nstuff)


#95% c.i. if assuming independence
(yls.upper2=yls.diff+1.96*sqrt(eyeresults$yls.nopair.type.var)/nstuff)
(yls.lower2=yls.diff-1.96*sqrt(eyeresults$yls.nopair.type.var)/nstuff)
```

#### Method of Tayob and Murray, 2015 (TM)

#### Method of Tayob and Murray, 2016 (TM2)
