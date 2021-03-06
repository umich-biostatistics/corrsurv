% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TM2.R
\name{TM2}
\alias{TM2}
\title{Tau-Restricted Mean Survival}
\usage{
TM2(X, delta, Tau, t, var_output = "proposed")
}
\arguments{
\item{X}{follow-up time of right censored data (values must be geq 0)}

\item{delta}{status indicator, 0=alive, 1=dead. (values must be 0,1)}

\item{Tau}{upper limit of integration (cannot be greater than largest follow-up
time, cannot be negative)}

\item{t}{start times of follow-up windows (default=seq(from=0, to=A-Tau,by=(A-Tau)/(b-1)),
must be of length b if both specified, largest value cannot be greater than A-Tau,
no repeats)}

\item{var_output}{Type of variance estimator. Options are c("proposed","independence",
"sandwich","all")}
}
\value{
A \code{list} object which contains
\itemize{
\item{Mean - vector containing sample estimates of overall tau-restricted mean survival in each group}
\item{Var - vector containing empirical variance of estimates of overall tau-restricted mean survival in each group}
\item{test_stat - test statistic of two-sample test}
\item{test_stat_p - p-value of two-sample test}
}
}
\description{
Estimate the tau-restricted mean survival across multiple follow-up intervals.
}
\details{
Estimate the tau-restricted mean survival as described in Tayob, N. and Murray, S., 2016.
Nonparametric restricted mean analysis across multiple follow-up intervals. Statistics
& probability letters, 109, pp.152-158.
}
\examples{
data("TM2data")

output = TM2(X = TM2data$X, delta = TM2data$delta, Tau = 12, 
             t = seq(from = 0, to = 24, by = 6), var_output = "all")
summary(output)
plot(output)

output = TM2(X = TM2data$X, delta = TM2data$delta, Tau = 12, 
             t = seq(from = 0, to = 24, by = 6))
summary(output)
plot(output)

output = TM2(X = TM2data$X, delta = TM2data$delta, Tau = 12, 
             t = seq(from = 0, to = 24, by = 6), var_output = "independence")
summary(output)
plot(output)

output = TM2(X = TM2data$X, delta = TM2data$delta, Tau = 12, 
             t = seq(from = 0, to = 24, by = 6), var_output = "sandwich")
summary(output)
plot(output)

}
\references{
Tayob, N. and Murray, S., 2016. Nonparametric restricted mean analysis
across multiple follow-up intervals. Statistics & probability letters, 109, pp.152-158.
}
\author{
Nabihah Tayob
}
