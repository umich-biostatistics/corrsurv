% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TM.R
\name{TM}
\alias{TM}
\title{Tayob and Murray two-sample recurrent events test}
\usage{
TM(X, delta, Z, Group, Tau, t, method = "average")
}
\arguments{
\item{X}{vector containing observed time to terminating event}

\item{delta}{vector that is 1 if terminating event is observed and 0 if patient is censored}

\item{Z}{array of times to recurrent events. Number of columns corresponds to
maximum number of recurrent events observed for a patient}

\item{Group}{vector indicating which group each patient belongs to}

\item{Tau}{upper limit of integration (cannot be greater than largest follow-up
time, cannot be negative)}

\item{t}{vector containing start times of follow-up windows chosen}

\item{method}{Choose "average" for mean difference test or "area" for area between
the RMRL curves}
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
Perform the Tayob and Murray two-sample recurrent events test described in
Tayob, N. and Murray, S., 2014. Nonparametric tests of treatment
effect based on combined endpoints for mortality and recurrent events. Biostatistics,
16(1), pp.73-83.
}
\examples{
data("TMdata")
head(TMdata)

N=nrow(TMdata) #Number of subjects
X=TMdata$X
delta=TMdata$delta
table(delta)
Z=array(cbind(TMdata$Z1,TMdata$Z2,TMdata$Z3,TMdata$Z4,TMdata$Z5,TMdata$Z6),c(N,6))
A=max(X) #length of study
Treatment=as.numeric(TMdata$Group==1) #1 if Case and 0 if Control
table(Treatment,delta)

srec.average = TM(X = X, delta = delta, Z = Z, Group = Treatment, Tau = A/2, 
                  t = seq(from = 0, to = A-A/2, by = A/4))
summary(srec.average)

srec.area = TM(X = X, delta = delta, Z = Z, Group = Treatment, Tau = A/2, 
               t = seq(from = 0, to = A-A/2, by = A/4), method = "area")
summary(srec.area)
plot(srec.area)

}
\references{
Tayob, N. and Murray, S., 2014. Nonparametric tests of treatment
effect based on combined endpoints for mortality and recurrent events. Biostatistics,
16(1), pp.73-83.
}
\author{
Nabihah Tayob
}
