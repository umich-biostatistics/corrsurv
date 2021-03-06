% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TM.R
\name{format_data_ourmethod}
\alias{format_data_ourmethod}
\title{Format observed data for TM (Tayob and Murray method)}
\usage{
format_data_ourmethod(X, delta, Z, t)
}
\arguments{
\item{X}{vector containing observed time to terminating event}

\item{delta}{vector that is 1 if terminating event is observed and 0 if patient
is censored}

\item{Z}{array of times to recurrent events. Number of columns corresponds to maximum
number of recurrent events observed for a patient}

\item{t}{vector containing start times of follow-up windows chosen}
}
\value{
a \code{list} containing
\itemize{
\item{ID - vector containing numeric identifying numbers for each patient corresponding to row numbers of input data}
\item{X_tk - vector containing observed time to combined end-point for each follow-up window}
\item{delta_tk - vector containing combined end-point event indicator. 1 if combined-endpoint is observed and 0 if censored}
\item{t_k - vector containing start times of corresponding follow-up window}
\item{n - number of patients in dataset}
}
}
\description{
Formats the data for the method described in Tayob, Murray 2015. This function
is for internal use by TM().
}
\references{
Tayob, N. and Murray, S., 2014. Nonparametric tests of treatment
effect based on combined endpoints for mortality and recurrent events. Biostatistics,
16(1), pp.73-83.
}
\author{
Nabihah Tayob
}
