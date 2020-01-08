

#' Early Treatment Diabetic Retinopathy Study
#' 
#' Enrolled 3711 patients with mild-to-severe nonproliferative or early proliferative 
#' diabetic retinopathy in both eyes from April 1980 to July 1985 - One eye per patient 
#' was randomized to early photocoagulation - The other eye deferred photocoagulation 
#' until detection of high-risk proliferative retinopathy. 
#' Survival endpoint: Time to severe visual loss (defined) as visual acuity < 5/00 
#' at two consecutive visits) subject to censoring eventually detected a beneﬁt with early photocoagulation 
#' after 9 years of follow-up, only 5.5% of patients had experienced the event (lots of censoring)
#' 
#' @format A \code{data.frame} with 3711 observations (rows) and 6 variables (columns):
#' \describe{
#'   \item{trt1}{numeric, indicator for treatment group 1}
#'   \item{x1}{numeric, time to event for treatment group 1}
#'   \item{delta1}{numeric, 1 if subject died, 0 if censored}
#'   \item{trt2}{numeric, indicator for treatment group 2}
#'   \item{x2}{numeric, time to event for treatment group 2}
#'   \item{delta2}{numeric, 1 if subject died, 0 if censored}
#' }
"pairdata"


#' Repeated time-to-event data with two treatment groups
#' 
#' Example data set for the paper: Tayob, N. and Murray, S., 2014. Nonparametric tests 
#' of treatment effect based on combined endpoints for mortality and recurrent events. 
#' Biostatistics, 16(1), pp.73-83.
#' 
#' The data set contains repeated measures time-to-event data with two treatment 
#' groups and a single continuous covariate. Summary of the data set:
#' 
#' @format A \code{data.frame} with 200 observations (rows) and 10 variables (columns):
#' \describe{
#'   \item{Subject_ID}{numeric, unique subject identifier}
#'   \item{X}{numeric, continuous covariate}
#'   \item{delta}{numeric, 1 if subject died, 0 if censored}
#'   \item{Group}{numeric, treatment group}
#'   \item{Z1}{numeric, time to first follow-up}
#'   \item{Z2}{numeric, time to second follow-up}
#'   \item{Z3}{numeric, time to third follow-up}
#'   \item{Z4}{numeric, time to fourth follow-up}
#'   \item{Z5}{numeric, time to fifth follow-up}
#' }
"TMdata"


#' Repeated time-to-event data with censoring
#' 
#' Example data set for the paper: Tayob, N. and Murray, S., 2016. Nonparametric 
#' restricted mean analysis across multiple follow-up intervals. Statistics & probability 
#' letters, 109, pp.152-158.
#' 
#' @format A \code{data.frame} with 100 observations (rows) and 2 variables (columns):
#' \describe{
#'   \item{X}{numeric, time to follow-up or event}
#'   \item{delta}{numeric, 1 if subject died, 0 if censored}
#' }
"TM2data"
