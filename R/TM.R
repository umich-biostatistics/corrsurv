

#' Helper function for column sums of upper trianglular matrix
#' 
#' @param X matrix
#' @author Nabihah Tayob
#' 
#' @references Tayob, N. and Murray, S., 2014. Nonparametric tests of treatment 
#' effect based on combined endpoints for mortality and recurrent events. Biostatistics, 
#' 16(1), pp.73-83.

sum_function = function(X) {
  temp = array(X, c(length(X), length(X)))
  gdata::lowerTriangle(temp, diag = F) = 0
  apply(temp, 2, sum)
}


#' Format observed data for TM (Tayob and Murray method)
#' 
#' Formats the data for the method described in Tayob, Murray 2015. This function
#' is for internal use by TM().
#' 
#' @param X vector containing observed time to terminating event
#' @param delta vector that is 1 if terminating event is observed and 0 if patient 
#' is censored
#' @param Z array of times to recurrent events. Number of columns corresponds to maximum 
#' number of recurrent events observed for a patient
#' @param t vector containing start times of follow-up windows chosen
#'
#' @return a \code{list} containing
#' \itemize{
#'   \item{ID - vector containing numeric identifying numbers for each patient corresponding to row numbers of input data}
#'   \item{X_tk - vector containing observed time to combined end-point for each follow-up window}
#'   \item{delta_tk - vector containing combined end-point event indicator. 1 if combined-endpoint is observed and 0 if censored}
#'   \item{t_k - vector containing start times of corresponding follow-up window}
#'   \item{n - number of patients in dataset}
#' }
#' 
#' @author Nabihah Tayob
#' 
#' @references Tayob, N. and Murray, S., 2014. Nonparametric tests of treatment 
#' effect based on combined endpoints for mortality and recurrent events. Biostatistics, 
#' 16(1), pp.73-83.

format_data_ourmethod = function(X, delta, Z, t) {
  
  Z_star = Z #rename data
  n = length(X)
  b = length(t)
  
  X_tk = array(NA, c(n, b))
  delta_tk = array(NA, c(n, b))
  
  for(k in 1:b) {
    Z_star_k = Z_star - t[k]
    Z_star_k[Z_star_k < 0] = NA #recurrent events observed before time t[k]
    
    X_k = X - t[k]
    delta_k = delta
    delta_k[X_k < 0] = 0 #terminating events observed before time t[k]->patient censored at 0
    X_k[X_k < 0] = 0
    
    X_Z_star_k = apply(cbind(Z_star_k, X_k), 1, min, na.rm = T) #time to next recurrent event/terminating event
    X_tk[,k] = X_Z_star_k
    
    delta_k2 = as.numeric(X_Z_star_k < X_k) #=1 if event observed was a recurrent event
    delta_k2[X_Z_star_k == X_k] = delta_k[X_Z_star_k == X_k] #=delta for terminating event
    delta_tk[,k] = delta_k2
  }
  
  t_k = rep(t, times = rep(n, length(t)))
  X_tk = array(X_tk, c(n*b, 1))
  delta_tk = array(delta_tk, c(n*b, 1))
  ID = rep(seq(1:n), b)
  list(ID = ID, X_tk = X_tk, delta_tk = delta_tk, t_k = t_k, n = n)
  
}


#' Estimate and variance for Tayob and Murray
#' 
#' Get estimate and variance of estimate for method described in Tayob, Murray 2015. 
#' This function is used internally when TM() is called.
#'
#' @param X_km vector containing observed time to combined end-point for each follow-up window
#' @param delta_km vector containing combined end-point event indicator. 1 if 
#' combined-endpoint is observed and 0 if censored
#' @param Tau length of follow-up intervals of interest
#' @param t vector containing start times of follow-up windows chosen
#'
#' @return a \code{list} containing
#' \itemize{
#'   \item{mean - estimate of overall tau restricted mean survival}
#'   \item{var - empirical variance estimate of mean}
#' }
#' 
#' @author Nabihah Tayob
#' 
#' @references Tayob, N. and Murray, S., 2014. Nonparametric tests of treatment 
#' effect based on combined endpoints for mortality and recurrent events. Biostatistics, 
#' 16(1), pp.73-83.

get_mu_hat_star_tau = function(X_km, delta_km, Tau, t) {
  
  b = length(t)
  n = length(X_km) / b
  X = array(X_km, c(n, b))
  delta = array(delta_km, c(n, b))
  
  observed_events = sort(X*delta)
  T = unique(c(observed_events[observed_events <= Tau], Tau))
  M = length(T) - 1
  T_array = round(t(array(T[1:M], c(M, b))), 6)
  
  dN_i = array(NA, c(n, M))
  Y_i = array(NA, c(n, M))
  
  for(j in 1:n) {
    temp1 = array(round(X[j,], 6), c(b, M))
    temp2 = array(delta[j,], c(b, M))
    dN_i[j,] = apply((temp1 == T_array & temp2 == 1), 2, sum) 
    Y_i[j,] = apply(temp1 >= T_array, 2, sum)
  }
  
  dN = apply(dN_i, 2, sum)
  Y = apply(Y_i, 2, sum)
  
  time_int = T[2:(M+1)] - T[1:M]
  temp = array(dN/Y, c(M, M))
  gdata::lowerTriangle(temp, diag = F) = 0
  CH = apply(temp, 2, sum)
  S_hat = exp(-CH)
  mean = sum(time_int*S_hat)
  
  q = dN / Y
  z_i_q = t(t(dN_i - t(q*t(Y_i)))*(n/Y))
  temp3 = t(apply(z_i_q, 1, sum_function))
  time_int_array = t(array(time_int, c(M, n)))
  S_hat_array = t(array(S_hat, c(M, n)))
  z_i_ATau = apply(time_int_array*S_hat_array*temp3, 1, sum)
  williams_var = var(z_i_ATau) / n
  
  list(mean = mean, var = williams_var)
  
}


#' Calculate RMRL
#' 
#' Area under RMRL and variance of area under RMRL. 
#' RMRL is the restricted mean residual life function evaluated at times defined by t.
#' This function is used internally when TM() is called.
#'
#' @param data_format output of function "format_data_ourmethod"
#' @param Tau length of follow-up intervals of interest
#' @param t vector containing start times of follow-up windows chosen
#' 
#' @return a \code{list} containing
#' \itemize{
#'   \item{RMRL - restricted mean residual life function evaluated at times defined by t}
#'   \item{area_under_RMRL - area under RMRL function evaluated at times t}
#'   \item{var_area_under_RMRL - empirical variance of estimate area_under_RMRL}
#' }
#' 
#' @author Nabihah Tayob
#'
#' @references Tayob, N. and Murray, S., 2014. Nonparametric tests of treatment 
#' effect based on combined endpoints for mortality and recurrent events. Biostatistics, 
#' 16(1), pp.73-83.

RMRL_function = function(data_format, Tau, t) {
  
  b = length(t)
  data_tk = data.frame(data_format)
  t_k = data_format$t_k
  time = unique(sort(data_tk$X_tk*as.numeric(data_tk$X_tk <= Tau)))
  M_time = length(time)
  time_int_star = time[2:M_time] - time[1:(M_time-1)]
  RMRL = rep(NA, b)
  Z_i_tk = NULL
  
  for(i in 1:b) {
    data_tk_subset = subset(data_tk, t_k == t[i])
    n = nrow(data_tk_subset)
    X = data_tk_subset$X_tk
    delta = data_tk_subset$delta_tk
    
    #calculate components of variance of area under RMRL
    time_array = t(array(time, c(M_time, n)))
    dN_i_tk = array(as.numeric(X == time_array & delta == 1), c(n, M_time))
    Y_i_tk = array(as.numeric(X >= time_array), c(n, M_time))
    
    dN_tk = apply(dN_i_tk, 2, sum)
    Y_tk = apply(Y_i_tk, 2, sum)
    lambda_tk = dN_tk / Y_tk
    temp = array(lambda_tk, c(M_time, M_time))
    gdata::lowerTriangle(temp, diag = F) = 0
    S_tk_m = exp(-apply(temp, 2, sum))
    mean = sum(time_int_star*S_tk_m[1:(M_time-1)])
    RMRL[i] = mean
    
    temp1 = t(array(lambda_tk, c(M_time, n)))
    temp2 = t(array(Y_tk/n, c(M_time, n)))
    temp3 = (dN_i_tk - temp1*Y_i_tk)/temp2
    Z_i_tk_m = t(apply(temp3, 1, sum_function))
    temp4 = t(array(S_tk_m, c(M_time, n)))
    temp5 = (temp4*Z_i_tk_m)[,1:(M_time-1)]
    temp6 = t(array(time_int_star, c(M_time-1, n)))*temp5
    Z_i_tk = cbind(Z_i_tk, apply(temp6, 1, sum))
  }
  
  area_under_RMRL = sum((t[2:b] - t[1:(b-1)])*(RMRL[2:b] + RMRL[1:(b-1)])/2)
  
  #variance of area under RMRL
  Z_i = rep(0, nrow(Z_i_tk))
  
  for(i in 2:b) {
    Z_i = Z_i + (t[i] - t[i-1])*(Z_i_tk[,i] + Z_i_tk[,i-1])/2
  }
  var_area_under_RMRL = var(Z_i)/n
  
  list(RMRL = RMRL, area_under_RMRL = area_under_RMRL, 
       var_area_under_RMRL = var_area_under_RMRL)
  
}  


#' Tayob and Murray two-sample recurrent events test
#' 
#' Perform the Tayob and Murray two-sample recurrent events test described in 
#' Tayob, N. and Murray, S., 2014. Nonparametric tests of treatment 
#' effect based on combined endpoints for mortality and recurrent events. Biostatistics, 
#' 16(1), pp.73-83.
#' 
#' @param X vector containing observed time to terminating event
#' @param delta vector that is 1 if terminating event is observed and 0 if patient is censored
#' @param Z array of times to recurrent events. Number of columns corresponds to 
#' maximum number of recurrent events observed for a patient
#' @param Group vector indicating which group each patient belongs to
#' @param Tau upper limit of integration (cannot be greater than largest follow-up 
#' time, cannot be negative)
#' @param t vector containing start times of follow-up windows chosen
#' @param method Choose "average" for mean difference test or "area" for area between
#' the RMRL curves
#' 
#' @return A \code{list} object which contains
#' \itemize{
#'   \item{Mean - vector containing sample estimates of overall tau-restricted mean survival in each group}
#'   \item{Var - vector containing empirical variance of estimates of overall tau-restricted mean survival in each group}
#'   \item{test_stat - test statistic of two-sample test}
#'   \item{test_stat_p - p-value of two-sample test}
#' }
#' 
#' @examples
#' data("TMdata")
#' head(TMdata)
#' 
#' N=nrow(TMdata) #Number of subjects
#' X=TMdata$X
#' delta=TMdata$delta
#' table(delta)
#' Z=array(cbind(TMdata$Z1,TMdata$Z2,TMdata$Z3,TMdata$Z4,TMdata$Z5,TMdata$Z6),c(N,6))
#' A=max(X) #length of study
#' Treatment=as.numeric(TMdata$Group==1) #1 if Case and 0 if Control
#' table(Treatment,delta)
#' 
#' srec.average = TM(X = X, delta = delta, Z = Z, Group = Treatment, Tau = A/2, 
#'                   t = seq(from = 0, to = A-A/2, by = A/4))
#' summary(srec.average)
#' 
#' srec.area = TM(X = X, delta = delta, Z = Z, Group = Treatment, Tau = A/2, 
#'                t = seq(from = 0, to = A-A/2, by = A/4), method = "area")
#' summary(srec.area)
#' plot(srec.area)
#' 
#' @author Nabihah Tayob
#' 
#' @references Tayob, N. and Murray, S., 2014. Nonparametric tests of treatment 
#' effect based on combined endpoints for mortality and recurrent events. Biostatistics, 
#' 16(1), pp.73-83.

TM = function(X, delta, Z, Group, Tau, t, method = "average") {
  
  args = match.call()
  args$method = method
  
  Group_id = unique(Group)
  formatted_group1_data = format_data_ourmethod(X = X[Group == Group_id[1]],
                                                delta = delta[Group == Group_id[1]],
                                                Z = Z[Group == Group_id[1],], t)
  formatted_group2_data = format_data_ourmethod(X = X[Group == Group_id[2]],
                                                delta = delta[Group == Group_id[2]],
                                                Z = Z[Group == Group_id[2],], t)  
  
  if(method == "average") {
    results1 = get_mu_hat_star_tau(formatted_group1_data$X_tk, 
                                   formatted_group1_data$delta_tk, Tau = Tau, t = t)
    results2 = get_mu_hat_star_tau(formatted_group2_data$X_tk, 
                                   formatted_group2_data$delta_tk, Tau = Tau, t = t)
    difference = results1$mean - results2$mean
    sd_diff = sqrt(results1$var + results2$var)
    test_stat = abs(difference / sd_diff)
    test_stat_p = 2*(1 - pnorm(test_stat))
    mean12 = c(results1$mean, results2$mean)
    var12 = c(results1$var, results2$var)
  }
  
  if(method == "area") {
    results1 = RMRL_function(data_format = formatted_group1_data, Tau = Tau, t = t)
    results2 = RMRL_function(data_format = formatted_group2_data, Tau = Tau, t = t)
    difference = results1$area_under_RMRL - results2$area_under_RMRL
    sd_diff = sqrt(results1$var_area_under_RMRL + results2$var_area_under_RMRL)
    test_stat = abs(difference / sd_diff)
    test_stat_p = 2*(1 - pnorm(test_stat))
    mean12 = c(results1$area_under_RMRL, results2$area_under_RMRL)
    var12 = c(results1$var_area_under_RMRL, results2$var_area_under_RMRL)   
  }
  
  result = 
    list(
      Mean = mean12, 
      Var = var12, 
      test_stat = test_stat, 
      test_stat_p = test_stat_p,
      args = args,
      Group_id = Group_id,
      formatted_group1_data = formatted_group1_data,
      formatted_group2_data = formatted_group2_data,
      results1 = results1,
      results2 = results2
    )
  
  attr(result, "class") = "TM"
  
  return(
    result
  )

}


#' Print the summary of a TM object
#' 
#' @param object an object of class 'TM'
#' @param digits number of digits to round to after decimal
#' @param ... additional options

summary.TM = function(object, digits = max(3, getOption("digits") - 3), ...) {
  
  test_stat = object$test_stat
  test_stat_p = object$test_stat_p
  args = object$args
  t = eval(args$t)
  Tau = eval(args$Tau)
  method = eval(args$method)
  Group_id = object$Group_id
  formatted_group1_data = object$formatted_group1_data
  formatted_group2_data = object$formatted_group2_data
  results1 = object$results1
  results2 = object$results2
  
  test_stat_p_print = paste("=", round2(test_stat_p, digits = digits))
  #if(round(test_stat_p, 5) == 0) { test_stat_p_print = "<0.00001" }
  follow_up_windows_print = NULL
  
  for(k in 1:length(t)) {
    follow_up_windows_print = paste(follow_up_windows_print, 
                                    paste("[", t[k], ",", t[k] + Tau, "]", sep = ""), "     ", sep = "")
  }
  
  cat("\n")
  cat(" ************ Two-Sample Test for combined end-point across multiple follow-up windows ************ ")
  cat("\n")
  cat("\n")
  cat(" Reference Paper: Nonparametric Tests of Treatment Effect for a Recurrent")
  cat("\n")
  cat(" Event process that Terminates - N Tayob and S Murray")
  cat("\n")
  cat("\n")
  cat(" Group definitions:")
  cat("\n")
  cat(paste("  Group 1: Group = ", Group_id[1], sep = ""))
  cat("\n")
  cat(paste("  Group 2: Group = ", Group_id[2], sep = ""))
  cat("\n")
  cat("\n")  
  cat(" Sample Size: ")
  cat("\n")
  cat(paste("  Group 1: ", formatted_group1_data$n, sep = ""))
  cat("\n")
  cat(paste("  Group 2: ", formatted_group2_data$n, sep = ""))
  cat("\n")
  cat("\n")
  cat(" Follow-up windows:")
  cat("\n")
  cat(paste("  ", follow_up_windows_print, sep = ""))
  cat("\n")
  cat("\n")
  
  if(method == "average") {
    cat(" Sample Estimates: Average restricted mean survival across follow-up windows")
    cat("\n")
    cat(paste("  Group 1: ", round(results1$mean, digits = digits), sep = ""))
    cat("\n")
    cat(paste("  Group 2: ", round(results2$mean, digits = digits), sep = ""))
    cat("\n")
    cat("\n")
    cat("  Alternative hypothesis: True difference in the average restricted mean 
        survival across follow-up windows is not equal to 0")
    cat("\n")
  }
  
  if(method == "area") {
    cat(" Sample Estimates: Area under RMRL function evaluated for each follow-up window")
    cat("\n")
    cat(paste("  Group 1: ", round(results1$area_under_RMRL, digits = digits), sep = ""))
    cat("\n")
    cat(paste("  Group 2: ", round(results2$area_under_RMRL, digits = digits), sep = ""))
    cat("\n")
    cat("\n")
    cat(" Alternative hypothesis: True difference in the area under the RMRL functions 
        is not equal to 0")
    cat("\n")
  }
  cat(paste("   Test-statistic = ", round(test_stat, digits = digits), ", p-value ", test_stat_p_print, sep = ""))
  
}


#' Plot a TM object
#' 
#' @param x an object of class 'TM'
#' @param ... addtional options to be passed on

plot.TM = function(x, ...) {
  
  args = x$args
  method = eval(args$method)
  t = eval(args$t)
  Tau = eval(args$Tau)
  results1 = x$results1
  results2 = x$results2
  
  if(method != "area") { stop("Plotting is only implemented for the 'area' method") }
  
  RMRL1 = results1$RMRL
  RMRL2 = results2$RMRL
  xadj = Tau/50
  yadj1 = max(RMRL1, RMRL2)/20
  yadj2 = max(RMRL1, RMRL2)/50
  matplot(cbind(t, t), cbind(RMRL1, RMRL2), type = "l", xaxs = "i", xaxt = "n", 
          yaxs = "i", xlab = "Time", ylab= "RMRL", xlim = c(t[1] - xadj, t[length(t)] + xadj), 
          ylim = c(min(RMRL1, RMRL2) - yadj1, max(RMRL1, RMRL2) + yadj2))
  points(t, RMRL1, pch = 20, cex = 1, col = "black")
  points(t, RMRL2, pch = 20, cex = 1, col = "red")
  axis(1, at = t, labels = t)
  legend("bottomright", c("Group 1", "Group 2"), col = c("black", "red"), 
          lty = c(1, 2), bty = "n")
  
}


#' Format observed data for Andersen-Gill method
#'
#' @param X vector containing observed time to terminating event
#' @param delta vector that is 1 if terminating event is observed and 0 if patient is censored
#' @param Z_star FILL IN
#'
#' @return A \code{list} object which contains
#' \itemize{
#'   \item{ID - Patient ID}
#'   \item{T_start - start time}
#'   \item{T_stop - Stop time}
#'   \item{status - status indicator: 1 if event occured, 0 if censored}
#' }
#' 
#' @author Nabihah Tayob
#' 
#' @references Tayob, N. and Murray, S., 2014. Nonparametric tests of treatment 
#' effect based on combined endpoints for mortality and recurrent events. Biostatistics, 
#' 16(1), pp.73-83.

format_data_AG = function(X, delta, Z_star) {
  
  n = length(X)
  
  X_Z_star = cbind(Z_star, X)
  T_start = NULL
  T_stop = NULL
  status = NULL
  ID_AG = NULL
  for(i in 1:n) {
    T_stop_subject = as.numeric(na.omit(X_Z_star[i,]))
    n_reccurent_subject = length(T_stop_subject) - 1
    
    if(n_reccurent_subject > 0) {
      T_start_subject = c(0, T_stop_subject[1:n_reccurent_subject])
    }
    
    if(n_reccurent_subject == 0) {
      T_start_subject = 0
    }    
    status_subject = c(rep(1, n_reccurent_subject), delta[i])
    
    T_start = c(T_start, T_start_subject)
    T_stop = c(T_stop, T_stop_subject)
    status = c(status, status_subject)  
    ID_AG = c(ID_AG, rep(i, n_reccurent_subject + 1))
    
    if(length(T_start) != length(T_stop)) { print(i) }
  }
  
  return(
    list(
      ID = ID_AG,
      T_start = T_start, 
      T_stop = T_stop, 
      status = status
    )
  )
  
}


#' Andersen-Gill test
#' 
#' Andersen-Gill model based two-sample test.
#'
#' @param data1_format Data for sample 1 formatted using \code{format_data_AG}
#' @param data2_format Data for sample 2 formatted using \code{format_data_AG}
#' 
#' @return A \code{list} object which contains
#' \itemize{
#'   \item{test_stat_p - p-value for the test}
#' }
#' 
#' @author Nabihah Tayob
#'
#' @references Tayob, N. and Murray, S., 2014. Nonparametric tests of treatment 
#' effect based on combined endpoints for mortality and recurrent events. Biostatistics, 
#' 16(1), pp.73-83.

AG_test_stat = function(data1_format, data2_format) {
  
  n1 = length(data1_format$ID)
  treatment = c(rep(0, length(data1_format$ID)), rep(1, length(data2_format$ID)))
  T_start = c(data1_format$T_start, data2_format$T_start)
  T_stop = c(data1_format$T_stop, data2_format$T_stop)
  status = c(data1_format$status, data2_format$status)
  id = c(data1_format$ID, data2_format$ID + n1)
  ag_model = coxph(Surv(T_start, T_stop, status) ~ treatment + cluster(id))
  test_stat_p = 1 - pchisq(ag_model$wald.test, 1)
  
  return(
    list(
      test_stat_p = test_stat_p
    )
  )
  
}


#' Format observed data for Ghosh and Lin method
#'
#' @param X vector containing observed time to terminating event
#' @param delta vector that is 1 if terminating event is observed and 0 if patient 
#' is censored
#' @param Z_star FILL IN
#' @param time vector containing start times of follow-up windows chosen
#' 
#' @author Nabihah Tayob
#' 
#' @references Tayob, N. and Murray, S., 2014. Nonparametric tests of treatment 
#' effect based on combined endpoints for mortality and recurrent events. Biostatistics, 
#' 16(1), pp.73-83.

format_data_GL = function(X, delta, Z_star, time) {
  
  n_time = length(time)
  n = length(X)
  
  time_array = t(array(time, c(n_time, n)))
  delta_array = array(delta, c(n, n_time))
  dN_D_i = array(as.numeric(X == time_array & delta_array == 1), c(n, n_time))
  dN_D = apply(dN_D_i, 2, sum)
  Y_i = array(as.numeric(X >= time_array), c(n, n_time))
  Y = apply(Y_i, 2, sum)
  time_array2 <- t(array(time, c(n_time, ncol(Z_star))))

  get_dN_for_i = function(X1) {
    temp = (X1 == time_array2)
    apply(temp, 2, sum, na.rm = T)
  }
  
  dN_i = t(apply(Z_star, 1, get_dN_for_i))
  
  dN = apply(dN_i, 2, sum)
  dR_hat = dN/Y
  dch_hat = dN_D/Y
  dM_i = dN_i - Y_i*t(array(dR_hat, c(n_time, n)))
  dM_D_i = dN_D_i - Y_i*t(array(dch_hat, c(n_time, n)))
  
  #calculate ch estimate
  temp = t(array(dch_hat, c(n_time, n_time)))
  upperTriangle(temp, diag = F) = 0
  ch_hat = apply(temp, 1, sum)
  
  #calculate survival estimate
  s_hat = exp(-ch_hat)
  
  #calculate mu estimate
  dmu_hat = s_hat*dR_hat
  temp = t(array(dmu_hat, c(n_time, n_time)))
  upperTriangle(temp, diag = F) = 0
  mu_hat = apply(temp, 1, sum)  
  
  #calculate psi estimate
  temp = dM_D_i*t(array(n/Y, c(n_time, n)))
  temp2 = t(apply(temp, 1, sum_function))
  dpsi_i = dM_i*t(array(s_hat*n/Y, c(n_time, n))) - temp2*t(array(dmu_hat, c(n_time, n)))
  
  #calculate last observed event time
  tau = max(X*delta)
  
  return(
    list(
      Y = Y,
      tau = tau,
      dmu_hat = dmu_hat,
      dch_hat = dch_hat,
      dpsi_i = dpsi_i,
      dM_D_i = dM_D_i
    )
  )
  
}


#' Calculate test statistic for Ghosh and Lin method
#' 
#' @param p_GL p parameter
#' @param time time to event
#' @param data1_format Formatted data set 1
#' @param data2_format Formatted data set 2
#'
#' @return A \code{list} object which contains
#' \itemize{
#'   \item{test_stat_p - P-value for test statistic}
#' }
#' 
#' @author Nabihah Tayob
#' 
#' @references Tayob, N. and Murray, S., 2014. Nonparametric tests of treatment 
#' effect based on combined endpoints for mortality and recurrent events. Biostatistics, 
#' 16(1), pp.73-83.

GL_test_stat = function(p_GL, time, data1_format, data2_format) {
  
  Y1 = data1_format$Y
  Y2 = data2_format$Y
  tau1 = data1_format$tau
  tau2 = data2_format$tau
  dmu_hat1 = data1_format$dmu_hat
  dmu_hat2 = data2_format$dmu_hat
  dch_hat1 = data1_format$dch_hat
  dch_hat2 = data2_format$dch_hat
  dpsi_i1 = data1_format$dpsi_i
  dpsi_i2 = data2_format$dpsi_i
  dM_D_i1 = data1_format$dM_D_i
  dM_D_i2 = data2_format$dM_D_i
  
  n_time = length(time)
  
  n1 = nrow(dpsi_i1)
  n2 = nrow(dpsi_i2)
  
  tau = min(tau1, tau2)
  time_used = rep(1, length(time))
  time_used[time > tau] = 0
  
  #calculate test stats
  K_LR = Y1*Y2*(n1 + n2) / ((Y1 + Y2)*n1*n2)
  Q_LR = sum((K_LR*(dmu_hat1 - dmu_hat2))*time_used)
  Q_D = sum((K_LR*(dch_hat1 - dch_hat2))*time_used)
  Q_WCT = p_GL*Q_LR + (1 - p_GL)*Q_D
  
  #calculate variance of Q_WCT
  rho_hat1 = n1 / (n1 + n2)
  rho_hat2 = n2 / (n1 + n2)
  U_i1 = apply((p_GL*dpsi_i1 + (1 - p_GL)*dM_D_i1*t(array(n1/Y1, c(n_time, n1))))*
                 t(array(K_LR, c(n_time, n1))), 1, sum)
  U_i2 = apply((p_GL*dpsi_i2 + (1 - p_GL)*dM_D_i2*t(array(n2/Y2, c(n_time, n2))))*
                 t(array(K_LR, c(n_time, n2))), 1, sum)
  var_Q_WCT = rho_hat2*mean(U_i1^2) + rho_hat1*mean(U_i2^2)
  
  #calculate test statistic that has standard normal distribution under H0
  test_stat = sqrt(n1*n2 / (n1+n2))*Q_WCT / sqrt(var_Q_WCT)
  test_stat_p = 2*(1 - pnorm(abs(test_stat)))
  
  return(
    list(
      test_stat_p = test_stat_p
    )
  )
  
}

