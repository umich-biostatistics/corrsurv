##############################################################
#  R package: Estimate of the tau-restricted mean survival
##############################################################

# Tayob, N. and Murray, S., 2016. Nonparametric restricted mean analysis across multiple follow-up intervals. Statistics & probability letters, 109, pp.152-158.

# Arguments for main function
# X=follow-up time of right censored data (values must be geq 0)
# delta=status indicator, 0=alive, 1=dead. (values must be 0,1)
# Tau=upper limit of integration (cannot be greater than largest follow-up time, cannot be negative)
# A=study period (default=largest follow-up time,cannot be greater, cannot be negative, should be based on rmrl plots)
# b=number of follow-up windows used (default=floor(2*(A-Tau)/Tau)+1, must be an integer, must be geq 1)
# t=start times of follow-up windows (default=seq(from=0, to=A-Tau,by=(A-Tau)/(b-1)), must be of length b if both specified, largest value cannot be greater than A-Tau, no repeats)
# var_output=c("proposed","independence","sandwich","all")
# plot=Logical argument for whether RMRL function should be plotted
  # PARAMETERS FOR RMRL PLOT
  # alpha=significance level
  # conservative_index=minimum number of events after the start of the last interval
  # k=number of points for intergration within follow-up window [0,Tau]
  # n.sim=number of samples simulated to calculate confidence bands

#Load required packages
# library(survival)
# library(MASS)

#' Tau-Restricted Mean Survival
#' 
#' Estimate the tau-restricted mean survival across multiple follow-up intervals
#' 
#' Estimate the tau-restricted mean survival as described in Tayob, N. and Murray, S., 2016. 
#' Nonparametric restricted mean analysis across multiple follow-up intervals. Statistics 
#' & probability letters, 109, pp.152-158.
#' 
#' @param X follow-up time of right censored data (values must be geq 0)
#' @param delta status indicator, 0=alive, 1=dead. (values must be 0,1)
#' @param Tau upper limit of integration (cannot be greater than largest follow-up 
#' time, cannot be negative)
#' @param A study period (default=largest follow-up time,cannot be greater, cannot 
#' be negative, should be based on rmrl plots)
#' @param b number of follow-up windows used (default=floor(2*(A-Tau)/Tau)+1, must 
#' be an integer, must be geq 1)
#' @param t start times of follow-up windows (default=seq(from=0, to=A-Tau,by=(A-Tau)/(b-1)), 
#' must be of length b if both specified, largest value cannot be greater than A-Tau, 
#' no repeats)
#' @param var_output Type of variance estimator. Options are c("proposed","independence",
#' "sandwich","all")
#' @param
#' @param
#' 
#' @return A \code{list} object which contains
#' \itemize {
#'   \item{Mean}{vector containing sample estimates of overall tau-restricted mean survival in each group}
#'   \item{Var}{vector containing empirical variance of estimates of overall tau-restricted mean survival in each group}
#'   \item{test_stat}{test statistic of two-sample test}
#'   \item{test_stat_p}{p-value of two-sample test}
#' }
#' 
#' @examples
#' data("TMdata")
#' 
#' 
#' 
#' @author Nabihah Tayob
#' 
#' @references Tayob, N. and Murray, S., 2016. Nonparametric restricted mean analysis 
#' across multiple follow-up intervals. Statistics & probability letters, 109, pp.152-158.

TM2 = function(X, delta, Tau, t, var_output = "proposed") { 
  #CHECK FOR ERRORS IN INPUT
  
  if(sum(X < 0) > 0) {
    stop("Error in follow-up times input")
  }
  
  if(sum(delta < 0 | delta > 1 | floor(delta) != delta) > 0) {
    stop("Error in failure indicator input")
  }   
  
  if(Tau > max(X)) {
    stop("Upper limit of integration greater than largest observed follow-up time")
  }
  
  if(Tau < 0) {
    stop("Invalid upper limit of integration")
  }
  
  if(max(t) > (max(X) - Tau)) {
    stop("Largest start time greater than allowable")
  }
  
  if(length(unique(t)) != length(t)) {
    stop("Repeated values in start time vector")
  }
  
  if(var_output != "proposed" & var_output != "independence" & 
     var_output != "sandwich" & var_output != "all") {
    stop("Invalid variance type")
  }
  
  args = match.call()
  args$var_output = var_output
  
  #remove any subjects with missing data
  data = data.frame(X, delta)
  data_nomiss = na.omit(data)
  X_nomiss = data_nomiss$X
  delta_nomiss = data_nomiss$delta
  n = length(X_nomiss)
  b = length(t)
  Tau = Tau
  
  #PUT DATA IN CORRECT FORMAT
  X_tk = array(NA, c(n, b))
  delta_tk = array(NA, c(n, b))
  
  for(k in 1:b) {
    X_tk[,k] = (X_nomiss - t[k])*as.numeric(X_nomiss >= t[k])
    delta_tk[,k] = delta_nomiss*as.numeric(X_nomiss >= t[k])
  }
  
  X_km = array(X_tk, c(n*b, 1))
  delta_km = array(delta_tk, c(n*b, 1))
  
  results = get_independent_var(X_km, delta_km, Tau, t, n)
  
  proposed_variance = NULL
  sandwich_variance = NULL
  independent_variance = NULL
  
  Output = "****************************************************************************************"
  Output = rbind(Output, "Nonparametric estimation of restricted mean survival across multiple follow-up intervals")
  Output = rbind(Output, "****************************************************************************************")
  Output = rbind(Output, paste("Number of patients used in the analysis is ", n))
  t_for_printing = paste(t[1])
  
  for(k in 2:b) {
    t_for_printing = paste(t_for_printing, t[k])
  }
  
  Output = rbind(Output, paste("Start time of follow-up intervals:", t_for_printing))
  Output = rbind(Output, paste("Restricted Mean Survival Estimate=", round(results$mean, 4)))
  
  if(var_output == "proposed" | var_output == "all") {
    variance_results = get_williams_var(X_km, delta_km, X_tk, delta_tk, Tau, t, n)
    Output = rbind(Output, paste("Proposed Variance=", round(variance_results$var, 4)))
    proposed_variance = variance_results$var
  }
  
  if(var_output == "sandwich" | var_output == "all") {
    variance_results = get_sandwich_var(X_km, delta_km, X_tk, delta_tk, Tau, t, n)
    Output = rbind(Output, paste("Sandwich Variance=", round(variance_results$var, 4)))
    sandwich_variance = variance_results$var
  }  
  
  if(var_output == "independence" | var_output == "all") {
    Output = rbind(Output, paste("Independent Variance=", round(results$var, 4)))
    independent_variance = results$var
  }
  
  cat(Output,sep="\n")
  
  RMRL_output=NULL
  
  # if(plot == TRUE){
  #   RMRL_output = plot_RMRL(X = X_nomiss, delta = delta_nomiss, Tau = Tau, 
  #                           alpha = alpha, conservative_index = conservative_index, 
  #                           k = k, n.sim = n.sim)
  # }
  
  list(
    Mean = results$mean,  
    proposed_variance = proposed_variance, 
    sandwich_variance = sandwich_variance, 
    independent_variance = independent_variance, 
    RMRL_output = RMRL_output,
    args = args,
    plot_args = list(
        X = X_nomiss, delta = delta_nomiss, Tau = Tau, 
        alpha = alpha
      )
    )
}


#######################
#variance functions
#######################
get_independent_var = function(X_km, delta_km, Tau, t, n) {
  b = length(t)
  km_results = summary(survfit(Surv(X_km, delta_km) ~ 1))
  T = c(0, km_results$time[km_results$time <= Tau], Tau)
  dN = c(0, km_results$n.event[km_results$time <= Tau])
  Y = c(n*b,km_results$n.risk[km_results$time <= Tau])
  M = sum(as.numeric(dN > 0))
  time_int = T[2:(M+2)] - T[1:(M+1)]
  
  CH = rep(NA, M+1)
  var_CH = rep(NA, M+1)
  
  for(m in 1:(M+1)){
    CH[m] = sum((dN/Y)[1:m])
    var_CH[m] = sum((dN/(Y^2))[1:m])
  }
  
  var_CH_array = array(NA, c(M+1, M+1))
  
  for(m in 1:(M+1)) {
    var_CH_array[m:(M+1), m:(M+1)] = var_CH[m]
  }
  
  S_hat = exp(-CH)
  mean = sum(time_int*S_hat)
  ind_var = sum(array(time_int*S_hat, c(M+1, 1)) %*% 
                 array(time_int*S_hat, c(1, M+1))*var_CH_array)
  list(
    mean = mean, 
    var = ind_var
    )
}


#proposed variance
get_williams_var = function(X_km, delta_km, X_array, delta_array, Tau, t, n) {
  b = length(t)
  km_results = summary(survfit(Surv(X_km, delta_km) ~ 1))
  T = c(0, km_results$time[km_results$time <= Tau], Tau)
  dN_old = c(0, km_results$n.event[km_results$time <= Tau])
  Y_old = c(n*b, km_results$n.risk[km_results$time <= Tau])
  M = sum(as.numeric(dN_old > 0))
  time_int = T[2:(M+2)] - T[1:(M+1)]
  CH = rep(NA, M+1)
  
  for(m in 1:(M+1)) {
    CH[m] = sum((dN_old/Y_old)[1:m])
  }

  S_hat = exp(-CH)
  
  dN_i = array(NA, c(n, M+1))
  Y_i = array(NA, c(n, M+1))
  
  for(m in 1:(M+1)) {
    dN_i[,m] = apply(array(as.numeric(round(X_array, 8) == round(T[m], 8) & delta_array == 1), 
                           c(n, b)), 1, sum)
    Y_i[,m] = apply(array(as.numeric(round(X_array, 8) >= round(T[m], 8)), 
                          c(n, b)), 1, sum)
  }
  
  dN = apply(dN_i, 2, sum)
  Y = apply(Y_i, 2, sum)
  q = dN / Y
  z_i_q = t(t(dN_i - t(q*t(Y_i)))*(1/Y))
  
  z_i_S_old = rep(0, n)
  z_i_ATau = time_int[1]*z_i_S_old
  
  for(m in 2:(M+1)) {
    z_i_S_new = exp(-q[m])*z_i_S_old + S_hat[m]*z_i_q[,m]
    z_i_ATau = z_i_ATau + time_int[m]*z_i_S_new
    z_i_S_old = z_i_S_new
  }
  
  z_bar_ATau = sum(z_i_ATau) / n
  williams_var = n*sum((z_i_ATau - z_bar_ATau)^2) / (n-1)
  
  return(
    list(
      var = williams_var
    )
  )
}


#Robust Sandwich Variance 
get_sandwich_var = function(X_km, delta_km, X_array, delta_array, Tau, t, n) {
  b = length(t)
  km_results = summary(survfit(Surv(X_km, delta_km) ~ 1))
  T = c(0, km_results$time[km_results$time <= Tau], Tau)
  dN = c(0, km_results$n.event[km_results$time <= Tau])
  Y = c(n*b, km_results$n.risk[km_results$time <= Tau])
  M = sum(as.numeric(dN > 0))
  time_int = T[2:(M+2)] - T[1:(M+1)]
  
  dN_i_Tj = array(NA, c(n, M+1))
  Y_i_Tj = array(NA, c(n, M+1))
  
  for(m in 1:(M+1)) {
    dN_i_Tj[,m] = apply(array(as.numeric(round(X_array, 8) == round(T[m], 8) & delta_array == 1), 
                              c(n, b)), 1, sum)
    Y_i_Tj[,m] = apply(array(as.numeric(round(X_array, 8) >= round(T[m], 8)), 
                             c(n, b)), 1, sum)
  }
  
  Y_array = t(array(rep(Y, n), c(M+1, n)))
  dN_array = t(array(rep(dN, n), c(M+1, n)))
  cov_matrix = (dN_i_Tj*Y_array - Y_i_Tj*dN_array) / (Y_array^2)
  
  CH = rep(NA, M+1)
  var_CH = rep(NA, M+1)
  
  for(m in 1:(M+1)) {
    CH[m] = sum((dN/Y)[1:m])
    var_CH[m] = sum((apply(array(cov_matrix[,1:m], c(n, m)), 1, sum))^2)
  }
  
  var_CH_array = array(NA, c(M+1, M+1))
  
  for(m in 1:(M+1)) {
    var_CH_array[m:(M+1), m:(M+1)] = var_CH[m]
  }
  
  S_hat = exp(-CH)
  robust_var = sum(array(time_int*S_hat, c(M+1, 1)) %*% 
                     array(time_int*S_hat, c(1, M+1))*var_CH_array)
  return(
    list(
      var = robust_var
    )
  )
  
}


#######################
#RMRL functions
#######################
plot.RMRL = function(object, alpha = 0.05, conservative_index = 25, 
                     k = 500, n.sim = 1000, ...) {
  # INPUT
  # X=time to event
  # delta=event indicator
  # Tau=length of follow-up intervals of interest
  
  # pull these in
  
  X = args$X
  delta = args$delta
  Tau = args$Tau

  n = length(X)
  Y_X = rep(NA, n)
  
  for(l in 1:n) {
    Y_X[l] = sum(as.numeric(X >= X[l]))
  }
  
  CB_tau = sort(X*delta)[n - conservative_index]
  A = max(X) #end of study time
  time_look = seq(from = 0, to = min(CB_tau, A-Tau), by = Tau/2)
  
  cat("Times at which RMRL function is evaluated:")
  cat("\n")
  cat(time_look)
  
  output = array(NA, c(length(time_look), 2))
  
  for(l in 1:length(time_look)) {
    output[l,] = expectedlife_function(X, delta, Y_X, time_look[l], Tau, k)
  } 
  
  mean = output[,1]
  Kappa = confidence_band_kappa(X, delta, Y_X, Tau, k, time_look)
  chosen_q_alpha = get_q_alpha(Kappa, alpha, n.sim)
  
  upper_CB = mean + chosen_q_alpha / sqrt(n)
  lower_CB = mean - chosen_q_alpha / sqrt(n)
  
  adj_upper_CB = apply(cbind(upper_CB, rep(Tau, length(time_look))), 1, min)
  
  par(mar = c(4,4.4,2,2) + 0.1)
  plot(time_look, mean, type = "l", ylab = "RMRL", xlab = "Time", 
       ylim = c(min(lower_CB), max(upper_CB)), xaxs = "i", yaxs = "i", lwd = 2, 
       cex.lab = 1.5, cex.axis = 1.5)
  lines(time_look, lower_CB, lty=2, lwd=2)
  lines(time_look, adj_upper_CB, lty=2, lwd=2)
  
  return(
    list(
      RMRL = mean, 
      lower_CB = lower_CB, 
      upper_CB = upper_CB
    ) 
  )
}


prosper_function = function(X, delta, Y_X, t_in, a_in) {
  sum1 = sum(as.numeric(X > t_in & X <= (t_in + a_in))*delta/Y_X)
  exp(-sum1)
}


expectedlife_function = function(X, delta, Y_X, t_in, a, k) {
  n = length(X)
  #construct parameters for trapazoidal method
  a_i = seq(from = 0, to = a, length.out = k+1)
  d_i = c(a_i, a) - c(0, a_i)
  b_i = (d_i[1:(k+1)] + d_i[2:(k+2)]) / 2
  
  #get values for prosper function
  prosper_a_i = rep(NA, k+1)
  
  for(i in 1:(k+1)) {
    prosper_a_i[i] = prosper_function(X, delta, Y_X, t_in, a_i[i])
  }
  
  estimate = sum(b_i*prosper_a_i)
  
  array1 = array(as.numeric(X > t_in)*delta/(Y_X^2), 
                 c(n, k+1))*array(as.numeric(array(X, c(n, k+1)) <= 
                                               t(array(t_in + a_i, c(k+1, n)))), c(n, k+1))
  
  array2 = array(as.numeric(array(X, c(n, k+1)) <= 
                              t(array(t_in + a_i, c(k+1, n)))), c(n, k+1))
  inner_array = t(array1) %*% array2
  
  array3 = array(prosper_a_i*b_i, c(k+1, k+1))
  array4 = t(array(prosper_a_i*b_i, c(k+1, k+1)))
  array5 = array3*array4*inner_array
  variance = sum(array5)
  
  #estimate
  return(
    c(estimate,variance)
  )
}


confidence_band_kappa = function(X, delta, Y_X, a, k, t_k) { 
  n = length(X)
  a_i = seq(from = 0, to = a, length.out = k+1)
  d_i = c(a_i, a) - c(0, a_i)
  b_i = (d_i[1:(k+1)] + d_i[2:(k+2)]) / 2
  
  n_star = length(t_k)
  
  #get values for prosper function for each time
  prosper_a_i = array(NA, c(n_star, k+1))
  
  for(i in 1:n_star) {
    for(j in 1:(k+1)) {
      prosper_a_i[i,j] = prosper_function(X, delta, Y_X, t_k[i], a_i[j])
    }
  }
  
  kappa = array(NA, c(n_star, n_star))
  
  for(i_star in 1:n_star) {
    for(j_star in 1:n_star) {
      array1 = array(as.numeric(X > t_k[i_star])*as.numeric(X > t_k[j_star])*delta/(Y_X^2), c(n, k+1)) * 
        array(as.numeric(array(X, c(n, k+1)) <= t(array(t_k[i_star] + a_i, c(k+1, n)))), c(n, k+1))
      
      array2 = array(as.numeric(array(X, c(n, k+1)) <= t(array(t_k[j_star] + a_i, c(k+1, n)))), c(n, k+1))
      inner_array = t(array1) %*% array2
      
      array3 = array(prosper_a_i[i_star,]*b_i, c(k+1, k+1))
      array4 = t(array(prosper_a_i[j_star,]*b_i, c(k+1, k+1)))
      array5 = array3*array4*inner_array
      kappa[i_star,j_star] = n*sum(array5)
    }
  }
  return(kappa)
}

get_q_alpha = function(kappa, alpha, n_sim) {
  N = nrow(kappa) 
  simulated_B = mvrnorm(n_sim, rep(0, N), kappa)
  sim_datasets = apply(abs(simulated_B), 1, max)
  q_alpha = sort(sim_datasets)[n_sim*(1 - alpha)]
  return(q_alpha)
}

#####################
# Example - to add to documentation!
#####################

# #read in example data
# data=read.csv("example_data.csv")
# max(data$X)
# output=TM2(X=data$X,delta=data$delta,Tau=12,t=seq(from=0,to=24,by=6),var_output="all",plot=TRUE,conservative_index=10)
# 
# output=TM2(X=data$X,delta=data$delta,Tau=12,t=seq(from=0,to=24,by=6))
# # ****************************************************************************************
# #   Nonparametric estimation of restricted mean survival across multiple follow-up intervals
# # ****************************************************************************************
# #   Number of patients used in the analysis is  100
# # Start time of follow-up intervals: 0 6 12 18 24
# # Restricted Mean Survival Estimate= 10.9551
# # Proposed Variance= 0.0277
# 
# output=TM2(X=data$X,delta=data$delta,Tau=12,t=seq(from=0,to=24,by=6),var_output="independence")
# # ****************************************************************************************
# #   Nonparametric estimation of restricted mean survival across multiple follow-up intervals
# # ****************************************************************************************
# #   Number of patients used in the analysis is  100
# # Start time of follow-up intervals: 0 6 12 18 24
# # Restricted Mean Survival Estimate= 10.9551
# # Independent Variance= 0.0186
# 
# output=TM2(X=data$X,delta=data$delta,Tau=12,t=seq(from=0,to=24,by=6),var_output="sandwich")
# # ****************************************************************************************
# #   Nonparametric estimation of restricted mean survival across multiple follow-up intervals
# # ****************************************************************************************
# #   Number of patients used in the analysis is  100
# # Start time of follow-up intervals: 0 6 12 18 24
# # Restricted Mean Survival Estimate= 10.9551
# # Sandwich Variance= 0.022
