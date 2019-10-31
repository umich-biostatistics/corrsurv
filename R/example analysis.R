#FIRST RUN ALL CODE IN FILE: functions.R

#LOAD DATA
data=read.csv("example_data.csv")

#RUN ANALYSIS
N=nrow(data) #Number of subjects
X=data$X
delta=data$delta
table(delta)
Z=array(cbind(data$Z1,data$Z2,data$Z3,data$Z4,data$Z5,data$Z6),c(N,6))
A=max(X) #length of study
Treatment=as.numeric(data$Group==1) #1 if Case and 0 if Control
table(Treatment,delta)

two_sample_recurrent_events_test(X=X,delta=delta,Z=Z,Group=Treatment,Tau=A/2,t=seq(from=0,to=A-A/2,by=A/4))
# ************Two-Sample Test for combined end-point across multiple follow-up windows************
#   Reference Paper: Nonparametric Tests of Treatment Effect for a Recurrent Event process that Terminates- N Tayob and S Murray
# 
# Group definitions:
#   Group 1: Group= 0
# Group 2: Group= 1
# 
# Sample Size: 
#   Group 1:  100
# Group 2:  100
# 
# Follow-up windows:
#   [ 0 , 18 ]       [ 9 , 27 ]       [ 18 , 36 ]      
# 
# Sample Estimates: Average restricted mean survival across follow-up windows
# Group 1:  8.18453
# Group 2:  10.26689
# 
# Alternative hypothesis: True difference in the average restricted mean survival across follow-up windows is not equal to 0
# Test-statistic= 3.18359 , p-value = 0.00145$Mean
# [1]  8.184535 10.266895
# 
# $Var
# [1] 0.1937549 0.2340821
# 
# $test_stat
# [1] 3.183586
# 
# $test_stat_p
# [1] 0.001454627

two_sample_recurrent_events_test(X=X,delta=delta,Z=Z,Group=Treatment,Tau=A/2,t=seq(from=0,to=A-A/2,by=A/4),method="area",plot=TRUE)
# ************Two-Sample Test for combined end-point across multiple follow-up windows************
#   Reference Paper: Nonparametric Tests of Treatment Effect for a Recurrent Event process that Terminates- N Tayob and S Murray
# 
# Group definitions:
# Group 1: Group= 0
# Group 2: Group= 1
# 
# Sample Size: 
# Group 1:  100
# Group 2:  100
# 
# Follow-up windows:
#   [ 0 , 18 ]       [ 9 , 27 ]       [ 18 , 36 ]      
# 
# Sample Estimates: Area under RMRL function evaluated for each follow-up window
# Group 1:  147.8743
# Group 2:  185.3959
# 
# Alternative hypothesis: True difference in the area under the RMRL functions is not equal to 0
# Test-statistic= 3.05018 , p-value = 0.00229$Mean
# [1] 147.8743 185.3959
# 
# $Var
# [1] 68.40778 82.91802
# 
# $test_stat
# [1] 3.050175
# 
# $test_stat_p
# [1] 0.002287079


#coxph test
placebo_data=format_data_ourmethod(X=X[Treatment==0],delta=delta[Treatment==0],Z=Z[Treatment==0,],t=c(0))
active_data=format_data_ourmethod(X=X[Treatment==1],delta=delta[Treatment==1],Z=Z[Treatment==1,],t=c(0))
Time=c(placebo_data$X_tk,active_data$X_tk)
Delta=c(placebo_data$delta_tk,active_data$delta_tk)
X_treat=c(rep(0,length(placebo_data$X_tk)),rep(1,length(active_data$X_tk)))
summary(coxph(Surv(Time,Delta)~X_treat))
plot(survfit(Surv(Time,Delta)~X_treat))

#Anderson and Gill Method
placebo_data_AG=format_data_AG(X=X[Treatment==0],delta=delta[Treatment==0],Z_star=Z[Treatment==0,])
active_data_AG=format_data_AG(X=X[Treatment==1],delta=delta[Treatment==1],Z_star=Z[Treatment==1,])
placebo_data_AG_2=subset(data.frame(placebo_data_AG),placebo_data_AG$T_stop>placebo_data_AG$T_start)
active_data_AG_2=subset(data.frame(active_data_AG),active_data_AG$T_stop>active_data_AG$T_start)
AG_test_stat(placebo_data_AG_2,active_data_AG_2)

#Ghosh and Lin
time1=unique(sort(cbind(Z,X)))
placebo_data_GL=format_data_GL(X=X[Treatment==0],delta=delta[Treatment==0],Z_star=Z[Treatment==0,],time1)
active_data_GL=format_data_GL(X=X[Treatment==1],delta=delta[Treatment==1],Z_star=Z[Treatment==1,],time1)
GL_test_stat(0.5,time1,placebo_data_GL,active_data_GL)$test_stat_p

