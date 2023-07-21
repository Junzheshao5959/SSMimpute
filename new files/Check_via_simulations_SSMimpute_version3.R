library(MASS)
library(crayon) # for cat colors
library(lubridate) # for dates
library(mice) # for multiple imputation
library(forecast) # for most the time series functions: auto.arima
library(changepoint) # for finding change points
library(dlm) # for statespace model
library(dplyr) # for mutate
source("new files/Summary 2_helper functions.R")
source("new files/helper_multipleimputation.R")
source("new files/helper_SSM_version3.R")
source("new files/helper_SSMimpute_version3.R")

# Example 1: (fixed effects for all coefficients)
#            ss_param=list(inits=c(log(1)),m0=c(40,0.5,-1.5,-0.5,-1),C0=diag(rep(10^3),5),
#                          AR1_coeffi=NULL,rw_coeffi=NULL,v_cp_param=NULL,w_cp_param=NULL,max_iteration=50)
# Checked on 05/19/23
{
  set.seed(1)
  periods = 1100
  # parameters for c
  parameters_for_c=list(baseline=extend(5,periods),c=extend(0.25,periods),sd=extend(sqrt(1),periods))
  model_for_c=c("baseline","c")
  # parameters for x
  parameters_for_x=list(baseline=extend(5,periods),x=extend(0.5,periods),sd=extend(sqrt(1),periods))
  model_for_x=c("baseline","x")
  # parameters for y
  lag_x_on_y=2;lag_y_on_y=1;lag_c_on_y=1;
  treatment_effects=matrix(c(-1.5,-0.5),ncol=lag_x_on_y,byrow=T)
  model_for_y=c("baseline","c","x","y")
  type_for_x="continuous"
  parameters_for_y=list(baseline=as.matrix(extend(40,periods)),
                        y=as.matrix(extend(0.5,periods)),
                        x=extend_into_matrix(treatment_effects,periods),
                        c=as.matrix(extend(-1,periods)),
                        sd=extend(c(sqrt(1)),periods))
  colnames(parameters_for_y$x)=paste("t",0:(ncol(parameters_for_y$x)-1),sep="-")
  colnames(parameters_for_y$y)=paste("t",1:ncol(parameters_for_y$y),sep="-")
  colnames(parameters_for_y$c)=paste("t",0:(ncol(parameters_for_y$c)-1),sep="-")
  # all parameters
  sim_parameters=list(given.c=F,input.c=NULL,parameters_for_c=parameters_for_c,model_for_c=model_for_c,initial.c=NULL,
                      given.x=F,input.x=NULL,parameters_for_x=parameters_for_x,model_for_x=model_for_x,initial.x=NULL,type_for_x=type_for_x,
                      given.y=F,input.y=NULL,parameters_for_y=parameters_for_y,model_for_y=model_for_y,initial.y=NULL,
                      lag_y_on_y=lag_y_on_y,lag_x_on_y=lag_x_on_y,lag_c_on_y=lag_c_on_y,
                      interaction_pairs=NULL,n=sum(periods))
  simulation=sim_y_x_c_univariate(given.c=sim_parameters$given.c,input.c=sim_parameters$input.c,parameters_for_c=sim_parameters$parameters_for_c,model_for_c=sim_parameters$model_for_c,initial.c=sim_parameters$initial.c,
                                  given.x=sim_parameters$given.x,input.x=sim_parameters$input.x,parameters_for_x=sim_parameters$parameters_for_x,model_for_x=sim_parameters$model_for_x,initial.x=sim_parameters$initial.x,type_for_x=sim_parameters$type_for_x,
                                  given.y=sim_parameters$given.y,input.y=sim_parameters$input.y,parameters_for_y=sim_parameters$parameters_for_y,model_for_y=sim_parameters$model_for_y,initial.y=sim_parameters$initial.y,
                                  lag_y_on_y=sim_parameters$lag_y_on_y,lag_x_on_y=sim_parameters$lag_x_on_y,lag_c_on_y=sim_parameters$lag_c_on_y,
                                  interaction_pairs=sim_parameters$interaction_pairs,n=sim_parameters$n,
                                  printFlag=F)
  simulation$y=simulation$y[-(1:100)]
  simulation$x=simulation$x[-(1:100)]
  simulation$c=simulation$c[-(1:100)]
  data=data.frame(Date=Sys.Date()-days(length(simulation$y):1),y=simulation$y,x=simulation$x,c=simulation$c)
  rm(lag_c_on_y,lag_x_on_y,lag_y_on_y,model_for_c,model_for_x,model_for_y,periods,type_for_x)
  rm(treatment_effects,parameters_for_c,parameters_for_x,parameters_for_y)
  rm(sim_parameters,simulation)
  set.seed(1)
  missing_index=generate.missing_index(type = "MCAR", n=length(data$y), param = list(p=0.5))$missing_index
  data_mis_y=insert.missingness(data,variables="y",missing_index=missing_index)
  add_variable_param=list(list(lagged_param=list(variables=c("y","x","c"),param=c(1,1,1))))
  data_ssm_mis_y=add_variables_procedures(data_mis_y,param=add_variable_param,printFlag=T)
  data_ssm_mis_y=data_ssm_mis_y[-1,]
}
ss_param=list(inits=c(log(1)),m0=c(40,0.5,-1.5,-0.5,-1),C0=diag(rep(10^3),5),
              AR1_coeffi=NULL,rw_coeffi=NULL,v_cp_param=NULL,w_cp_param=NULL,max_iteration=50)
result_new=run.SSMimpute(data_ss_ori=data_ssm_mis_y,formula="y~y_1+x+x_1+c",ss_param_temp=ss_param,
                         initial_imputation_option="StructTS",estimate_convergence_cri=0.01,lik_convergence_cri=0.01,stepsize_for_newpart=1/3,max_iteration=100,
                         cpt_learning_param=list(cpt_method="mean",burnin=1/10,mergeband=20,convergence_cri=15),
                         cpt_initial_guess_option="ignore", cpt_merge_option = "unanimous",
                         dlm_imputation_option="smooth", dlm_option="smooth",
                         dlm_cpt_learning_option="smooth", bandwidth=20, cpt_V=5, m=5,seed=1,printFlag=T)
result_new$result_convergence
result_old=run.SSMimpute_unanimous_cpts(data_ss_ori=data_ssm_mis_y,formula=c("y_1","x","x_1","c"),ss_param_temp=ss_param,
                                        initial_imputation_option="StructTS",estimate_convergence_cri=0.01,lik_convergence_cri=0.01,stepsize_for_newpart=1/3,max_iteration=100,
                                        cpt_learning_param=list(cpt_method="mean",burnin=1/10,mergeband=20,convergence_cri=15),
                                        cpt_initial_guess_option="ignore",
                                        dlm_option="smooth",m=5,seed=1,printFlag=T)
result_old$result_convergence

# Example 6:   (periodic coefficient for X of 3 periods with unknown change points)
#              ss_param=list(inits=c(log(1)),m0=c(40,0.5,-1/-2,-0.5,1),C0=diag(rep(10^3),5),
#                            AR1_coeffi=NULL,rw_coeffi=NULL,v_cp_param=NULL,
#                            w_cp_param=list(list(variable="x",segments=3)),
#                            max_iteration=50)
# Example 6.1: (periodic coefficient for X of 3 periods with known change points)
#              ss_param=list(inits=c(log(1)),m0=c(40,0.5,-1/-2,-0.5,1),C0=diag(rep(10^3),5),
#                            AR1_coeffi=NULL,rw_coeffi=NULL, v_cp_param=NULL,
#                            w_cp_param=list(list(variable="x","segments"=3,changepoints=c(400,700),fixed_cpts=T)),
#                            max_iteration=50)
# Checked on 05/19/23
{
  set.seed(1)
  periods = c(400,300,300)
  # parameters for c
  parameters_for_c=list(baseline=extend(5,periods),c=extend(0.25,periods),sd=extend(sqrt(1),periods))
  model_for_c=c("baseline","c")
  # parameters for x
  parameters_for_x=list(baseline=extend(5,periods),x=extend(0.5,periods),sd=extend(sqrt(1),periods))
  model_for_x=c("baseline","x")
  # parameters for y
  lag_x_on_y=2;lag_y_on_y=1;lag_c_on_y=1;
  treatment_effects=matrix(c(-1,-0.5,-2,-0.5,-1,-0.5),ncol=lag_x_on_y,byrow=T)
  model_for_y=c("baseline","c","x","y")
  type_for_x="continuous"
  parameters_for_y=list(baseline=as.matrix(extend(40,periods)),
                        y=as.matrix(extend(0.5,periods)),
                        x=extend_into_matrix(treatment_effects,periods),
                        c=as.matrix(extend(-1,periods)),
                        sd=extend(c(sqrt(1)),periods))
  colnames(parameters_for_y$x)=paste("t",0:(ncol(parameters_for_y$x)-1),sep="-")
  colnames(parameters_for_y$y)=paste("t",1:ncol(parameters_for_y$y),sep="-")
  colnames(parameters_for_y$c)=paste("t",0:(ncol(parameters_for_y$c)-1),sep="-")
  # all parameters
  sim_parameters=list(given.c=F,input.c=NULL,parameters_for_c=parameters_for_c,model_for_c=model_for_c,initial.c=NULL,
                      given.x=F,input.x=NULL,parameters_for_x=parameters_for_x,model_for_x=model_for_x,initial.x=NULL,type_for_x=type_for_x,
                      given.y=F,input.y=NULL,parameters_for_y=parameters_for_y,model_for_y=model_for_y,initial.y=NULL,
                      lag_y_on_y=lag_y_on_y,lag_x_on_y=lag_x_on_y,lag_c_on_y=lag_c_on_y,
                      interaction_pairs=NULL,n=sum(periods))
  simulation=sim_y_x_c_univariate(given.c=sim_parameters$given.c,input.c=sim_parameters$input.c,parameters_for_c=sim_parameters$parameters_for_c,model_for_c=sim_parameters$model_for_c,initial.c=sim_parameters$initial.c,
                                  given.x=sim_parameters$given.x,input.x=sim_parameters$input.x,parameters_for_x=sim_parameters$parameters_for_x,model_for_x=sim_parameters$model_for_x,initial.x=sim_parameters$initial.x,type_for_x=sim_parameters$type_for_x,
                                  given.y=sim_parameters$given.y,input.y=sim_parameters$input.y,parameters_for_y=sim_parameters$parameters_for_y,model_for_y=sim_parameters$model_for_y,initial.y=sim_parameters$initial.y,
                                  lag_y_on_y=sim_parameters$lag_y_on_y,lag_x_on_y=sim_parameters$lag_x_on_y,lag_c_on_y=sim_parameters$lag_c_on_y,
                                  interaction_pairs=sim_parameters$interaction_pairs,n=sim_parameters$n,
                                  printFlag=F)
  simulation$y=simulation$y[-(1:100)]
  simulation$x=simulation$x[-(1:100)]
  simulation$c=simulation$c[-(1:100)]
  data=data.frame(Date=Sys.Date()-days(length(simulation$y):1),y=simulation$y,x=simulation$x,c=simulation$c)
  rm(lag_c_on_y,lag_x_on_y,lag_y_on_y,model_for_c,model_for_x,model_for_y,periods,type_for_x)
  rm(treatment_effects,parameters_for_c,parameters_for_x,parameters_for_y)
  rm(sim_parameters,simulation)
  plot(data$y)
  plot(data$x)
  plot(data$c)
  set.seed(1)
  missing_index=generate.missing_index(type = "MCAR", n=length(data$y), param = list(p=0.5))$missing_index
  data_mis_y=insert.missingness(data,variables="y",missing_index=missing_index)
  add_variable_param=list(list(lagged_param=list(variables=c("y","x","c"),param=c(1,1,1))))
  data_ssm_mis_y=add_variables_procedures(data_mis_y,param=add_variable_param,printFlag=T)
  data_ssm_mis_y=data_ssm_mis_y[-1,]
}
ss_param=list(inits=c(log(1)),m0=c(40,0.5,-1,-0.5,-1),C0=diag(rep(10^3),5),
              AR1_coeffi=NULL,rw_coeffi=NULL,v_cp_param=NULL,
              w_cp_param=list(list(variable="x",segments=3)),max_iteration=50)
ss_param=list(inits=c(log(1)),m0=c(40,0.5,-1,-0.5,-1),C0=diag(rep(10^3),5),
              AR1_coeffi=NULL,rw_coeffi=NULL,v_cp_param=NULL,
              w_cp_param=list(list(variable="x",segments=3,changepoints=c(300,600),fixed_cpts=F)),max_iteration=50)
result_new=run.SSMimpute(data_ss_ori=data_ssm_mis_y,formula="y~y_1+x+x_1+c",ss_param_temp=ss_param,
                         initial_imputation_option="StructTS",estimate_convergence_cri=0.01,lik_convergence_cri=0.01,stepsize_for_newpart=1/3,max_iteration=100,
                         cpt_learning_param=list(cpt_method="mean",burnin=1/10,mergeband=20,convergence_cri=15),
                         cpt_initial_guess_option="ignore", cpt_merge_option = "unanimous",
                         dlm_imputation_option="smooth", dlm_option="smooth",
                         dlm_cpt_learning_option="smooth", bandwidth=10, cpt_V=10, m=5,seed=1,printFlag=T)
result_new$result_convergence
result_new$estimated_cpts
result_old=run.SSMimpute_unanimous_cpts(data_ss_ori=data_ssm_mis_y,formula=c("y_1","x","x_1","c"),ss_param_temp=ss_param,
                                        initial_imputation_option="StructTS",estimate_convergence_cri=0.01,lik_convergence_cri=0.01,stepsize_for_newpart=1/3,max_iteration=100,
                                        cpt_learning_param=list(cpt_method="mean",burnin=1/10,mergeband=20,convergence_cri=15),
                                        cpt_initial_guess_option="ignore",
                                        dlm_option="smooth",m=5,seed=1,printFlag=T)
result_old$result_convergence
result_old$estimated_cpts

# Example 6.2: (periodic coefficient for X and X_1 of 3 uncertain known perios)
#              ss_param=list(inits=c(log(1)),m0=c(40,0.5,-1/-2,-0.5/0,1),C0=diag(rep(10^3),5),
#                            AR1_coeffi=NULL,rw_coeffi=NULL,
#                            v_cp_param=NULL,
#                            w_cp_param=list(list(variable="x","segments"=3,changepoints=c(300,600),fixed_cpts=F),
#                                            list(variable="x_1","segments"=3,changepoints=c(100,600),fixed_cpts=F)),
#                            max_iteration=max_iteration)
{
  set.seed(1)
  periods = c(200,200,300,300)
  # parameters for c
  parameters_for_c=list(baseline=extend(5,periods),c=extend(0.25,periods),sd=extend(sqrt(1),periods))
  model_for_c=c("baseline","c")
  # parameters for x
  parameters_for_x=list(baseline=extend(5,periods),x=extend(0.5,periods),sd=extend(sqrt(1),periods))
  model_for_x=c("baseline","x")
  # parameters for y
  lag_x_on_y=2;lag_y_on_y=1;lag_c_on_y=1;
  treatment_effects=matrix(c(-1,-0.5,-1,0,-2,0,-1,-0.5),ncol=lag_x_on_y,byrow=T)
  model_for_y=c("baseline","c","x","y")
  type_for_x="continuous"
  parameters_for_y=list(baseline=as.matrix(extend(40,periods)),
                        y=as.matrix(extend(0.5,periods)),
                        x=extend_into_matrix(treatment_effects,periods),
                        c=as.matrix(extend(-1,periods)),
                        sd=extend(c(sqrt(1)),periods))
  colnames(parameters_for_y$x)=paste("t",0:(ncol(parameters_for_y$x)-1),sep="-")
  colnames(parameters_for_y$y)=paste("t",1:ncol(parameters_for_y$y),sep="-")
  colnames(parameters_for_y$c)=paste("t",0:(ncol(parameters_for_y$c)-1),sep="-")
  # all parameters
  sim_parameters=list(given.c=F,input.c=NULL,parameters_for_c=parameters_for_c,model_for_c=model_for_c,initial.c=NULL,
                      given.x=F,input.x=NULL,parameters_for_x=parameters_for_x,model_for_x=model_for_x,initial.x=NULL,type_for_x=type_for_x,
                      given.y=F,input.y=NULL,parameters_for_y=parameters_for_y,model_for_y=model_for_y,initial.y=NULL,
                      lag_y_on_y=lag_y_on_y,lag_x_on_y=lag_x_on_y,lag_c_on_y=lag_c_on_y,
                      interaction_pairs=NULL,n=sum(periods))
  simulation=sim_y_x_c_univariate(given.c=sim_parameters$given.c,input.c=sim_parameters$input.c,parameters_for_c=sim_parameters$parameters_for_c,model_for_c=sim_parameters$model_for_c,initial.c=sim_parameters$initial.c,
                                  given.x=sim_parameters$given.x,input.x=sim_parameters$input.x,parameters_for_x=sim_parameters$parameters_for_x,model_for_x=sim_parameters$model_for_x,initial.x=sim_parameters$initial.x,type_for_x=sim_parameters$type_for_x,
                                  given.y=sim_parameters$given.y,input.y=sim_parameters$input.y,parameters_for_y=sim_parameters$parameters_for_y,model_for_y=sim_parameters$model_for_y,initial.y=sim_parameters$initial.y,
                                  lag_y_on_y=sim_parameters$lag_y_on_y,lag_x_on_y=sim_parameters$lag_x_on_y,lag_c_on_y=sim_parameters$lag_c_on_y,
                                  interaction_pairs=sim_parameters$interaction_pairs,n=sim_parameters$n,
                                  printFlag=F)
  simulation$y=simulation$y[-(1:100)]
  simulation$x=simulation$x[-(1:100)]
  simulation$c=simulation$c[-(1:100)]
  data=data.frame(Date=Sys.Date()-days(length(simulation$y):1),y=simulation$y,x=simulation$x,c=simulation$c)
  rm(lag_c_on_y,lag_x_on_y,lag_y_on_y,model_for_c,model_for_x,model_for_y,periods,type_for_x)
  rm(treatment_effects,parameters_for_c,parameters_for_x,parameters_for_y)
  rm(sim_parameters,simulation)
  plot(data$y)
  plot(data$x)
  plot(data$c)
  set.seed(1)
  missing_index=generate.missing_index(type = "MCAR", n=length(data$y), param = list(p=0.5))$missing_index
  data_mis_y=insert.missingness(data,variables="y",missing_index=missing_index)
  add_variable_param=list(list(lagged_param=list(variables=c("y","x","c"),param=c(1,1,1))))
  data_ssm_mis_y=add_variables_procedures(data_mis_y,param=add_variable_param,printFlag=T)
  data_ssm_mis_y=data_ssm_mis_y[-1,]
}
ss_param=list(inits=c(log(1)),m0=c(40,0.5,-1,-0.5,-1),C0=diag(rep(10^3),5),
              AR1_coeffi=NULL,rw_coeffi=NULL,v_cp_param=NULL,
              w_cp_param=list(list(variable="x",segments=3),
                              list(variable="x_1",segments=3)),max_iteration=50)
ss_param=list(inits=c(log(1)),m0=c(40,0.5,-1,-0.5,-1),C0=diag(rep(10^3),5),
              AR1_coeffi=NULL,rw_coeffi=NULL,v_cp_param=NULL,
              w_cp_param=list(list(variable="x",segments=3,changepoints=c(300,600),fixed_cpts=F),
                              list(variable="x_1",segments=3,changepoints=c(100,600),fixed_cpts=F)),max_iteration=50)
result_new=run.SSMimpute(data_ss_ori=data_ssm_mis_y,formula="y~y_1+x+x_1+c",ss_param_temp=ss_param,
                         initial_imputation_option="StructTS",estimate_convergence_cri=0.01,lik_convergence_cri=0.01,stepsize_for_newpart=1/3,max_iteration=100,
                         cpt_learning_param=list(cpt_method="mean",burnin=1/10,mergeband=20,convergence_cri=15),
                         cpt_initial_guess_option="ignore", cpt_merge_option = "separate",
                         dlm_imputation_option="smooth", dlm_option="smooth",
                         dlm_cpt_learning_option="smooth", bandwidth=20, cpt_V=5, m=5,seed=1,printFlag=T)
result_new$result_convergence
result_new$estimated_cpts
result_old=run.SSMimpute_partial_converged_cpts(data_ss_ori=data_ssm_mis_y,formula=c("y_1","x","x_1","c"),ss_param_temp=ss_param,
                                                initial_imputation_option="StructTS",estimate_convergence_cri=0.01,lik_convergence_cri=0.01,stepsize_for_newpart=1/3,max_iteration=100,
                                                cpt_learning_param=list(cpt_method="mean",burnin=1/10,mergeband=20,convergence_cri=15),
                                                cpt_initial_guess_option="ignore",dlm_imputation_option="smooth",
                                                dlm_cpt_learning_option="smooth",dlm_option="smooth",
                                                bandwidth=20,cpt_V=5, m=5,seed=1,printFlag=T)
result_old$result_convergence
result_old$estimated_cpts



