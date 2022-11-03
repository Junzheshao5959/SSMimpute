## code to prepare `data_complex` dataset goes here
set.seed(123)
periods = c(500,300,300)
# parameters for c
parameters_for_c=list(baseline=extend(5,periods),
                      c=extend(0.25,periods),
                      sd=extend(sqrt(10),periods))
model_for_c=c("baseline","c")
# parameters for x
parameters_for_x=list(baseline=extend(5,periods),
                      x=extend(0.5,periods),
                      sd=extend(sqrt(10),periods))
model_for_x=c("baseline","x")
# parameters for y
lag_x_on_y=2;treatment_effects=matrix(c(-1,-0.5,-2.5,-0.5,-1,-0.5),ncol=lag_x_on_y,byrow=T)
model_for_y=c("baseline","c","x","y")
type_for_x="continuous"
parameters_for_y=list(baseline=extend(40,periods),
                      y=extend(0.5,periods),
                      x=extend_into_matrix(treatment_effects,periods),
                      c=sim_randomwalk(baseline=1.5,sd=0.05,periods = c(1100) ),
                      sd=extend(c(sqrt(.1)),periods))
colnames(parameters_for_y$x)=paste("t",0:(ncol(parameters_for_y$x)-1),sep="-")
# all parameters
sim_parameters=list(given.c=F,input.c=NULL,parameters_for_c=parameters_for_c,model_for_c=model_for_c,initial.c=NULL,
                    given.x=F,input.x=NULL,parameters_for_x=parameters_for_x,model_for_x=model_for_x,initial.x=NULL,type_for_x=type_for_x,
                    given.y=F,input.y=NULL,parameters_for_y=parameters_for_y,model_for_y=model_for_y,initial.y=NULL,lag_x_on_y=lag_x_on_y,
                    interaction_pairs=NULL,n=sum(periods))
set.seed(1)
simulation=sim_y_x_single_c(given.c=sim_parameters$given.c,input.c=sim_parameters$input.c,parameters_for_c=sim_parameters$parameters_for_c,model_for_c=sim_parameters$model_for_c,initial.c=sim_parameters$initial.c,
                            given.x=sim_parameters$given.x,input.x=sim_parameters$input.x,parameters_for_x=sim_parameters$parameters_for_x,model_for_x=sim_parameters$model_for_x,initial.x=sim_parameters$initial.x,type_for_x=sim_parameters$type_for_x,
                            given.y=sim_parameters$given.y,input.y=sim_parameters$input.y,parameters_for_y=sim_parameters$parameters_for_y,model_for_y=sim_parameters$model_for_y,initial.y=sim_parameters$initial.y,lag_x_on_y=sim_parameters$lag_x_on_y,
                            interaction_pairs=sim_parameters$interaction_pairs,n=sim_parameters$n,
                            printFlag=T)
simulation$y=simulation$y[-(1:100)]
simulation$x=simulation$x[-(1:100)]
simulation$c=simulation$c[-(1:100)]
data=data.frame(Date=Sys.Date()-days(length(simulation$y):1),y=simulation$y,x=simulation$x,c=simulation$c)
data$y_1=c(NA,data$y[-nrow(data)])
data$x_1=c(NA,data$x[-nrow(data)])
print(head(data))
data_complex=data[-1,]
usethis::use_data(data_complex, overwrite = TRUE)
