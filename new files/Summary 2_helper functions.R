# File: 2021.05.21 helper functions summaries
# Date: 2021.05.21
# Content: contain all function in "/Users/xiaoxuancai/Dropbox (Personal)/MHealthPsychSummerProject2020/Xiaoxuan_Cai/Rcode/discarded codes/helper_summary.R"
#          contain all function in "/Users/xiaoxuancai/Documents/analysis/discard_Rfunctions/helper_traditional_timeseries_analysis.R"
###################################################################
##        simulate intertwined time series of y, x and c         ##
###################################################################
# Goal: simulate time series data, when y, x and c are all entangled together with possible one interaction term
#       y ~ baseline + (y_1) + (x) + (x_1,...x_k) + (c) + noise
#       x ~ baseline + (y_1) + (x_1) + (c) + noise
#       c ~ baseline + (y_1) + (x_1) + (c_1) + noise
#       used for "20210122_simulation_all.R" ("/Users/xiaoxuancai/Dropbox (Personal)/MHealthPsychSummerProject2020/Xiaoxuan_Cai/Rcode")
#                "20210122_exploration for the delay period of tx effect(rewrite 20201015 version).R" ("/Users/xiaoxuancai/Dropbox (Personal)/MHealthPsychSummerProject2020/Xiaoxuan_Cai/20201014_Exploration for delay and duration of tx effect")
{
  # [extend] extend a single number, or a vector to a length of the same length of the time points; unchanged is given the full length vector
  extend=function(parameter,periods){
    if(length(parameter)==1){
      result=rep(parameter,sum(periods))
    }else if(length(parameter)==length(periods)){
      result=rep(parameter,periods)
    }else if(length(parameter)==sum(periods)){
      return=parameter
    }else{
      stop(red("The format of the parameter is not scalar, or of same length of #periods, or of same length of timline!\n"))
    }
    return(result)
  } # transform input to a vector of the same length of the timeline
  # [extend_into_matrix] extend one line or a small matrix of #periods of lines into a full-length matrix, which is of the same length of the time points
  extend_into_matrix=function(parameters,periods){
    if(is.vector(parameters)==T){
      result=matrix(rep(parameters,sum(periods))
                    ,ncol=length(parameters),nrow=sum(periods),
                    byrow=T)
    }else if(is.matrix(parameters)==T){
      if(length(periods)==nrow(parameters)){
        result=matrix(NA,ncol=ncol(parameters),nrow=sum(periods))
        for(j in 1:ncol(parameters)){
          result[,j]=extend(parameters[,j],periods)
        }
      }else{
        stop("The periods in parameter and ``periods'' should agree!\n")
      }
    }else{
      stop("The input should be either a vector or a matrix!\n")
    }  
    return(result)
  }
  # generate time series with tangled y, x and c
  #   current y ~ previous y + current x + multiple past x + current c + interaction
  #   current x ~ previous y + previous x + current c
  #   current c ~ previous y + previous x + previous c
  sim_y_x_single_c=function(given.c=NULL,input.c=NULL,parameters_for_c=NULL,model_for_c=NULL,initial.c=NULL,
                            given.x=F,input.x=NULL,parameters_for_x=NULL,model_for_x=NULL,initial.x=NULL,type_for_x,
                            given.y=F,input.y=NULL,parameters_for_y=NULL,model_for_y=NULL,initial.y=NULL,lag_x_on_y=NULL,interaction_pairs=NULL,
                            n,printFlag=F){
    # Goal: generate time series with tangled y, x and c
    #   current y ~ previous y + current x + multiple past x + current c + interaction
    #   current x ~ previous y + previous x + current c
    #   current c ~ previous y + previous x + previous c
    
    # Check covariates: if need to generate -> GENERATE.C=T
    #                   if given -> GENERATE.C=F
    #                   if error -> stop the process
    if(given.c==T){ # if c is given, take it
      if(!is.null(input.c) & length(input.c) == n & is.null(parameters_for_c) & is.null(initial.c) & is.null(model_for_c)){
        # specifing "given.c==T" means that covariate will be given, rather than being generated
        #   check (i) input for C exists
        #         (ii) no speficiation of stuff used for generating covariates: "parameters_for_c","model_for_c", and "initial.c"
        GENERATE.C=F
        c=input.c
        if(printFlag){cat("The covariates are given!\n")}
      }else{
        if(is.null(input.c)) stop("Covariates information is needed, when we specify covariates as given!\n")
        if(!is.null(parameters_for_c)) stop("Coefficients to generate covariates are not needed, as covariates are given!\n")
        if(!is.null(initial.c)) stop("Initial covariate value is not needed, as covariates are given!\n")
        if(!is.null(model_for_c)) stop("Models for covariate generation is not needed, as covariates are given!\n")
        if(length(input.c) != n) stop("The length of input covariates is wrong!\n")
      }
    }else if(given.c==F){ # if c is not given, generate C
      # specifing "given.c==F" means that covariate need to be generated
      #   check (i) no input for C
      #         (ii) speficiation of stuff used for generating covariates: "parameters_for_c","model_for_c", "initial.c"(optional)
      if(is.null(input.c) & !is.null(parameters_for_c) & !is.null(model_for_c)){
        if(printFlag){cat("The generation of covariates depend on: ",paste(paste("previous",model_for_c),collapse=", "),".\n",sep="")}
        # check (iii) if the dimention of input variables agrees with that of the coefficients
        #       (iv)  baseline coefficient and sd of noise must exist 
        #       (v)   auto-correlation coefficient is optional
        if(length(model_for_c) == (length(parameters_for_c)-1) & 
           "sd" %in% names(parameters_for_c) &
           "baseline" %in% names(parameters_for_c)
        ){
          GENERATE.C=T
          c=rep(NA,n)
          noise_c=rnorm(n,sd=parameters_for_c$sd)
        }else if(length(model_for_c) != (length(parameters_for_c)-1)){
          stop("The dimention of coefficients for covariate generation doesn't agreen with its model!\n")
        }else{
          stop("Baseline or sd of noise is missing!\n")
        }
      }else{
        if(!is.null(input.c)) stop("Covariate values are not needed, when covariates are specified to be generated!\n")
        if(is.null(model_for_c)) stop("Models for covariates are needed to generate covariates!\n")
        if(is.null(parameters_for_c)) stop("Coefficients are needed to generate covariates!\n")
      }
    }else if(given.c=="null"){
      if(length(grep("^c$",model_for_y))==0 & length(grep("^c$",model_for_x))==0){
        GENERATE.C=F
      }else{
        if(length(grep("^c$",model_for_y))!=0) stop("Covariates are need for the generation of outcome. Please give a input or ask for generation!\n")
        if(length(grep("^c$",model_for_x))!=0) stop("Covariates are need for the generation of treatment. Please give a input or ask for generation!\n")
      }
    }
    
    # Check treatments: if treatment is given:
    #                   if treatment is randomized: only depend on baseline
    #                   if treatment depends on previous treatment, previous outcome, and current covariate
    if(!(type_for_x %in% c("continuous","binary","categorical"))){
      stop("The specified categorial of x can only choose continuous/binary/categorical!")
    }
    if(given.x==T){ # if x is given, just take it
      if(!is.null(input.x) & length(input.x) == n & is.null(parameters_for_x) & is.null(initial.x) & is.null(model_for_x)){
        GENERATE.X=F
        x=input.x
        if(printFlag){cat("The treatments are given! ")}
        if(type_for_x=="continuous"){
          if(printFlag){cat("The given x is continuous variable.\n")}
        }else if(type_for_x=="binary"){
          if(all(names(table(x)) %in% c("0","1"))){
            if(printFlag){cat("The given x is binary variable.\n")}
          }else{
            stop("Given binary x only allow values of 0/1.")
          }
        }else if(type_for_x=="categorical"){
          if(all(as.numeric(names(table(x)))%% 1 ==0)){
            if(printFlag){cat("The given x is categorical variable.")}
          }else{
            stop("Given categorical x are not all integers.")
          }
        }else{
          stop("The specified categorial of x can only choose continuous/binary/categorical!")
        }
      }else{
        if(is.null(input.x)) stop("Treatment information is needed, when we specify treatments as given!\n")
        if(!is.null(parameters_for_x)) stop("Coefficients to generate treatment are not needed, as treatments are given!\n")
        if(!is.null(initial.x)) stop("Initial treatment value is not needed, as treatments are given!\n")
        if(!is.null(model_for_x)) stop("Models for treatment generation is not needed, as treatments are given!\n")
        if(length(input.x) != n) stop("The length of input treatment is wrong!\n")
      }
    }else if(given.x==F){ # if x is not given, generate X
      # specifing "given.x==F" means that treatments need to be generated
      #   check (i) no input for X
      #         (ii) speficiation of stuff used for generating treatments: "parameters_for_x","model_for_x", "initial.x"(optional)
      if(is.null(input.x) & !is.null(parameters_for_x) & !is.null(model_for_x)){
        if(length(model_for_x) == (length(parameters_for_x)-1) & 
           "sd" %in% names(parameters_for_x) &
           "baseline" %in% names(parameters_for_x)
        ){
          if(printFlag){cat("The generation of treatments depend on: ", paste(c(unlist(lapply(grep("x|y",model_for_x,value=T),function(a){paste("previous",a)})),
                                                                                paste("current",grep("^c$",model_for_x,value=T))),
                                                                              collapse=", "),".\n",sep="")}
          GENERATE.X=T
          x=rep(NA,n)
          noise_x=rnorm(n,sd=parameters_for_x$sd)
          # cat("We generate",type_for_x,"variables for x!\n")
          if(type_for_x=="continuous"){
            if(printFlag){cat("We generate",type_for_x,"variables for x!\n")}
          }else if(type_for_x=="binary"){
            logit=p=rep(NA,n)
          }else if(type_for_x=="categorical"){
            # to be filled when it's necessary
          }else{
            stop("The specified categorial of x can only choose continuous/binary/categorical!")
          }
          
        }else if(length(model_for_x) != (length(parameters_for_x)-1)){
          stop("The dimention of coefficients for treatment generation doesn't agreen with its model!\n")
        }else{
          stop("Baseline or sd of noise is missing!\n")
        }
      }else{
        if(!is.null(input.x)) stop("Treatment values are not needed, when treatments are specified to be generated!\n")
        if(is.null(model_for_x)) stop("Models for treatment are needed to generate treatments!\n")
        if(is.null(parameters_for_x)) stop("Coefficients are needed to generate treatments!\n")
      }
    }
    
    if(given.y==T){ # outcome is given (almost never used, except for checking code errors)
      if(!is.null(input.y)  & length(input.y) == n & is.null(parameters_for_y) & is.null(initial.y) & is.null(model_for_y) & is.null(lag_x_on_y)){
        GENERATE.Y=F
        y=input.y
        if(printFlag){cat("The outcomes are given!\n")}
      }else{
        if(is.null(input.y)) stop("Outcome information is needed, when we specify outcomes as given!\n")
        if(!is.null(parameters_for_y)) stop("Coefficients to generate treatment are not needed, as outcomes are given!\n")
        if(!is.null(initial.y)) stop("Initial outcome value is not needed, as outcomes are given!\n")
        if(!is.null(model_for_y)) stop("Models for treatment generation is not needed, as outcomes are given!\n")
        if(!is.null(lag_x_on_y)) stop("#lag is not needed, as outcomes are given!\n")
        if(length(input.y) != n) stop("The length of input outcome is wrong!\n")
      }
    }else if(given.y==F){ # if outcome is not given, generate Y
      # specifing "given.y==F" means that outcomes need to be generated
      #   check (i) no input for Y
      #         (ii) speficiation of stuff used for generating treatments: "parameters_for_y","model_for_y", "lag_x_on_y","initial.y"(optional)
      if(is.null(input.y) & !is.null(parameters_for_y) & !is.null(model_for_y) & !is.null(lag_x_on_y) & lag_x_on_y>=0){
        if(length(model_for_y) == (length(parameters_for_y)-1) & 
           "sd" %in% names(parameters_for_y) &
           "baseline" %in% names(parameters_for_y)){
          
          if( is.matrix(parameters_for_y$x) & ncol(parameters_for_y$x)==lag_x_on_y & all.equal(colnames(parameters_for_y$x),paste("t",1:ncol(parameters_for_y$x)-1,sep="-")) ){
            if(printFlag){cat("The generation of outcomes depend on: ", paste(c(unlist(lapply(grep("y",model_for_y,value=T),function(a){paste("previous",a)})),
                                                                                paste("current",grep("x|c",model_for_y,value=T))),
                                                                              collapse=", "),".\n",sep="")}
            
            GENERATE.Y=T
            y=rep(NA,n)
            noise_y=rnorm(n,sd=parameters_for_y$sd)
            # logit=p=rep(NA,n)
          }else{
            if(!is.matrix(parameters_for_y$x)) stop("parameter_for_y$x should be a matrix, even if only one row.")
            if(ncol(parameters_for_y$x)!=lag_x_on_y) stop("The dimention of treatment coefficients doesn't agree with #lags!\n")
            if(!all.equal(colnames(parameters_for_y$x),paste("t",1:ncol(parameters_for_y$x)-1,sep="-"))) stop("The the column names of treatment effect matrix!\n")
          }
        }else if (length(model_for_y) != (length(parameters_for_y)-1)){
          stop("The dimention of coefficients for treatment generation doesn't agreen with its model!\n")
        }else{
          stop("Baseline or sd of noise is missing!\n")
        }
      }else{
        if(!is.null(input.y)) stop("Outcomes are not needed, when outcomes are specified to be generated!\n")
        if(is.null(model_for_y)) stop("Models for outcomes are needed to generate outcomess!\n")
        if(is.null(parameters_for_y)) stop("Coefficients are needed to generate outcomess!\n")
        if(is.null(lag_x_on_y)) stop("The period of treatment relavent is not specified!\n")
        if(lag_x_on_y<0) stop("The period of treatment relavent should >=0!\n")
      }
    }
    
    # a small function to help extract information from interaction_pairs
    interaction_info_table=function(interaction){
      if(length(interaction)==2 & all(str_detect(interaction,"^(x|y|c)$|^(x|y|c)\\-\\d$"))){
        a=as.data.frame(str_split(interaction,"\\-",simplify = T))
        if(ncol(a)==1){a=cbind(a,rep(NA,2)) }
        colnames(a)=c("variable","shift")
        a$variable=as.character(a$variable)
        a$shift=as.numeric(as.character(a$shift))
        a$shift[is.na(a$shift)]=0
        a$shift[grep("[\\-]",interaction)]=-a$shift[grep("\\-",interaction)]
        
        if("y" %in% a$variable){
          if(a$shift[a$variable=="y"]==0) stop("Current value of y hasn't been generated, and should be in the interaction term.")
        }
        return(a)
      }else{
        if(length(interaction)!=2) stop("The length of each interaction line should be 2.")
        if(!all(str_detect(interaction,"^(x|y|c)$|^(x|y|c)[\\-]\\d$"))) stop("Variables should follow x|y|c [-] number pattern.")
      }
    }
    
    
    if(!is.null(interaction_pairs)){
      if("interaction" %in% names(parameters_for_y) & length(parameters_for_y$interaction)==n){
        GENERATE.interaction=T
        interaction_info=interaction_info_table(interaction_pairs)
        if(!all(dim(interaction_info)==c(2,2))){
          stop("Something wrong with the generated interaction_info")
        }
        interaction=rep(NA,n)
      }else{
        if(!("interaction" %in% names(parameters_for_y))) stop("interaction coefficient is missing in parameter_for_y.")
        if(length(parameters_for_y$interaction)!=n) stop("the length of interaction coefficients differ from n.")
      }
    }else{
      GENERATE.interaction=F
    }
    
    # Main part
    for(i in 1:n){
      if(GENERATE.C){
        if(i==1 & !is.null(initial.c)){
          c[i]=initial.c
        }else if(i==1 & is.null(initial.c)){
          c[i]=parameters_for_c$baseline[i]+noise_c[i]
        }else if(i>1){
          names=names(parameters_for_c)[grep("x|y|c",names(parameters_for_c))]
          c[i]=parameters_for_c$baseline[i]+
            sum(unlist(lapply(names,function(name){parameters_for_c[[name]][i]*get(name)[i-1]})))+noise_c[i]
        }
      }
      
      if(GENERATE.X){
        if(i==1 & !is.null(initial.x)){
          x[i]=initial.x
        }else if(i==1 & is.null(initial.x)){
          names_c=names(parameters_for_x)[grep("^c$",names(parameters_for_x))]
          if(type_for_x=="continuous"){
            x[i]=parameters_for_x$baseline[i]+
              sum(unlist(lapply(names_c,function(name){parameters_for_x[[name]][i]*get(name)[i]})))+noise_x[i]
          }else if(type_for_x=="binary"){
            logit_temp=parameters_for_x$baseline[i]+
              sum(unlist(lapply(names_c,function(name){parameters_for_x[[name]][i]*get(name)[i]})))+noise_x[i]
            p_temp=exp(logit_temp)/(1+exp(logit_temp))
            u=runif(1)
            x[i]=(u<p_temp)*1
            logit[i]=logit_temp
            p[i]=p_temp
          }
        }else if(i>1){
          names_xy=names(parameters_for_x)[grep("x|y",names(parameters_for_x))]
          names_c=names(parameters_for_x)[grep("^c$",names(parameters_for_x))]
          
          if(type_for_x=="continuous"){
            x[i]=parameters_for_x$baseline[i]+
              sum(unlist(lapply(names_xy,function(name){parameters_for_x[[name]][i]*get(name)[i-1]})))+
              sum(unlist(lapply(names_c,function(name){parameters_for_x[[name]][i]*get(name)[i]})))+
              noise_x[i]
          }else if(type_for_x=="binary"){
            logit_temp=parameters_for_x$baseline[i]+
              sum(unlist(lapply(names_xy,function(name){parameters_for_x[[name]][i]*get(name)[i-1]})))+
              sum(unlist(lapply(names_c,function(name){parameters_for_x[[name]][i]*get(name)[i]})))+
              noise_x[i]
            p_temp=exp(logit_temp)/(1+exp(logit_temp))
            u=runif(1)
            x[i]=(u<p_temp)*1
            #cat("i=",i,"p=",p_temp,"u=",u,"x[i]=",x[i],"\n")
            logit[i]=logit_temp
            p[i]=p_temp
          }
        }
      }
      
      if(GENERATE.Y){
        if(i<lag_x_on_y & !is.null(initial.y)){
          y[i]=initial.y[i]
        }else if(i<lag_x_on_y & is.null(initial.y)){
          # generate before the cutoff
          if(i==1){
            names_x=names(parameters_for_y)[grep("x",names(parameters_for_y))]
            names_c=names(parameters_for_y)[grep("^c$",names(parameters_for_y))]
            if(GENERATE.interaction){
              # since i=1, for all interaction terms, they should be current values 
              if(all(interaction_info$shift==0)){
                interaction[i]=prod(unlist(lapply(1:2, function(l){get(interaction_info$variable[l])[interaction_info$shift[l]+i]})))
              }else{
                interaction[i]=NA
              }
            }
            y[i]=parameters_for_y$baseline[i]+
              sum(unlist(lapply(names_c,function(name){parameters_for_y[[name]][i]*get(name)[i]})))+
              sum(unlist(lapply(names_x,function(name){parameters_for_y[[name]][i,1]*get(name)[i]})))+
              ifelse(GENERATE.interaction,ifelse(!is.na(interaction[i]),interaction[i]*parameters_for_y$interaction[i],0),0)+
              noise_y[i]
          }else{
            names_y=names(parameters_for_y)[grep("y",names(parameters_for_y))]
            names_x=names(parameters_for_y)[grep("x",names(parameters_for_y))]
            names_c=names(parameters_for_y)[grep("^c$",names(parameters_for_y))]
            if(GENERATE.interaction){
              if(all(i+interaction_info$shift<=i) & all(i+interaction_info$shift>0)){
                interaction[i]=prod(unlist(lapply(1:2, function(l){get(interaction_info$variable[l])[interaction_info$shift[l]+i]})))
              }else{
                interaction[i]=NA
              }
            }
            y[i]=parameters_for_y$baseline[i]+
              sum(unlist(lapply(names_y,function(name){parameters_for_y[[name]][i]*get(name)[i-1]})))+
              sum(unlist(lapply(names_c,function(name){parameters_for_y[[name]][i]*get(name)[i]})))+
              sum(unlist(lapply(names_x,function(name){parameters_for_y[[name]][i,1:i]*rev(get(name)[1:i])})))+
              ifelse(GENERATE.interaction,ifelse(!is.na(interaction[i]),interaction[i]*parameters_for_y$interaction[i],0),0)+
              noise_y[i]
          }
        }else if(i>(lag_x_on_y-1)){
          # generate since the cutoff
          names_y=names(parameters_for_y)[grep("y",names(parameters_for_y))]
          names_x=names(parameters_for_y)[grep("x",names(parameters_for_y))]
          names_c=names(parameters_for_y)[grep("^c$",names(parameters_for_y))]
          if(GENERATE.interaction){
            if(all(i+interaction_info$shift<=i) & all(i+interaction_info$shift>0)){
              interaction[i]=prod(unlist(lapply(1:2, function(l){get(interaction_info$variable[l])[interaction_info$shift[l]+i]})))
            }else{
              interaction[i]=NA
            }
          }
          
          y[i]=parameters_for_y$baseline[i]+
            sum(unlist(lapply(names_y,function(name){parameters_for_y[[name]][i]*get(name)[i-1]})))+
            sum(unlist(lapply(names_c,function(name){parameters_for_y[[name]][i]*get(name)[i]})))+
            sum(unlist(lapply(names_x,function(name){parameters_for_y[[name]][i,1:lag_x_on_y]*rev(get(name)[ ((i-lag_x_on_y)+1):i])})))+
            ifelse(GENERATE.interaction,ifelse(!is.na(interaction[i]),interaction[i]*parameters_for_y$interaction[i],0),0)+
            noise_y[i]
          # cat("i=",i,parameters_for_y$baseline[i]," ",sum(unlist(lapply(names_y,function(name){parameters_for_y[[name]][i]*get(name)[i-1]})))," ",
          #     sum(unlist(lapply(names_c,function(name){parameters_for_y[[name]][i]*get(name)[i]})))," ",
          #     sum(unlist(lapply(names_x,function(name){parameters_for_y[[name]][i,1:lag_x_on_y]*rev(get(name)[((i-lag_x_on_y)+1):i])}))),"\n")
          
        }
      }
    }
    
    # save result
    if(type_for_x=="binary" & GENERATE.X==T){
      if(given.c=="null"){
        return(list(x=x,y=y,logit=logit,p=p))
      }else{
        return(list(c=c,x=x,y=y,logit=logit,p=p))
      }
    }else{
      if(given.c=="null"){
        return(list(x=x,y=y))
      }else{
        return(list(c=c,x=x,y=y))
      }
    }
    
  }
  # [sim_continuous_AR1] generate one single AR(1) time series
  sim_continuous_AR1=function(baseline,sd,rho,periods){
    # Generate a single auto-correlated continous variable: 
    #       can be used to generate gradually changed coefficients
    
    # Either baseline, sd and rho can be periodic-specific, independently
    #   -> single number or vector
    if(length(periods)==1 & is.numeric(periods)){
      # one period, and the #time specified by this period is >=2
      if( periods>=2 & length(baseline)==1 & length(sd)==1 & length(rho)==1 & rho>=0 & sd>0 ){
        cat("Only one period!")
        n=periods
        noise=rnorm(n,sd=sd)
        result=rep(NA,n)
        result[1]=baseline+noise[1]
        for(i in 2:n){
          result[i]=baseline*(1-rho)+rho*result[i-1]+noise[i]
        }
      }else{
        stop("The length of baseline, sd, rho should all be 1!\n")
      }
    }else if(length(periods)>1){
      # check the length of input parameter if they are 1 or same as the # periods
      if(length(baseline)!=1 & length(baseline) != length(periods)){
        stop("The length of baseline parameter is wrong!\n")
      }else if(length(sd)!=1 & length(sd)!=length(periods)){
        stop("The length of sd parameter is wrong!\n")
      }else if(length(rho)!=1 & length(rho)!=length(periods)){
        stop("The length of rho parameter is wrong!\n")
      }
      cat("There are",length(periods),"periods!\n")
      # for the first period
      noise_temp=rnorm(periods[1],sd=ifelse(length(sd)==1,sd,sd[1]))
      c_temp=rep(NA,periods[1])
      c_temp[1]=ifelse(length(baseline)==1,baseline,baseline[1])+noise_temp[1]
      for(i in 2:periods[1]){
        c_temp[i]=(1-ifelse(length(rho)==1,rho,rho[1]))*ifelse(length(baseline)==1,baseline,baseline[1])+
          ifelse(length(rho)==1,rho,rho[1])*c_temp[i-1]+noise_temp[i]
      }
      result=c_temp
      # for the second to the last periods
      for(j in 2:length(periods)){
        noise_temp=rnorm(periods[j],sd=ifelse(length(sd)==1,sd,sd[j]))
        c_temp=rep(NA,periods[j])
        c_temp[1]=(1-ifelse(length(rho)==1,rho,rho[j]))*ifelse(length(baseline)==1,baseline,baseline[j])+
          ifelse(length(rho)==1,rho,rho[j])*result[length(result)]+noise_temp[1]
        for(k in 2:periods[j]){
          c_temp[k]=(1-ifelse(length(rho)==1,rho,rho[j]))*ifelse(length(baseline)==1,baseline,baseline[j])+
            ifelse(length(rho)==1,rho,rho[j])*c_temp[k-1]+noise_temp[k]
        }
        result=c(result,c_temp)
      }
    }else{
      stop("The length of periods shoule be >=1!\n")
    }
    
    return(result)
  } # single auto-correlated continous variable
  # [sim_randomwalk] generase a random walk on top of a baseline
  sim_randomwalk=function(baseline,sd,periods){
    # Generate a single random-walk variable: 
    #       can be used to generate non-stationary changed coefficients
    
    #   -> single number or vector
    if(length(periods)==1 & is.numeric(periods)){
      # one period, and the #time specified by this period is >=2
      if( periods>=2 & length(baseline)==1 & length(sd)==1 & sd>0 ){
        cat("Only one period!")
        n=periods
        noise=cumsum(rnorm(n,sd=sd)) 
        result=baseline+noise
      }else{
        stop("The length of baseline, sd, rho should all be 1!\n")
      }
    }else if(length(periods)>1){
      # check the length of input parameter if they are 1 or same as the # periods
      #      baseline * n, sd = 1
      #      baseline * 1, sd = n
      #      baseline * n, sd = n
      if(length(baseline)!=1 & length(baseline) != length(periods)){
        stop("The length of baseline parameter is wrong!\n")
      }else if(length(sd)!=1 & length(sd)!=length(periods)){
        stop("The length of sd parameter is wrong!\n")
      }else if( length(baseline)==1 & length(sd)==1 ){
        stop("The length of both baseline and sd are 1! Period of 1 should be sufficient!\n")
      }
      cat("There are",length(periods),"periods!\n")
      # first, check sd 
      if(length(sd)>1){
        noise=c()
        for(i in 1:length(sd)){
          noise_temp=rnorm(periods[i],sd=sd[i])
          noise=c(noise,noise_temp)
        }
      }
      # second, add baseline
      if(length(baseline)>1){
        base=extend(baseline,periods)
      }
      result=base+cumsum(noise)
    }else{
      stop("The length of periods shoule be >=1!\n")
    }
    
    return(result)
    
  }
}

###################################################################
##           generate derived variables for data set             ##
###################################################################
check.one_obs_each_day=function(dates){
  if(mean(as.numeric(table(dates))==rep(1,length(as.numeric(table(dates)))))==1){
    return(T)
  }else{
    return(F)
  }
} # check a series of dates only happen for once
# Goal: add variables to long-format data, when there is at most one obs each day
# adding more variables
# 1. calculate past average over for possibly multiple variables
# 2. get lagged variables for possibly multiple variables
# 3. get dummy variables
# 4. create interactions between two variables
average.past_values=function(longdata,variables,param,na.rm=F,printFlag=T){
  # check if 1) <Date> is in the dataset
  #          2) at most one obs each day
  #          3) check if the variables requested are in the long data
  #          4) check the specification of param
  if(!("Date" %in% colnames(longdata))) stop("Variable [Date] is required for the data!\n")
  if(!check.one_obs_each_day(longdata$Date)){
    stop("More than one obs each day.")
  }
  # check if the variables requested are in the long data
  if(!all(variables %in% colnames(longdata))){
    stop("Selected variables are not in the data set")
  }
  # check the specification of param
  if(!(all(param>1) & length(param) == length(variables))){
    stop("''param'' specification is incorrect.")
  }
  
  if(printFlag){
    cat(red("Calculate average past values of variables!\n"))
    cat(blue("=======================================\n"))
    cat(blue("Data, varaible, and param check:\n"))
    alltime = colnames(longdata)[which(colnames(longdata) %in% c("EST.time","UTC.time","Date","Time"))]
    cat(blue(" -> 1) The dataset contains the following ''time'' relevant varibles: ",paste(alltime,collapse = ", "),". \n",sep=""))
    cat(blue(" -> 2) At most one obs each day.\n"))
    cat(blue(" -> 3) Selected variables are inside data set.\n"))
    cat(blue(" -> 4) ''param'' specification is correct.\n"))
  }
  
  longdata=longdata[order(longdata$Date),]
  
  # # calculate the average into the past
  # past_average=function(k,value){
  #   result=rep(NA,length(value))
  #   for(i in 1:k){
  #     result[i]=mean(value[1:i])
  #   }
  #   for(i in (k+1):length(value)){
  #     result[i]=mean(value[(i-k+1):i])
  #   }
  #   return(result)
  # }
  # calculate the average into the past
  past_average=function(k,value,na.rm){
    result=rep(NA,length(value))
    for(i in 1:k){
      result[i]=mean(value[1:i],na.rm=na.rm)
    }
    for(i in (k+1):length(value)){
      result[i]=mean(value[(i-k+1):i],na.rm=na.rm)
    }
    return(result)
  }
  temp=as.data.frame(matrix(NA,nrow=nrow(longdata),ncol=length(variables)))
  for(i in 1:length(variables)){
    temp[,i]=past_average(param[i],longdata[,variables[i]],na.rm=na.rm)
  }
  result=cbind(longdata$Date,temp)
  colnames(result)=c("Date",paste(variables,"_ave",sep=""))
  if(printFlag){
    cat("Calculation for average past values complete.\n")
    cat(blue("=======================================\n"))
  }
  # return the result 
  return(result)
}
calculate.lagged_values=function(longdata,variables,param,printFlag=T){
  # check if 1) <Date> is in the dataset
  #          2) at most one obs each day
  #          3) check if the variables requested are in the long data
  #          4) check the specification of param
  if(!("Date" %in% colnames(longdata))) stop("Variable [Date] is required for the data!\n")
  
  if(!check.one_obs_each_day(longdata$Date)){
    stop("More than one obs each day.")
  }
  # check if the variables requested are in the long data
  if(!all(variables %in% colnames(longdata))){
    stop("Selected variables are not in the data set.")
  }
  # check the specification of param
  if(!(all(param>0) & length(param) == length(variables))){
    stop("''param'' specification is incorrect.")
  }
  if(printFlag){
    cat(red("Calculate lagged values of variables!\n"))
    cat(blue("=======================================\n"))
    cat(blue("Data, varaible, and param check:\n"))
    alltime = colnames(longdata)[which(colnames(longdata) %in% c("EST.time","UTC.time","Date","Time"))]
    
    cat(blue(" -> 1) The dataset contains the following ''time'' relevant varibles: ",paste(alltime,collapse = ", "),". \n",sep=""))
    cat(blue(" -> 2) At most one obs each day.\n"))
    cat(blue(" -> 3) Selected variables are inside data set.\n"))
    cat(blue(" -> 4) ''param'' specification is correct.\n"))
  }
  
  longdata=longdata[order(longdata$Date),]
  
  # calculate the lagged variables
  shiftback=function(a,shift,name){
    m=length(a)
    result=matrix(NA,nrow=m,ncol=shift)
    result[,1]=c(NA,a[1:(m-1)])
    if(shift>1){
      for(i in 2:shift){
        result[,i]=c(NA,result[1:(m-1),i-1])
      }
    }
    colnames(result)=paste(name,1:shift,sep="_")
    return(as.data.frame(result))
  }
  result=shiftback(longdata[,variables[1]],shift=param[1],name=variables[1])
  if(length(variables)>1){
    for(i in 2:length(variables)){
      temp=shiftback(longdata[,variables[i]],shift=param[i],name=variables[i])
      result=cbind(result,temp)
    }
  }
  result=cbind(longdata$Date,result)
  colnames(result)[1]="Date"
  if(printFlag){
    cat("Calculation of lagged variables complete!\n")
    cat(blue("=======================================\n"))
  }
  return(result)
}
create.dummy_values=function(longdata,variables,printFlag=T){
  # check if 1) <Date> is in the dataset
  #          2) at most one obs each day
  #          3) check if the variables requested are in the long data
  #          4) check the specification of param
  if(!("Date" %in% colnames(longdata))) stop("Variable [Date] is required for the data!\n")
  if(!check.one_obs_each_day(longdata$Date)){
    stop("More than one obs each day.")
  }
  # check if the variables requested are in the long data
  if(!all(variables %in% colnames(longdata))){
    stop("Selected variables are not in the data set.")
  }
  
  if(printFlag){
    cat(red("Calculate dummy variables!\n"))
    cat(blue("=======================================\n"))
    cat(blue("Data, varaible, and param check:\n"))
    alltime = colnames(longdata)[which(colnames(longdata) %in% c("EST.time","UTC.time","Date","Time"))]
    cat(blue(" -> 1) The dataset contains the following ''time'' relevant varibles: ",paste(alltime,collapse = ", "),". \n",sep=""))
    cat(blue(" -> 2) At most one obs each day.\n"))
    cat(blue(" -> 3) Selected variables are inside data set.\n"))
  }
  
  longdata=longdata[order(longdata$Date),]
  
  # check the specification of param
  for(k in 1:length(variables)){
    if(all(as.numeric(names(table(longdata[,variables[k]]))) %% 1 ==0)){
      if(printFlag){
        cat(blue("The levels for variable ",variables[k]," is: ",paste(names(table(longdata[,variables[k]])),collapse = ", "),"\n",sep=""))
      }
    }else{
      stop("The levels of variable",variables[k],"is not a integer!")
    }
  }
  
  create.dummy=function(value,name){
    temp=dummy_cols(value)
    temp=temp[,-1]
    colnames(temp)=gsub(".data_",paste(name,".",sep=""),colnames(temp))
    return(temp)
  }
  
  result=create.dummy(value=longdata[,variables[1]],name=variables[1])
  if(length(variables)>1){
    for(i in 2:length(variables)){
      temp=create.dummy(longdata[,variables[i]],name=variables[i])
      result=cbind(result,temp)
    }
  }
  result=cbind(longdata$Date,result)
  colnames(result)[1]="Date"
  result=result[,-grep(".NA",colnames(result))]
  
  if(printFlag){
    cat("Calculation of dummy variables complete!\n")
    cat(blue("=======================================\n"))
  }
  return(result)
}
create.interaction_dummy_values=function(longdata,variable1,variable2,var_type=NULL,printFlag=T){
  
  # check if 1) <Date> is in the dataset
  #          2) at most one obs each day
  #          3) check if the variables requested are in the long data
  #          4) check the specification of param
  if(!("Date" %in% colnames(longdata))) stop("Variable [Date] is required for the data!\n")
  if(!check.one_obs_each_day(longdata$Date)){
    stop("More than one obs each day.")
  }
  # check if the variables requested are in the long data
  if(!all(c(variable1,variable2) %in% colnames(longdata))){
    stop("Selected variables are not in the data set.\n")
  }
  # check the length of two data are the same!
  if(!all(c(length(longdata[,variable1]),length(longdata[,variable2]))==nrow(longdata))){
    stop("The length of the two variables don't match!")
  }
  longdata=longdata[order(longdata$Date),]
  
  if(printFlag){
    cat(red("Calculate interaction dummy variables for [",variable1,"] and [", variable2,"].\n",sep=""))
    cat(blue("=======================================\n"))
    cat(blue("Data, varaible, and param check:\n"))
    alltime = colnames(longdata)[which(colnames(longdata) %in% c("EST.time","UTC.time","Date","Time"))]
    cat(blue(" -> 1) The dataset contains the following ''time'' relevant varibles: ",paste(alltime,collapse = ", "),". \n",sep=""))
    cat(blue(" -> 2) At most one obs each day.\n"))
    cat(blue(" -> 3) Selected variables are inside data set.\n"))
  }
  
  # check the type of two variables
  if(is.null(var_type)){
    var_type=rep(NA,2)
    if(all(as.numeric(names(table(longdata[,variable1]))) %% 1 ==0)){
      var_type[1]="categorical"
      if(printFlag){
        cat(blue("The levels for variable ",variable1," is: ",paste(names(table(longdata[,variable1])),collapse = ", "),"\n",sep=""))
      }
    }else{
      var_type[1]="continuous"
      var_range=c(min(range(longdata[,variable1])),max(range(longdata[,variable1])))
      if(printFlag){
        cat(blue("The range for variable ",variable1," is: ",paste(var_range,collapse = ", "),"\n",sep=""))
      }
    }
    if(all(as.numeric(names(table(longdata[,variable2]))) %% 1 ==0)){
      var_type[2]="categorical"
      if(printFlag){
        cat(blue("The levels for variable ",variable2," is: ",paste(names(table(longdata[,variable2])),collapse = ", "),"\n",sep=""))
      }
    }else{
      var_type[2]="continuous"
      var_range=c(min(range(longdata[,variable2],na.rm =T)),max(range(longdata[,variable2],na.rm =T)))
      if(printFlag){
        cat(blue("The range for variable ",variable2," is: ",paste(var_range,collapse = ", "),"\n",sep=""))
      }
    }
  }
  if(printFlag){
    cat("The variable types are: ",paste(var_type,collapse = ", "),".\n",sep="")
  }
  if(length(var_type)!=2) stop("The length of ''var_type'' is wrong!\n")
  
  if(all(var_type=="categorical")){
    variable1.level=as.character(names(table(longdata[,variable1])))
    variable2.level=as.character(names(table(longdata[,variable2])))
    levels=list(variable1=variable1.level,variable2=variable2.level)
    result=longdata$Date
    for(i in 2:length(levels[[1]])) for(j in 2:length(levels[[2]])){
      if(printFlag){
        cat("Generate interactions for [",variable1,"]=",levels[[1]][i],", [",variable2,"]=",levels[[2]][j],"\n",sep="")
      }
      temp=rep(0,length(longdata[,variable1]))
      temp[longdata[,variable1]==levels[[1]][i] & longdata[,variable2]==levels[[2]][j]]=1
      result=data.frame(result,temp)
      colnames(result)[ncol(result)]=paste(paste(variable1,".",levels[[1]][i],sep=""),"_and_",paste(variable2,".",levels[[2]][j],sep=""),sep="")
    }
  }else if(all(var_type=="continuous")){
    result=longdata$Date
    temp=longdata[,variable1]*longdata[,variable2]
    result=data.frame(result,temp)
    colnames(result)[ncol(result)]=paste(variable1,"_and_",variable2,sep="")
  }else if(all(var_type %in% c("categorical","continuous"))){
    result=longdata$Date
    index=which(var_type=="categorical")
    levels=as.character(names(table(longdata[,get(paste("variable",2,sep=""))])))
    if(index==1){ # the first variable is categorical
      for(m in 2:length(levels)){
        if(printFlag){
          cat("Generate interactions for [",variable2,"] and [",variable1,"]=",levels[m],"\n",sep="")
        }
        temp=as.numeric(longdata[,variable1]==levels[m])*longdata[,variable2]
        result=data.frame(result,temp)
        colnames(result)[ncol(result)]=paste(variable1,".",levels[m],"_and_",variable2,sep="")
      }
    }else if(index==2){# the second variable is categorical
      for(m in 2:length(levels)){
        if(printFlag){
          cat("Generate interactions for [",variable1,"] and [",variable2,"]=",levels[m],"\n",sep="")
        }
        temp=as.numeric(longdata[,variable2]==levels[m])*longdata[,variable1]
        result=data.frame(result,temp)
        #colnames(result)[ncol(result)]=paste(variable2,".",levels[m],"_and_",variable1,sep="")
        colnames(result)[ncol(result)]=paste(variable1,"_and_",variable2,".",levels[m],sep="")
      }
    }
  }else{
    stop("''var_type'' should be either continuous or categorical!\n")
  }
  colnames(result)[1]="Date"
  return(result)
}
calculate.consecutive_values=function(longdata,variables,seq,selected_values,printFlag=T){
  # Example: variables=c("Social_in_Person","Social_Digitally")
  #          seq=list(Social_in_Person=c(2,3),Social_Digitally=c(2,3))
  #          selected_values=list(Social_in_Person=c(0,0),Social_Digitally=c(0,0))
  
  # check if 1) <Date> is in the dataset
  #          2) at most one obs each day
  #          3) check if the variables requested are in the long data
  #          4) check the specification of param
  if(!("Date" %in% colnames(longdata))) stop("Variable [Date] is required for the data!\n")
  if(!check.one_obs_each_day(longdata$Date)){
    stop("More than one obs each day.")
  }
  # check if the variables requested are in the long data
  if(!all(variables %in% colnames(longdata))){
    stop("Selected variables are not in the data set")
  }
  # check the specification of param
  if(length(selected_values)!=length(variables) | length(seq)!=length(variables)){
    stop("The length of [selected_values] or [seq] differs from that of [variables].")
  }
  for(i in 1:length(selected_values)){
    if(!all(c(names(selected_values)[i],names(seq)[i])==variables[i])) stop("Variable names don't match.")
    if(!all(seq[[i]]>1)) stop("The selected number of consecution should be greater than 1.")
  }
  
  longdata=longdata[order(longdata$Date),]
  
  if(printFlag){
    cat(red("Calculate consecutive variables!\n"))
    cat(blue("=======================================\n"))
    cat(blue("Data, varaible, and param check:\n"))
    alltime = colnames(longdata)[which(colnames(longdata) %in% c("EST.time","UTC.time","Date","Time"))]
    cat(blue(" -> 1) The dataset contains the following ''time'' relevant varibles: ",paste(alltime,collapse = ", "),". \n",sep=""))
    cat(blue(" -> 2) At most one obs each day.\n"))
    cat(blue(" -> 3) Selected variables are inside data set.\n"))
  }
  
  consecutive_value=function(one_data,k,value){
    result=rep(NA,length(one_data))
    for(i in k:length(one_data)){
      result[i]=all(one_data[(i-k+1):i]==value)*1
    }
    return(result)
  }
  consecutive_values=function(one_data,seq,value,name){
    result=matrix(NA,ncol=length(seq),nrow=length(one_data))
    for(j in 1:length(seq)){
      result[,j]=consecutive_value(one_data,seq[j],value[j])
    }
    result=as.data.frame(result)
    colnames(result)=paste(name,".",value,"_multi_",seq,sep="")
    
    return(result)
  }
  
  result=longdata[,"Date"]
  for(k in 1:length(selected_values)){
    temp=consecutive_values(longdata[,variables[k]],seq=seq[[k]],value=selected_values[[k]],variables[k])
    result=cbind(result,temp)
  }
  colnames(result)[1]="Date"
  if(printFlag){
    cat("Calculation of consecutive values complete!\n")
    cat(blue("=======================================\n"))
  }
  return(result)
}
add.varaible_into_longformat_data=function(longdata,average_param=NA,dummy_param=NA,lagged_param=NA,interaction_param=NA,consecutive_param=NA,printFlag=T){
  
  # check if 1) <Date> is in the dataset
  #          2) at most one obs each day
  if(class(longdata)!="data.frame"){stop("the data should be of a data frame.")}
  if(!("Date" %in% colnames(longdata))) stop("Variable [Date] is required for the data!\n")
  if(!check.one_obs_each_day(longdata$Date)){
    stop("More than one obs each day.")
  }
  longdata=longdata[order(longdata$Date),]
  
  if(printFlag){
    # check if one obs at most each day
    cat(blue("=======================================\n"))
    cat(blue("Data, varaible, and param check:\n"))
    alltime = colnames(longdata)[which(colnames(longdata) %in% c("EST.time","UTC.time","Date","Time"))]
    cat(blue(" -> 1) The dataset contains the following ''time'' relevant varibles: ",paste(alltime,collapse = ", "),". \n",sep=""))
    cat(blue(" -> 2) At most one obs each day.\n"))
    cat(blue("=======================================\n"))
  }
  
  # for each category and the requested variable, add in varaibles
  result=longdata
  if(!all(is.na(average_param))){
    if(is.null(average_param$na.rm)){
      cat("na.rm=False")
      average_param[["na.rm"]]=F
    }
    part_average=average.past_values(longdata,variables=average_param$variables,param=average_param$param,na.rm=average_param$na.rm,printFlag=printFlag)
    result=merge(result,part_average,by="Date")
  }
  if(!all(is.na(dummy_param))){
    part_dummy=create.dummy_values(longdata,variables=dummy_param,printFlag=printFlag)
    result=merge(result,part_dummy,by="Date")
  }
  if(!all(is.na(lagged_param))){
    part_lagged=calculate.lagged_values(longdata,variables=lagged_param$variables,param=lagged_param$param,printFlag=printFlag)
    result=merge(result,part_lagged,by="Date")
  }
  if(!all(is.na(interaction_param))){
    part_interaction=create.interaction_dummy_values(longdata,variable1=interaction_param$variable1,
                                                     variable2=interaction_param$variable2,
                                                     var_type=interaction_param$var_type,printFlag=printFlag)
    result=merge(result,part_interaction,by="Date")
  }
  if(!all(is.na(consecutive_param))){
    part_consecutive=calculate.consecutive_values(longdata,consecutive_param$variables,
                                                  seq=consecutive_param$seq,
                                                  selected_values=consecutive_param$selected_values,printFlag=printFlag)
    result=merge(result,part_consecutive,by="Date")
  }
  return(result)
}
# given series of command in  "add.varaible_into_longformat_data" (calculate past average,get lagged variables,get dummy variables,create interactions)
#   -> output the final long format data
# Examples for "param"
# param=list(list(average_param=list(variables=c("Physically_Active"),param=7),
#                 lagged_param=list(variables=c("Social_in_Person","Social_Digitally","Sleeping","outgoing_calls","outgoing_texts","interview","negative_total"),param=c(rep(1,5),2,2))),
#            list(dummy_param=c("Social_in_Person","Social_Digitally","Social_in_Person_1","Social_Digitally_1","Sleeping"),
#                 lagged_param=list(variables=c("Physically_Active_ave"),param=1),
#                 interaction_param=list(variable1="Social_in_Person",variable2="Social_Digitally",var_type=NULL))
#            )
add_variables_procedures=function(longdata,param,printFlag=T){
  collection=c("average_param","dummy_param","lagged_param","interaction_param","consecutive_param")
  result=longdata
  for(i in 1:length(param)){
    if(!all(names(param[[i]]) %in% collection)){
      stop("Error in parameter specification!")
    }
    param_temp=param[[i]]
    missed=collection[!(collection %in% names(param[[i]]))]
    for(j in 1:length(missed)){
      param_temp[length(param_temp)+j]=NA
    }
    names(param_temp)[(length(param[[i]])+1):(length(param[[i]])+length(missed))]=missed
    result=add.varaible_into_longformat_data(result,average_param=param_temp$average_param,
                                             lagged_param=param_temp$lagged_param,
                                             dummy_param=param_temp$dummy_param,
                                             interaction_param=param_temp$interaction_param,
                                             consecutive_param=param_temp$consecutive_param,printFlag=printFlag)
    
  }
  return(result)
}


#########################################################
##            small functions for unknown use          ##
#########################################################
# Impute missing values of all columns with NA, whose imputation method include:
#     "mean","median","mode","locf","nocb"
#     "simple":simple moving average
#     "exponential":exponential weighted moving average
#     "linear":linear interpolation
#     "spline":spline interpolation
#     "stine": stine interpolation
#     "StructTS": structural model & kalman sommthing
#     "auto.arima": ARIMA state space representation & Kalman smoothing
impute.all_columns=function(longdata,method){
  # longdata: for each column, if there is a missing value, impute using the according method
  #                            if there is no missing values, skip
  # Caution: 1) the longdata may contains many variables relevant to time, eg. "EST.time","Time","UTC.time","end","start","est.time","time","utc.time",
  #             which do not make sense to impute
  #          2) the longdata may contain missing days with no rows at all, or several rows for the sameday, which may be problematic for analysis
  # The purpose of this function is just do imputation, regardless of possible problems in Caution.
  
  vars_time=grep("time|Date|start|end",colnames(longdata),value=T)
  if(length(vars_time)>0){
    cat("The data set has the following time-relavant columns: ",paste(vars_time,collapse = ", "),".\n",sep="")
    if(length(vars_time)==1){
      columns_time=data.frame(longdata[,vars_time])
      colnames(columns_time)=vars_time
    }else{
      columns_time=longdata[,vars_time]
    }
    for(var_time in vars_time){
      temp=columns_time[,var_time]
      if(sum(diff(temp[!is.na(temp)])<0)!=0){cat("    The time column",var_time,"is not properly ordered.\n")}
    }
    if("Date" %in% vars_time){
      if(check.one_obs_each_day(longdata$Date)){cat("    [Date] column has at most one obs each day\n")}else{cat("    [Date] column has more than one obs for some days.\n")}
      if((diff(range(unique(longdata$Date)))+1)==nrow(longdata)){cat("    [Date] do not miss any day in the range of the observation.\n")}else{cat("    [Date] do not have certain days in the range of observation.\n")}
    }
    longdata_notimepart=longdata[,-grep("time|Date|start|end",colnames(longdata))]
  }else{
    cat("The data set doesn't have any time-relevant columns.\n")
    longdata_notimepart=longdata
  }

  method_collection=c("mean","median","mode","locf","nocb","simple","exponential","linear","spline","stine","StructTS","auto.arima")
  if(!all(method %in% method_collection)){
    stop("The selected imputations method are not supported.\nPlease select from",paste(method_collection,collapse = ", "))
  }
  if(length(method)==1){
    method=rep(method,ncol(longdata_notimepart))
  }else if(length(method)!=ncol(longdata_notimepart)){
    stop("The length of ''method'' don't equal to the length of time-irrelavant columns.\n")
  }
  
  impute.single_column=function(temp_column,method){
    if(sum(is.na(temp_column))==0){
      return(temp_column)
    }else{
      if(method %in% c("mean","median","mode")){
        temp_column_new=na.mean(temp_column, option = method)
      }else if(method %in% c("locf","nocb")){
        temp_column_new=na.locf(temp_column, option = method)
      }else if(method %in% c("simple","exponential")){
        temp_column_new=na.ma(temp_column, option = method)
      }else if(method %in% c("linear","spline","stine")){
        temp_column_new=na.interpolation(temp_column, option = method)
      }else if(method %in% c("StructTS","auto.arima")){
        temp_column_new=na.kalman(temp_column, model = method)
      }
      return(temp_column_new)
    }
  }
  result=longdata_notimepart
  for(k in 1:ncol(longdata_notimepart)){
    cat("Imput",colnames(result)[k],"using",method[k],"imputation.\n")
    result[,k]=impute.single_column(longdata_notimepart[,k],method=method[k])
  }
  if(length(vars_time)>0){
    result=data.frame(columns_time,result)
  }
  return(result)
  
}
{
  # take in long data -> make every colomn as time series format -> impute the data
  # imputation method: "Last observation carried forward", 
  #                    "Next observation carried backward",
  #                    "Linear interpolation"
  #                    "Spline interpolation"
  #                    "Median"
  # each column needs specification of imputation
  # impute.long_into_ts=function(longdata,method){
  #   # Check the variables: 1. only one "time" relevant varible [Date] is needed; all others should be removed
  #   #                      2. [Date] should be put at the first column
  #   cat(blue("Caution: only ''Date'' is allowed and should be put at the first column!\n"))
  #   time_not_necessary=which(colnames(longdata) %in% c("EST.time","Time","UTC.time","end","start","est.time","time","utc.time"))
  #   if(length(time_not_necessary)>0){
  #     cat(blue("Those unnecessary ''time'' relevant variables [", paste(colnames(longdata)[time_not_necessary],collapse = ","), "] are also included!\n"),sep="")
  #     longdata=longdata[,-time_not_necessary]
  #   }
  #   time_not_necessary=which(colnames(longdata) %in% c("EST.time","Time","UTC.time","end","start","est.time","time","utc.time"))
  #   if(length(time_not_necessary)==0){
  #     cat(blue("Unnecessary ''time'' relevant variables are removed.\n"))
  #   }
  #   if(!("Date" %in% colnames(longdata))){
  #     stop("Variable ''Date'' is needed for the dataset!\n")
  #   }else if(colnames(longdata)[1]!="Date"){
  #     stop("Variable ''Date'' should be at the first column!\n")
  #   }else{
  #     cat(blue("Variable ''Date'' exists at the first column!\n"))
  #   }
  #   
  #   # method should be specified for each variable, including outcome, exposure and confounders
  #   if((ncol(longdata)-1) == length(method) & all(method %in% c("median","locf","nocb","linear","spline","none"))){
  #     cat("The length of ''method'' equals the number of variables!\n")
  #     cat("The format of ''method'' is correct!\n")
  #   }else{
  #     stop("The length or format of ''method'' is wrong!\n")
  #   }
  #   # Next check imputation method
  #   result=longdata[,1]
  #   for(k in 2:length(colnames(longdata))){
  #     temp=ts(longdata[,k],frequency=7)
  #     if(method[k-1]=="median"){
  #       temp_new=na.mean(temp, option = "median")
  #     }else if(method[k-1]=="locf"){
  #       temp_new=na.locf(temp, option = "locf")
  #     }else if(method[k-1]=="nocb"){
  #       temp_new=na.locf(temp, option = "nocb") 
  #     }else if(method[k-1]=="linear"){
  #       temp_new=na.interpolation(temp, option = "linear")
  #     }else if(method[k-1]=="spline"){
  #       temp_new=na.interpolation(temp, option = "spline")
  #     }else if(method[k-1]=="none"){
  #       temp_new=temp
  #     }
  #     result=cbind(result,temp_new)
  #   }
  #   colnames(result)=colnames(longdata)
  #   return(result)
  # }  
} # replaced by [impute.all_columns]

# move the column with the name of [column] to be first column in the dataset
move_to_first=function(data,column){
  if(class(data)!="data.frame"){"The data, whose column to be moved, should be a data frame."}
  data=data[,c(which(colnames(data)==column),which(colnames(data)!=column))]
  return(data)
}
# move the variable [column] to the second
move_to_second=function(data,column){
  if(class(data)!="data.frame"){"The data, whose column to be moved, should be a data frame."}
  first=colnames(data)[1]
  data=data[,c(which(colnames(data)==first),which(colnames(data)==column),which(!(colnames(data) %in% c(first,column))))]
  return(data)
}
# find the # consecutive NA, to help learn the missing pattern
count_consecutive_NA=function(a){
  test=is.na(a)
  end_locations=cumsum(rle(test)$lengths)[rle(test)$values]
  count=rle(test)$lengths[rle(test)$values]
  begin_locations=end_locations-count+1
  return(data.frame(begin_locations=begin_locations,
                    end_locations=end_locations,
                    count= count))
}
{
# # this function is originally created for traditional time series analys.
# #     but it is no longer necessary, as we directly load the well-cleaned data, with all missing days filled and errors corrected
# # take out variables of interest [questions] into wide format, fill in missing days as NA, and also calculate the total score
# wideformat=function(df,questions,var_to_total=NULL){
#   # check if [var_to_total] <= [questions]
#   if(is.list(var_to_total)==T){
#     cat(green("[var_to_total] is a list -> multiple total scores will be caulated!\n"))
#     addrow=length(var_to_total)
#     # check the values, whose sum is to be calculated, are in the list of questions
#     if(all(unlist(lapply(1:2,function(i) all(var_to_total[[i]] %in% questions))))){
#       cat(green("Check! [var_to_total] <= [questions]\n"))
#     }else{
#       stop("[var_to_total] includes variables not in the [questions]")
#     }
#   }else if(is.vector(var_to_total)==T){
#     cat(green("[var_to_total] is a vector -> only one total score will be caulated!\n"))
#     addrow=1
#     if(mean(var_to_total %in% questions)==1){
#       cat(green("Check! [var_to_total] <= [questions]\n"))
#     }else{
#       stop("[var_to_total] includes variables not in the [questions]")
#     }
#   }
#   # for time series data, check if at most one obs each day
#   cat("=======================================\n")
#   result = df[df$question.text %in% questions,]
#   if(check.one_obs_each_day(unique(result[,c("Date","Time")])[,1])){
#     cat("At most one observation each day!\n")
#   }else{
#     stop("We require at most one obs each day!")
#   }
#   if(all(order(result$EST.time)==order(result$Date,result$Time))){
#     cat("The order by EST.time agrees with the order by Date+Time -> at most one obs each day!\n")
#   }else{
#     stop("We require at most one obs each day!")
#   }
#   result = result[with(result, order(EST.time,Date,Time,question.text)),]
#   target_total=day(days(max(range(result$Date))-min(range(result$Date))))+1
#   cat("The time range is: ",format(min(range(result$Date)),"%Y/%m/%d"),"->",format(max(range(result$Date)),"%Y/%m/%d"),
#       " (",target_total," days).\n",sep="")
#   cat("However, Only ",length(unique(result$Date))," days have survey answers.\n",sep="")
#   cat("=======================================\n")
#   cat(blue("Expand missing days horizontally begin!\n"))
#   
#   # impute missing dates and expand data horizontally
#   result = result[,c("Date","question.text","answer_number")] %>%
#     complete(Date = seq.Date(min(Date), max(Date), by="day"), question.text) %>%
#     complete(Date = seq.Date(min(Date), max(Date), by="day"), question.text) %>%
#     spread(Date,answer_number) %>% as.data.frame
#   rownames(result)=result[,1]
#   result=result[,-1]  
#   
#   # calcualte some total scores if required
#   if(!is.null(var_to_total)){
#     cat(blue("Calculate total score begin!"))
#     if(is.list(var_to_total)==T){
#       cat(blue(" -> Calculate total score for",paste(names(var_to_total),collapse =" and "),"mood.\n"))
#       total=matrix(NA,nrow=length(var_to_total),ncol=ncol(result))
#       for(k in 1:length(var_to_total)){
#         total[k,]=colSums(result[var_to_total[[k]],])
#       }
#       total=as.data.frame(total)
#       rownames(total)=paste(names(var_to_total),"_total",sep="")
#       colnames(total)=colnames(result)
#       result=rbind(result,total)
#     }else if(is.vector(var_to_total)==T){
#       cat(blue(" -> Calculate total score for only one variable.\n"))
#       total=colSums(result[var_to_total,])
#       result=rbind(result,total)
#       rownames(result)[nrow(result)]="total"
#     }
#   }
#   
#   # check result
#   if(mean(dim(result)==c(length(questions)+addrow,target_total))==1){
#     cat(green("The dimention of the result is correct.\n"))
#   }else{
#     stop("The dimention of the result is wrong!")
#   }
#   if(length(unique(colnames(result))) == target_total ){
#     cat(green("The final length of observation time is correct.\n"))
#   }else{
#     stop("The final length of observation is wrong!")
#   }
#   cat("After expansion, we have ",length(unique(colnames(result)))," days of survey answers.\n",sep="")
#   
#   return(result)
# }
} # no longer needed

# plot the different ways of imputation for a single ts_obj
plot.single_ts_under_different_imputation=function(ts_obj,methods,by_region=T){
  if(!(class(ts_obj) %in% c("ts","numeric"))){stop("the object to be imputed should be a ts or a vector.")}
  method_collection=c("mean","median","mode","locf","nocb","simple","exponential","linear","spline","stine","StructTS","auto.arima")
  if(!all(methods %in% method_collection)){
    stop("The selected imputations method are not supported.\nPlease select from",paste(method_collection,collapse = ", "))
  }
  impute.single_column=function(temp_column,method){
    if(sum(is.na(temp_column))==0){
      return(temp_column)
    }else{
      if(method %in% c("mean","median","mode")){
        temp_column_new=na.mean(temp_column, option = method)
      }else if(method %in% c("locf","nocb")){
        temp_column_new=na.locf(temp_column, option = method)
      }else if(method %in% c("simple","exponential")){
        temp_column_new=na.ma(temp_column, option = method)
      }else if(method %in% c("linear","spline","stine")){
        temp_column_new=na.interpolation(temp_column, option = method)
      }else if(method %in% c("StructTS","auto.arima")){
        temp_column_new=na.kalman(temp_column, model = method)
      }
      return(temp_column_new)
    }
  }
  
  result=list()
  par(mfrow=c(ceiling((length(methods)+1)/2),2))
  for(one_method in methods){
    result[[one_method]]=impute.single_column(ts_obj,one_method)
    invisible(plot(result[[one_method]],main=paste(one_method,"imputation")))
  }
  
  if(by_region){
    par(mfrow=c(2,1))
    col=c("red","blue","green","brown")
    for(i in 0:floor(length(ts_obj)/100)){
      invisible(plot(ts_obj[seq((i*100+1),(i+1)*100,1)],cex=1,type="l",ylab=paste(i*100+1,"to",(i+1)*100)))
      for(k in 1:length(result)){
        invisible(points(result[[k]][seq((i*100+1),(i+1)*100,1)],col=k+1,cex=(length(result)-k+1)/4))
      }
      legend("topright",legend=c("original",methods),col=1:(length(result)+1),lty=rep(1,length(result)),bty="n")
    }
  }
  return(result)
}
# help to grab variable names from SSModel (KFAS package)
#      including proper reformatting of factered variables
#      for example "factor(Social_in_Person)1"->"Social_in_Person.1"
grab_X=function(model){
  if(class(model)!="SSModel"){stop("Model should be a SSModel.")}
  gsub("\\)",".",gsub("factor\\(*","",colnames(model$T)))[-1]
}

# show some more readable about KFAS SSmodel
summary.SSModel=function(SSModel_subject){
  if(class(SSModel_subject)!="SSModel"){
    stop("The input subject is not a SSModel subject!")
  }
  cat("==========================================================\n")
  cat("The time-varying F matrix for the observation equation is", dim(SSModel_subject$Z)[1], "*", dim(SSModel_subject$Z)[2], "for",dim(SSModel_subject$Z)[3],"timepoints!\n") 
  cat("One example for the F matrix shows as:\n")
  print(SSModel_subject$Z[,,1])
  cat("\n")
  cat("==========================================================\n")
  cat("The variance of the observation equation, V, is:",SSModel_subject$H[,,1],"\n")
  cat("==========================================================\n")
  cat("The time-invariance G matrix for the system equation is:\n")
  temp=SSModel_subject$T[,,1]
  colnames(temp)=NULL
  print(temp)
  cat("==========================================================\n")
  cat("The variance matrix for the system equation is:\n")
  SSModel_subject$Q[,,1][is.na(SSModel_subject$Q[,,1])]=-99
  temp2=SSModel_subject$R[,,1] %*% SSModel_subject$Q[,,1] %*% t(SSModel_subject$R[,,1])
  colnames(temp2)=NULL
  print(round(temp2,4))
  cat("==========================================================\n")
  cat("The initial values of states are:\n")
  print(SSModel_subject$a1)
  cat("The initial values of the variance of states are:\n")
  temp3=10^7*SSModel_subject$P1inf+SSModel_subject$P1
  colnames(temp3)=NULL
  print(round(temp3,4))
}
summary.dlmFiltered=function(dlmFilteredsubject){
  if(class(dlmFilteredsubject)!="dlmFiltered"){
    stop("The input subject is not a dlmFiltered subject!")
  }
  
  cat("==========================================================\n")
  cat("The FF matrix is:",dlmFilteredsubject$mod$FF,"\n")
  if(!is.null(dlmFilteredsubject$mod$JFF)){
    cat("The JFF matrix is:",dlmFilteredsubject$mod$JFF,"\n")
    cat("Therefore, the time-invariate component is locates at:",which(dlmFilteredsubject$mod$JFF ==0),"column.\n")
    if((dim(dlmFilteredsubject$mod$F)[2]-sum(dlmFilteredsubject$mod$JFF ==0)) != dim(dlmFilteredsubject$mod$X[,dlmFilteredsubject$mod$JFF])[2]){
      stop("The specification for the dimention of time-varying JFF is wrong!")
    }
    cat("The time-varying component is", dim(dlmFilteredsubject$mod$F)[1], "*", dim(dlmFilteredsubject$mod$F)[2]-sum(dlmFilteredsubject$mod$JFF ==0), "for",dim(dlmFilteredsubject$mod$X[,dlmFilteredsubject$mod$JFF])[1], "timepoints!\n")
    cat("One example for the time-varying F matrix shows as:\n")
    print(dlmFilteredsubject$mod$X[1,dlmFilteredsubject$mod$JFF])
    cat("\n")
  }
  cat("==========================================================\n")
  cat("The variance of the observation equation, V, is:",dlmFilteredsubject$mod$V,"\n")
  cat("==========================================================\n")
  cat("The time-invariance G matrix for the system equation is:\n")
  print(dlmFilteredsubject$mod$GG)
  cat("==========================================================\n")
  cat("The variance matrix for the system equation is:\n")
  print(round(dlmFilteredsubject$mod$W,4))
  cat("==========================================================\n")
  cat("The initial values of states are:\n")
  print(dlmFilteredsubject$mod$m0)
  cat("The initial values of the variance of states are:\n")
  print(dlmFilteredsubject$mod$C0)
  cat("==========================================================\n")
  cat("The filtered estimation of the states E[theta(t)|y_{1:t}] is in m.\n")
  cat("The one-step-head predication of the states E[theta(t)|y_{1:t-1}] is in a.\n")
  cat("The one-step-head predication of the signals E[F*theta(t)|y_{1:t-1}] is in f.\n")
}
refit_changepoint=function(dlmFilteredsubjuct,W_learning_option,W_option_values,
                           V_learning_option=NULL,V_option_values=NULL,V_learning_category=NULL,
                           bandwidth=5,replace_value=10,alternation_W=NULL,alternation_V=NULL){
  if(class(dlmFilteredsubjuct) != "dlmFiltered") stop("Error! input subject should be of dlmFiltered class.")
  # W_learning_option = n_fragment
  # W_learning_option = penalty
  # W_learning_option = given
  if(!(W_learning_option %in% c("n_fragment","penalty","given"))) stop("Error in W_learning_option! choose in ''n_fragment'',''penalty'', and ''given''.")
  # V_learning_option = NULL
  # V_learning_option = n_fragment
  # V_learning_option = penalty
  # V_learning_option = given
  if(!is.null(V_learning_option)){
    if(!(V_learning_option %in% c("n_fragment","penalty","given"))) stop("Error in V_learning_option! choose in ''n_fragment'',''penalty'', and ''given''.")
    if(!(V_learning_category %in% c("mean","variance"))) stop("Error in V_learning_category! choose in ''mean'' or ''variance''.")
  }
  initialmodel=dlmFilteredsubjuct$mod
  n_row=nrow(initialmodel$X)
  n_covariates=ncol(initialmodel$W)
  if(ncol(initialmodel$W) != (ncol(initialmodel$X)+1)) stop("Error in dimention due to intercept term.")
  if(length(W_option_values) != ncol(initialmodel$W)) stop("Error in the length of W_option_values.")
  if(is.null(alternation_W)) alternation_W=1
  if(is.null(alternation_V)) alternation_V=1
  
  # fix JW
  start=ncol(initialmodel$X)+1 # start from the next column of X 
  end=ncol(initialmodel$X)+ncol(initialmodel$W) # end at the start of X plus # (intercept+covariates)
  JW(initialmodel)=diag(start:end)
  
  # add columns of JW to X
  temp=rep(diag(initialmodel$W/alternation_W),n_row)
  if(length(temp)==n_row*n_covariates){
    temp=matrix(temp,nrow=n_row,ncol=n_covariates)
  }else{
    stop("Error in making columns of changing W in X.")
  }
  
  ylab=c("intercept",colnames(initialmodel$X))
  cptm_result=list()
  
  if(W_learning_option=="n_fragment"){
    # find change points given # pices we what
    for(i in 1:n_covariates){
      a=as.numeric(dlmFilteredsubjuct$m[,i])
      #plot(a)
      cptm_stationary <- cpt.mean(a,penalty="Manual", Q=W_option_values[i], method="BinSeg")
      cptm_result[[ylab[i]]]=cptm_stationary 
      cat(ylab[i],":",cpts(cptm_stationary),"\n")
      invisible(plot(cptm_stationary,ylab=ylab[i]))
      selected=unlist(lapply(cpts(cptm_stationary),function(u){(u-bandwidth):(u+bandwidth)}))
      selected=unique(selected[selected>0 & selected<n_row])
      temp[selected,i]=replace_value
    }
  }else if(W_learning_option=="penalty"){
    # find change points given penalty
    for(i in 1:n_covariates){
      a=as.numeric(dlmFilteredsubjuct$m[,i])
      cptm_stationary <- cpt.mean(a,penalty="Manual", pen.value =W_option_values[i], method="PELT")
      cptm_result[[ylab[i]]]=cptm_stationary 
      invisible(plot(cptm_stationary,ylab=ylab[i]))
      cat(ylab[i],":",cpts(cptm_stationary),"\n")
      selected=unlist(lapply(cpts(cptm_stationary),function(u){(u-bandwidth):(u+bandwidth)}))
      selected=unique(selected[selected>0 & selected<n_row])
      temp[selected,i]=replace_value
    }
  }else if(W_learning_option=="given"){
    temp[unique(unlist(lapply(W_option_values,function(u){(u-bandwidth):(u+bandwidth)}))),]=replace_value
  }
  X(initialmodel)=cbind(initialmodel$X,temp)
  colnames(initialmodel$X)=c(ylab[-1],paste(ylab,"W",sep="_"))
  initialmodel=dlm(initialmodel)
  
  
  
  V(initialmodel)=V(initialmodel)*alternation_V
  if(!is.null(V_learning_option)){
    cat("The V also vary over time, according to the observed outcome!\n")
    JV(initialmodel)=ncol(initialmodel$X)+1
    tempv=rep(initialmodel$V,n_row)
    b=as.numeric(dlmFilteredsubjuct$y)
    ylabv=c(colnames(initialmodel$X),"V")
    
    if("mean" %in% V_learning_category){
      if(V_learning_option=="n_fragment"){
        cptm_stationary <- cpt.mean(b,penalty="Manual", Q=V_option_values, method="BinSeg")
        cat("change points for V in mean:",cpts(cptm_stationary),"\n")
        invisible(plot(cptm_stationary,ylab="change points for V in mean"))
        selectedv=unlist(lapply(cpts(cptm_stationary),function(u){(u-bandwidth):(u+bandwidth)}))
      }else if(V_learning_option=="penalty"){
        cptm_stationary <- cpt.mean(b,penalty="Manual", pen.value=V_option_values, method="PELT")
        cat("change points for V in mean:",cpts(cptm_stationary),"\n")
        invisible(plot(cptm_stationary,ylab="change points for V in mean"))
        selectedv=unlist(lapply(cpts(cptm_stationary),function(u){(u-bandwidth):(u+bandwidth)}))
      }else if(V_learning_option=="given"){
        cat("change points for V in mean:",V_option_values,"\n")
        selectedv=unlist(lapply(V_option_values,function(u){(u-bandwidth):(u+bandwidth)}))
      }
      selectedv=unique(selected[selectedv>0 & selectedv<n_row])
      tempv[selectedv]=replace_value
    }
    
    if("variance" %in% V_learning_category){
      if(V_learning_option=="n_fragment"){
        cptm_variance <- cpt.var(b,penalty="Manual", Q=V_option_values, method="BinSeg")
        cat("change points for V in variance:",cpts(cptm_variance),"\n")
        invisible(plot(cptm_variance,ylab="change points for V in variance"))
        cptm_result[["V"]]=cptm_stationary 
        selectedv=unlist(lapply(cpts(cptm_variance),function(u){(u-bandwidth):(u+bandwidth)}))
      }else if(V_learning_option=="penalty"){
        cptm_variance <- cpt.var(b,penalty="Manual", pen.value=V_option_values, method="PELT")
        cat("change points for V in variance:",cpts(cptm_variance),"\n")
        invisible(plot(cptm_variance,ylab="change points for V in variance"))
        cptm_result[["V"]]=cptm_stationary 
        selectedv=unlist(lapply(cpts(cptm_variance),function(u){(u-bandwidth):(u+bandwidth)}))
      }else if(V_learning_option=="given"){
        cat("change points for V in variance:",V_option_values,"\n")
        selectedv=unlist(lapply(V_option_values,function(u){(u-bandwidth):(u+bandwidth)}))
      }
      selectedv=unique(selected[selectedv>0 & selectedv<n_row])
      tempv[selectedv]=replace_value
    }
    
    X(initialmodel)=cbind(initialmodel$X,tempv)
    colnames(initialmodel$X)=ylabv
    initialmodel=dlm(initialmodel)
  }
  
  
  dlmout=dlmFilter(dlmFilteredsubjuct$y,initialmodel)
  return(list(dlmFiltered_out=dlmout,changepoint_plots=cptm_result))
  
}

periodic.tx=function(n,tx_period,p){
  x=rep(NA,n);
  t=1
  while(t<(n+1)){
    x[t:min((t+tx_period-1),n)]=(runif(1)<p)*1
    t=t+tx_period 
  }
  return(x)
}

{
# shift.by.time=function(a,shift,name=NULL){
#   m=length(a)
#   if(shift==1){
#     result=c(NA,a[1:(m-1)])
#   }else if(shift>1){
#     result=matrix(NA,nrow=m,ncol=shift)
#     result[,1]=c(NA,a[1:(m-1)])
#     for(i in 2:shift){
#       result[,i]=c(NA,result[1:(m-1),i-1])
#     }
#     colnames(result)=paste(name,1:shift,sep="_")
#   }
#   return(as.data.frame(result))
# }
} # already replaced by another function
{
# get.dlmoutput=function(n_variables,model,data){
#   summary.SSModel(model)
#   model_fit=fitSSM(model, inits =rep(0,n_variables+1), method = "BFGS") # n_variables + one NA in H
#   cat(blue("Initial fit for the NA parameters complete!\n"))
#   
#   X=data[,grab_X(model)]
#   model_2=function(u){
#     dlmModReg(X,dV=exp(u[1]),dW=exp(u[2:(n_variables+1)])) #  coefficients + baseline
#   }
#   init.parm=log(c(model_fit$model$H,diag(as.matrix(model_fit$model$Q[,,1]))))
#   model_fit_2=dlmMLE(model$y,init.parm,model_2) # estimate V and W
#   cat(blue("Final fit for the NA parameters complete!\n"))
#   model_fitted_2=model_2(model_fit_2$par)
#   model_out=dlmFilter(model$y,model_fitted_2)
#   summary.dlmFiltered(model_out)
#   return(model_out)
# }
}# this function is not truly necessary, as both package gives the same result
######################################################################
###               weekly concentrated analysis                     ###
######################################################################
average_over_every_n_days=function(data,n_group,na.rm=T){
  average_over_every_n_days_each=function(test,n_group){
    rowMeans(matrix(c(test,rep(NA, n_group-(length(test) %% n_group))),ncol=n_group,byrow=T),na.rm = na.rm)
  }
  
  if(!is.data.frame(data) & !is.matrix(data)) stop("The input data should be matrix or dataframe.")
  if(n_group%%1 != 0) stop("n_group should be a integer")
  if(n_group<2) stop("n_group should be at least 2 days")
  
  nrow=ceiling(nrow(data)/n_group)
  result=as.data.frame(matrix(NA,nrow=nrow,ncol=ncol(data)))
  for(i in 1:ncol(data)){
    temp=average_over_every_n_days_each(data[,i],n_group)
    if(length(temp)==nrow){
      result[,i]=average_over_every_n_days_each(data[,i],n_group)
    }else{
      stop("The length of averaged variables are wrong")
    }
  }
  if(!is.null(colnames(data))){
    colnames(result)=colnames(data)
  }
  return(result)
}

#########################################################
##               plot dlm or KFS objects               ##
#########################################################
{
# plot_one_coefficient, as it plots	one-step-ahead predictions of states but not filtered estimates of states
# plot_one_coefficient=function(out,k,truth,ylab){
#   plot(1:nrow(out$a),out$a[,k],
#        ylab=ylab,xlab="time",type="l",cex.lab=1)
#   abline(h=truth,col="red")
#   lines(out$a[,k]+1.65*sqrt(out$P[k,k,]),col="blue")
#   lines(out$a[,k]-1.65*sqrt(out$P[k,k,]),col="blue")
# }
} # wrong
plot.KFS_each=function(KFS_out,k,range){
  sd_temp=sqrt(KFS_out$Ptt[k,k,])
  plot(KFS_out$att[range,k],ylab=colnames(KFS_out$model$T)[k],type="l")
  lines(KFS_out$att[range,k]+1.65*sd_temp[range],col="blue")
  lines(KFS_out$att[range,k]-1.65*sd_temp[range],col="blue")
  abline(h=0,lty=3,col="red")
  cat(colnames(KFS_out$att)[k],": ",KFS_out$att[max(range),k]," (",KFS_out$att[max(range),k]-1.65*sd_temp[max(range)],",",KFS_out$att[max(range),k]+1.65*sd_temp[max(range)],").\n",sep="")
}
plot.KFS=function(KFS_out,selection=NULL,range=NULL){
  if(class(KFS_out)!="KFS"){stop("The result is not a KFS object.")}
  if(is.null(selection)){selection=1:ncol(KFS_out$att)}
  if(is.null(range)){range=1:nrow(KFS_out$att)}
  if(!all(selection %in% 1:ncol(KFS_out$att))){stop("The selected column overflow the #variables")}
  for(id in selection){
    plot.KFS_each(KFS_out,id,range=range)
  }
}
plot.dlmFiltered_each=function(dlmFiltered_subject,k,benchmark,option,range=NULL,plotrange=2){
  if(class(dlmFiltered_subject)!="dlmFiltered"){
    stop("The input should be of class ''dlmFiltered''.")
  }
  if(!(option %in% c("filtered_state", "one-step-ahead_state"))){
    stop("The option for plots could only be [filtered_state], [one-step-ahead_state].") 
  }
  
  if(is.null(range)) range=1:nrow(dlmFiltered_subject$m)
  if(option=="filtered_state"){
    temp=dlmFiltered_subject$m[range,k]
    sd=unlist(lapply(range,function(i){sqrt(diag(dlmSvd2var(dlmFiltered_subject$U.C[[i]],dlmFiltered_subject$D.C[i,]))[k])}))
  }else if(option=="one-step-ahead_state"){
    temp=dlmFiltered_subject$a[,k]
    sd=unlist(lapply(1:nrow(dlmFiltered_subject$a),function(i){sqrt(diag(dlmSvd2var(dlmFiltered_subject$U.R[[i]],dlmFiltered_subject$D.R[i,]))[k])}))
  }
  
  ylab=c("intercept",dimnames(dlmFiltered_subject$mod$X)[[2]])[k]
  
  width=(max(temp)-min(temp))*plotrange
  rule1=range(max(temp)+width,min(temp)-width)
  rule2=range(c(temp+1.65*sd,temp-1.65*sd))
  yrange=range(min(max(rule1),max(rule2)),max(min(rule1),min(rule2)))
  
  plot(temp,type="l",ylab=ylab,xlab="time",cex.lab=1,bty = "n",main=option,xaxt='n',ylim=yrange)
  axis(side=1,at=1:length(temp),labels=as.character(range-1))
  abline(h=benchmark,col="red")
  lines(temp+1.65*sd,col="brown",lty=2)
  lines(temp-1.65*sd,col="brown",lty=2)
}
{
  # plot.dlmFiltered_each=function(dlmFiltered_subject,k,benchmark,option){
  #   if(class(dlmFiltered_subject)!="dlmFiltered"){
  #     stop("The input should be of class ''dlmFiltered''.")
  #   }
  #   if(!(option %in% c("filtered_state", "one-step-ahead_state"))){
  #     stop("The option for plots could only be [filtered_state], [one-step-ahead_state].") 
  #   }
  #   
  #   if(option=="filtered_state"){
  #     temp=dlmFiltered_subject$m[,k]
  #     sd=unlist(lapply(1:nrow(dlmFiltered_subject$m),function(i){sqrt(diag(dlmSvd2var(dlmFiltered_subject$U.C[[i]],dlmFiltered_subject$D.C[i,]))[k])}))
  #   }else if(option=="one-step-ahead_state"){
  #     temp=dlmFiltered_subject$a[,k]
  #     sd=unlist(lapply(1:nrow(dlmFiltered_subject$a),function(i){sqrt(diag(dlmSvd2var(dlmFiltered_subject$U.R[[i]],dlmFiltered_subject$D.R[i,]))[k])}))
  #   }
  #   
  #   ylab=c("intercept",dimnames(dlmFiltered_subject$mod$X)[[2]])[k]
  #   plot(temp,type="l",ylab=ylab,xlab="time",cex.lab=1,bty = "n",main=option)
  #   abline(h=benchmark,col="red")
  #   lines(temp+1.65*sd,col="brown",lty=2)
  #   lines(temp-1.65*sd,col="brown",lty=2)
  # }
}
plot.dlmSmoothed_each=function(dlmSmoothed_subject,k,benchmark,option="smoothed_state",ylab,range=NULL,plotrange=2){
  if(class(dlmSmoothed_subject)!="list"){
    stop("The input should be of class ''dlmFiltered''.")
  }
  if(option != "smoothed_state"){
    stop("The option for plots could only be [smoothed_state].") 
  }
  if(is.null(range)) range=1:nrow(dlmFiltered_subject$m)
  
  temp=dlmSmoothed_subject$s[range,k]
  sd=unlist(lapply(range,function(i){sqrt(diag(dlmSvd2var(dlmSmoothed_subject$U.S[[i]],dlmSmoothed_subject$D.S[i,]))[k])}))
  ylab=c("intercept",dimnames(dlmFiltered_subject$mod$X)[[2]])[k]
  
  width=(max(temp)-min(temp))*plotrange
  rule1=range(max(temp)+width,min(temp)-width)
  rule2=range(c(temp+1.65*sd,temp-1.65*sd))
  yrange=range(min(max(rule1),max(rule2)),max(min(rule1),min(rule2)))
  
  plot(temp,type="l",ylab=ylab,xlab="time",cex.lab=1,bty = "n",main=option,ylim=yrange,xaxt='n')
  axis(side=1,at=1:length(temp),labels=as.character(range-1))
  abline(h=benchmark,col="red")
  lines(temp+1.65*sd,col="brown",lty=2)
  lines(temp-1.65*sd,col="brown",lty=2)
  
}
{
  # plot.dlmSmoothed_each=function(dlmSmoothed_subject,k,benchmark,option="smoothed_state",ylab){
  #   if(class(dlmSmoothed_subject)!="list"){
  #     stop("The input should be of class ''dlmFiltered''.")
  #   }
  #   if(option != "smoothed_state"){
  #     stop("The option for plots could only be [smoothed_state].") 
  #   }
  #   temp=dlmSmoothed_subject$s[,k]
  #   sd=unlist(lapply(1:nrow(dlmSmoothed_subject$s),function(i){sqrt(diag(dlmSvd2var(dlmSmoothed_subject$U.S[[i]],dlmSmoothed_subject$D.S[i,]))[k])}))
  #   range(c(temp+1.65*sd,temp-1.65*sd))
  #   
  #   plot(temp,type="l",ylab=ylab,xlab="time",cex.lab=1,bty = "n",main=option,ylim=range(c(temp+1.65*sd,temp-1.65*sd)))
  #   abline(h=benchmark,col="red")
  #   lines(temp+1.65*sd,col="brown",lty=2)
  #   lines(temp-1.65*sd,col="brown",lty=2)
  # }
  # plot.dlm=function(dlmFiltered_subject,benchmark,option){
  #   if(class(dlmFiltered_subject)!="dlmFiltered"){
  #     stop("The input should be of class ''dlmFiltered''.")
  #   }
  #   if(!all(option %in% c("filtered_state", "one-step-ahead_state","one-step-ahead_signal","smoothed_state"))){
  #     stop("The option for plots could only be [filtered_state], [one-step-ahead_state], [one-step-ahead_signal], [smoothed_state].") 
  #   }
  #   
  #   if("one-step-ahead_signal" %in% option){
  #     par(mfrow=c(1,1))
  #     plot(dlmFiltered_subject$y,type="l",bty = "n",xlab="time",ylab="Outcome",main="Observed versus One-step-ahead signal",)
  #     points(dlmFiltered_subject$f,type="l",col="blue",lty=2)
  #   }
  #   if("filtered_state" %in% option){
  #     par(mfrow=c(2,1))
  #     for(k in 1:length(dlmFiltered_subject$mod$m0)){
  #       plot.dlmFiltered_each(dlmFiltered_subject,k,benchmark[k],option="filtered_state")
  #     }
  #   }
  #   if("one-step-ahead_state" %in% option){
  #     par(mfrow=c(2,1))
  #     for(k in 1:length(dlmFiltered_subject$mod$m0)){
  #       plot.dlmFiltered_each(dlmFiltered_subject,k,benchmark[k],option="one-step-ahead_state")
  #     }
  #   }
  #   if("smoothed_state" %in% option){
  #     dlmSmoothed_subject=dlmSmooth(dlmFiltered_subject)
  #     ylab=c("intercept",dimnames(dlmFiltered_subject$mod$X)[[2]])
  #     par(mfrow=c(2,1))
  #     for(k in 1:length(dlmFiltered_subject$mod$m0)){
  #       plot.dlmSmoothed_each(dlmSmoothed_subject,k,benchmark[k],option="smoothed_state",ylab=ylab[k])
  #     }
  #   }
  #   
  # }
}
plot.dlm=function(dlmFiltered_subject,benchmark,option=NULL,range=NULL,selection=NULL){
  # option: [filtered_state], [one-step-ahead_state], [one-step-ahead_signal], [smoothed_state]
  # range: which range on x axis to plot
  # selection: select a variable to plot
  
  if(class(dlmFiltered_subject)!="dlmFiltered"){
    stop("The input should be of class ''dlmFiltered''.")
  }
  
  if(is.null(selection)){
    if(is.null(option)) stop("''option'' is required!")
    if(!all(option %in% c("filtered_state", "one-step-ahead_state","one-step-ahead_signal","smoothed_state"))){
      stop("The option for plots could only be [filtered_state], [one-step-ahead_state], [one-step-ahead_signal], [smoothed_state].") 
    }
    if("one-step-ahead_signal" %in% option){
      par(mfrow=c(1,1))
      if(is.null(range)){
        signal_range=1:nrow(dlmFiltered_subject$y)
      }else{
        signal_range=range
      }
      drawrange=range(c(dlmFiltered_subject$y[signal_range],dlmFiltered_subject$f[signal_range]))
      plot(dlmFiltered_subject$y[signal_range],type="l",bty = "n",ylim=drawrange,xlab="time",ylab="Outcome",main="Observed versus One-step-ahead signal",xaxt='n')
      axis(side=1,at=1:length(dlmFiltered_subject$y[signal_range]),labels=as.character(signal_range))
      points(dlmFiltered_subject$f[signal_range],type="l",col="blue",lty=2)
    }
    if("filtered_state" %in% option){
      par(mfrow=c(2,1))
      for(k in 1:length(dlmFiltered_subject$mod$m0)){
        plot.dlmFiltered_each(dlmFiltered_subject,k,benchmark[k],option="filtered_state",range=range)
      }
    }
    if("one-step-ahead_state" %in% option){
      par(mfrow=c(2,1))
      for(k in 1:length(dlmFiltered_subject$mod$m0)){
        plot.dlmFiltered_each(dlmFiltered_subject,k,benchmark[k],option="one-step-ahead_state",range=range)
      }
    }
    if("smoothed_state" %in% option){
      dlmSmoothed_subject=dlmSmooth(dlmFiltered_subject)
      ylab=c("intercept",dimnames(dlmFiltered_subject$mod$X)[[2]])
      par(mfrow=c(2,1))
      for(k in 1:length(dlmFiltered_subject$mod$m0)){
        plot.dlmSmoothed_each(dlmSmoothed_subject,k,benchmark[k],option="smoothed_state",ylab=ylab[k],range=range)
      }
    }
  }else if(all(selection %in% c("intercept",dimnames(dlmFiltered_subject$mod$X)[[2]]))){
    if(!is.null(option)) stop("''option'' is not necessary!")
    cat("Plot prediction, state estimation, exposure for ",paste(selection,collapse = ", "),".",sep="")
    ylab=c("intercept",dimnames(dlmFiltered_subject$mod$X)[[2]])
    for(i in 1:length(selection)){
      id=which(ylab==selection[i])
      par(mfrow=c(3,1))
      # prediction
      if(is.null(range)){
        signal_range=1:nrow(dlmFiltered_subject$y)
      }else{
        signal_range=range
      }
      drawrange=range(c(dlmFiltered_subject$y[signal_range],dlmFiltered_subject$f[signal_range]))
      plot(dlmFiltered_subject$y[signal_range],type="l",bty = "n",ylim=drawrange,xlab="time",ylab="Outcome",main="Observed versus One-step-ahead signal",xaxt='n')
      axis(side=1,at=1:length(dlmFiltered_subject$y[signal_range]),labels=as.character(signal_range))
      points(dlmFiltered_subject$f[signal_range],type="l",col="blue",lty=2)
      # state estimation
      plot.dlmFiltered_each(dlmFiltered_subject,id,benchmark[i],option="filtered_state",range=range)
      # exposure
      if(id!=1){
        plot(dlmFiltered_subject$mod$X[signal_range,id-1],type="l",xlab="time",ylab=selection[i],main="exposure/covariate",bty = "n",xaxt='n')
        axis(side=1,at=1:length(dlmFiltered_subject$mod$X[signal_range,id-1]),labels=as.character(signal_range))
      }
    }
  }else{
    stop("Error in ''selection''.")
  }
  
}

# compare plots from "KFS", "dlmFiltered" and "dlmSmoothed"
# contents for comparisions to choose: [filtered_state], [one-step-ahead_state], [one-step-ahead_signal], [smoothed_state]
# KFS: [filtered_state], [one-step-ahead_state], [smoothed_state]
# dlmFiltered: [filtered_state], [one-step-ahead_state]
# dlmSmoothed: [smoothed_state]
# 
#   KFS versus KFS: [filtered_state], [one-step-ahead_state], [smoothed_state]
#                   [filtered_vs_one-step-ahead_state], [filtered_vs_smoothed_state], [one-step-ahead_vs_smoothed_state]
#   KFS versus dlmFiltered: [filtered_state], [one-step-ahead_state]
#   dlmFiltered versus dlmFiltered: [filtered_state], [one-step-ahead_state], [one-step-ahead_signal]
#   KFS versus dlmSmoothed: [smoothed_state]
#   dlmFiltered versus dlmSmoothed: [smoothed_state]
#   dlmSmoothed versus dlmSmoothed: [smoothed_state]
compare.plots=function(subject1,subject2,option){
  # check the class of subjects
  if(!all(c(class(subject2),class(subject1)) %in% c("KFS","dlmFiltered","list"))){
    stop("Error in the class of subjects -> must be ''KFS'',''dlmFiltered'',''list''.")
  }
  
  if(class(subject1)=="KFS" & class(subject2)=="KFS"){
    cat("Plots comparisions between two KFS subjects.\n")
    
    # KFS: [filtered_state], [one-step-ahead_state], [smoothed_state]
    if(!(option %in% c("filtered_state","one-step-ahead_state","smoothed_state",
                       "filtered_vs_one-step-ahead_state","filtered_vs_smoothed_state","one-step-ahead_vs_smoothed_state"))){
      stop("Error in KFS option -> must be [filtered_state], [one-step-ahead_state], [smoothed_state],[filtered_vs_one-step-ahead_state], [filtered_vs_smoothed_state], [one-step-ahead_vs_smoothed_state].")
    }
    if(option %in% c("filtered_vs_one-step-ahead_state","filtered_vs_smoothed_state","one-step-ahead_vs_smoothed_state")){
      cat("Caution: the order in option applies to the order of the subjects!\n")
    }
    if(subject1$dims$m!=subject2$dims$m){
      stop("The number of covariates used by 2 subjects don't match.\n")
    }
    
    # check the dimention of covariates
    covariate1=gsub("\\(Intercept.","intercept",gsub("\\)",".",gsub("factor\\(*","",names(subject1$model$Z[,,1]))))
    covariate2=gsub("\\(Intercept.","intercept",gsub("\\)",".",gsub("factor\\(*","",names(subject2$model$Z[,,1]))))
    cat("The covariates used in the first KFS subject are:",paste(covariate1,collapse = ", "),"\n")
    cat("The covariates used in the second KFS subject are:",paste(covariate2,collapse = ", "),"\n")
    if(!all(covariate1==covariate2)){
      stop("The covariates used for the two models differ!")
    }else{
      cat("The covariates used by two models match.\n")
    }
    
    par(mfrow=c(2,1))
    for(i in 1:subject1$dims$m){
      if(option=="filtered_state"){
        plot(subject1$att[,i],ylab=covariate1[i],main="filtered state (KFS1)")
        plot(subject2$att[,i],ylab=covariate2[i],main="filtered state (KFS2)")
      }else if(option=="one-step-ahead_state"){
        plot(subject1$a[,i],ylab=covariate1[i],main="one-step-ahead state (KFS1)")
        plot(subject2$a[,i],ylab=covariate2[i],main="one-step-ahead state (KFS2)")
      }else if(option=="smoothed_state"){
        plot(subject1$alphahat[,i],ylab=covariate1[i],main="smoothed state (KFS1)")
        plot(subject2$alphahat[,i],ylab=covariate2[i],main="smoothed state (KFS2)")
      }else if(option=="filtered_vs_one-step-ahead_state"){
        plot(subject1$att[,i],ylab=covariate1[i],main="filtered state (KFS1)")
        plot(subject2$a[,i],ylab=covariate2[i],main="one-step-ahead state (KFS2)")
      }else if(option=="filtered_vs_smoothed_state"){
        plot(subject1$att[,i],ylab=covariate1[i],main="filtered state (KFS1)")
        plot(subject2$alphahat[,i],ylab=covariate2[i],main="smoothed state (KFS2)")
      }else if(option=="one-step-ahead_vs_smoothed_state"){
        plot(subject1$a[,i],ylab=covariate1[i],main="one-step-ahead state (KFS1)")
        plot(subject2$alphahat[,i],ylab=covariate2[i],main="smoothed state (KFS2)")
      }
    }
  }else if(class(subject1)=="KFS" & class(subject2)=="dlmFiltered"){ # ---------- case 2 -------------
    cat("Plots comparisions between a KFS subject and a dlmFiltered subject.\n")
    
    # KFS: [filtered_state], [one-step-ahead_state], [smoothed_state]
    # dlmFiltered: [filtered_state], [one-step-ahead_state]
    if(!(option %in% c("filtered_state","one-step-ahead_state",
                       "filtered_vs_one-step-ahead_state","one-step-ahead_vs_filtered_state","smoothed_vs_filtered_state","smoothed_vs_one-step-ahead_state"))){
      stop("Error in KFS versus dlmFiltered option -> must be [filtered_state], [one-step-ahead_state], [filtered_vs_one-step-ahead_state], [one-step-ahead_vs_filtered_state], [smoothed_vs_filtered_state], [smoothed_vs_one-step-ahead_state].")
    }
    if(subject1$dims$m!=length(subject2$mod$m0)){
      stop("The number of covariates used by 2 subjects don't match.\n")
    }
    
    # check the dimention of covariates
    covariate1=gsub("\\(Intercept.","intercept",gsub("\\)",".",gsub("factor\\(*","",names(subject1$model$Z[,,1]))))
    if(subject2$mod$JFF[,1]==0){
      covariate2=c("intercept",dimnames(subject2$mod$X)[[2]])
    }else{
      covariate2=dimnames(subject2$mod$X)[[2]]
    }
    cat("The covariates used in the first KFS subject are:",paste(covariate1,collapse = ", "),"\n")
    cat("The covariates used in the second dlmFiltered subject are:",paste(covariate2,collapse = ", "),"\n")
    if(!all(covariate1==covariate2)){
      stop("The covariates used for the two models differ!")
    }else{
      cat("The covariates used by two models match.\n")
    }
    
    par(mfrow=c(2,1))
    for(i in 1:subject1$dims$m){
      if(option=="filtered_state"){
        plot(subject1$att[,i],ylab=covariate1[i],main="filtered state (KFS)")
        plot(subject2$m[,i],type="l",ylab=covariate2[i],main="filtered state (dlmFiltered)")
      }else if(option=="one-step-ahead_state"){
        plot(subject1$a[,i],ylab=covariate1[i],main="one-step-ahead state (KFS)")
        plot(subject2$a[,i],type="l",ylab=covariate2[i],main="one-step-ahead state (dlmFiltered)")
      }else if(option=="filtered_vs_one-step-ahead_state"){
        plot(subject1$att[,i],ylab=covariate1[i],main="filtered state (KFS)")
        plot(subject2$a[,i],type="l",ylab=covariate2[i],main="one-step-ahead state (dlmFiltered)")
      }else if(option=="one-step-ahead_vs_filtered_state"){
        plot(subject1$a[,i],ylab=covariate1[i],main="one-step-ahead state (KFS)")
        plot(subject2$m[,i],type="l",ylab=covariate2[i],main="filtered state (dlmFiltered)")
      }else if(option=="smoothed_vs_filtered_state"){
        plot(subject1$alphahat[,i],ylab=covariate1[i],main="smoothed state (KFS)")
        plot(subject2$m[,i],type="l",ylab=covariate2[i],main="filtered state (dlmFiltered)")
      }else if(option=="smoothed_vs_one-step-ahead_state"){
        plot(subject1$alphahat[,i],ylab=covariate1[i],main="smoothed state (KFS)")
        plot(subject2$a[,i],type="l",ylab=covariate2[i],main="one-step-ahead state (dlmFiltered)")
      }
    }
  }else if(class(subject2)=="KFS" & class(subject1)=="dlmFiltered"){
    stop("Please reverse the position of two subjects. KFS should go first!")
  }else if(class(subject1)=="dlmFiltered" & class(subject2)=="dlmFiltered"){
    cat("Plots comparisions between two dlmFiltered subjects.\n")
    # dlmFiltered: [filtered_state], [one-step-ahead_state]
    if(!(option %in% c("filtered_state","one-step-ahead_state",
                       "filtered_vs_one-step-ahead_state","one-step-ahead_vs_filtered_state"))){
      stop("Error in dlmFiltered versus dlmFiltered option -> must be [filtered_state], [one-step-ahead_state], [filtered_vs_one-step-ahead_state], [one-step-ahead_vs_filtered_state].")
    }
    if(length(subject1$mod$m0)!=length(subject2$mod$m0)){
      stop("The number of covariates used by 2 subjects don't match.\n")
    }
    
    # check the dimention of covariates
    if(subject1$mod$JFF[,1]==0){
      covariate1=c("intercept",dimnames(subject1$mod$X)[[2]][1:(length(subject1$mod$m0)-1)])
    }else{
      covariate1=dimnames(subject1$mod$X)[[2]][1:(length(subject1$mod$m0)-1)]
    }
    if(subject2$mod$JFF[,1]==0){
      covariate2=c("intercept",dimnames(subject2$mod$X)[[2]][1:(length(subject2$mod$m0)-1)])
    }else{
      covariate2=dimnames(subject2$mod$X)[[2]][1:(length(subject2$mod$m0)-1)]
    }
    cat("The covariates used in the first dlmFiltered subject are:",paste(covariate1,collapse = ", "),"\n")
    cat("The covariates used in the second dlmFiltered subject are:",paste(covariate2,collapse = ", "),"\n")
    if(!all(covariate1==covariate2)){
      stop("The covariates used for the two models differ!")
    }else{
      cat("The covariates used by two models match.\n")
    }
    
    par(mfrow=c(2,1))
    for(i in 1:length(subject1$mod$m0)){
      if(option=="filtered_state"){
        plot(subject1$m[,i],type="l",ylab=covariate1[i],main="filtered state (dlmFiltered1)")
        plot(subject2$m[,i],type="l",ylab=covariate2[i],main="filtered state (dlmFiltered2)")
        #points(subject2$m[,i],type="l",col="blue")
      }else if(option=="one-step-ahead_state"){
        plot(subject1$a[,i],type="l",ylab=covariate1[i],main="one-step-ahead state (dlmFiltered1)")
        plot(subject2$a[,i],type="l",ylab=covariate2[i],main="one-step-ahead state (dlmFiltered2)")
      }else if(option=="filtered_vs_one-step-ahead_state"){
        plot(subject1$m[,i],type="l",ylab=covariate1[i],main="filtered state (dlmFiltered1)")
        plot(subject2$a[,i],type="l",ylab=covariate2[i],main="one-step-ahead state (dlmFiltered2)")
      }else if(option=="one-step-ahead_vs_filtered_state"){
        plot(subject1$a[,i],type="l",ylab=covariate1[i],main="one-step-ahead state (dlmFiltered1)")
        plot(subject2$m[,i],type="l",ylab=covariate2[i],main="filtered state (dlmFiltered2)")
      }
    }
  }else if(class(subject1)=="KFS" & class(subject2)=="list"){
    cat("Plots comparisions between a KFS subject and a dlmSmoothed subject.\n")
    
    # KFS: [filtered_state], [one-step-ahead_state], [smoothed_state]
    # dlmSmoothed: [smoothed_state]
    if(!(option %in% c("smoothed_state","filtered_vs_smoothed_state","one-step-ahead_vs_smoothed_state"))){
      stop("Error in KFS versus dlmSmoothed option -> must be [smoothed_state], [filtered_vs_smoothed_state], [one-step-ahead_vs_smoothed_state].")
    }
    if(subject1$dims$m!=dim(subject2$s)[2]){
      stop("The number of covariates used by 2 subjects don't match.\n")
    }
    # check the dimention of covariates
    covariate1=gsub("\\(Intercept.","intercept",gsub("\\)",".",gsub("factor\\(*","",names(subject1$model$Z[,,1]))))
    par(mfrow=c(2,1))
    for(i in 1:subject1$dims$m){
      if(option=="smoothed_state"){
        plot(subject1$alphahat[,i],ylab=covariate1[i],main="smoothed state (KFS)")
        plot(subject2$s[,i],ylab=covariate1[i],type="l",main="smoothed state (dlmSmoothed)")
      }else if(option=="filtered_vs_smoothed_state"){
        plot(subject1$att[,i],ylab=covariate1[i],main="filtered state (KFS)")
        plot(subject2$s[,i],ylab=covariate1[i],type="l",main="smoothed state (dlmSmoothed)")
      }else if(option=="one-step-ahead_vs_smoothed_state"){
        plot(subject1$a[,i],ylab=covariate1[i],main="one-step-ahead state (KFS)")
        plot(subject2$s[,i],ylab=covariate1[i],type="l",main="smoothed state (dlmSmoothed)")
      }
    }
  }else if(class(subject2)=="KFS" & class(subject1)=="list"){
    stop("Please reverse the position of two subjects. KFS should go first!")
  }else if(class(subject1)=="dlmFiltered" & class(subject2)=="list"){
    cat("Plots comparisions between a dlmFiltered subject and a dlmSmoothed subject.\n")
    # dlmFiltered: [filtered_state], [one-step-ahead_state]
    # dlmSmoothed: [smoothed_state]
    if(!(option %in% c("filtered_vs_smoothed_state","one-step-ahead_vs_smoothed_state"))){
      stop("Error in KFS versus dlmSmoothed option -> must be [filtered_vs_smoothed_state], [one-step-ahead_vs_smoothed_state].")
    }
    if(length(subject1$mod$m0)!=dim(subject2$s)[2]){
      stop("The number of covariates used by 2 subjects don't match.\n")
    }
    # check the dimention of covariates
    if(subject1$mod$JFF[,1]==0){
      covariate1=c("intercept",dimnames(subject2$mod$X)[[2]])
    }else{
      covariate1=dimnames(subject2$mod$X)[[2]]
    }
    
    par(mfrow=c(2,1))
    for(i in 1:length(subject1$mod$m0)){
      if(option=="filtered_vs_smoothed_state"){
        plot(subject1$m[,i],type="l",ylab=covariate1[i],main="filtered state (dlmFiltered)")
        plot(subject2$s[,i],ylab=covariate1[i],type="l",main="smoothed state (dlmSmoothed)")
      }else if(option=="one-step-ahead_vs_smoothed_state"){
        plot(subject1$a[,i],type="l",ylab=covariate1[i],main="one-step-ahead state (dlmFiltered1)")
        plot(subject2$s[,i],ylab=covariate1[i],type="l",main="smoothed state (dlmSmoothed)")
      }
    }
  }else if(class(subject2)=="dlmFiltered" & class(subject1)=="list"){
    stop("Please reverse the position of two subjects. dlmFiltered should go first!")
  }else if(class(subject1)=="list" & class(subject2)=="list"){
    cat("Plots comparisions between two dlmSmoothed subjects.\n")
    if(!(option %in% c("smoothed_state"))){
      stop("Error in dlmSmoothed option -> must be [smoothed_state].")
    }
    if(dim(subject1$s)[2]!=dim(subject2$s)[2]){
      stop("The number of covariates used by 2 subjects don't match.\n")
    }
    for(i in 1:subject1$dims$m){
      plot(subject1$s[,i],type="l",main="smoothed state (dlmSmoothed1)")
      plot(subject2$s[,i],type="l",main="smoothed state (dlmSmoothed2)")
    }
  }
  
}



