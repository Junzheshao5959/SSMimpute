## code to prepare `data_stationary` dataset goes here
library(lubridate)
{
  extend=function(parameter,periods){
    if(length(parameter)==1){
      result=rep(parameter,sum(periods))
    }else if(length(parameter)==length(periods)){
      result=rep(parameter,periods)
    }else if(length(parameter)==sum(periods)){
      result=parameter
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


  # generate time series with tangled y, x and c
  #   current y ~ previous y + current x + multiple past x + current c + interaction
  #   current x ~ previous y + previous x + current c
  #   current c ~ previous y + previous x + previous c
  sim_y_x_single_c=function(given.c=NULL,input.c=NULL,parameters_for_c=NULL,model_for_c=NULL,initial.c=NULL,
                            given.x=F,input.x=NULL,parameters_for_x=NULL,model_for_x=NULL,initial.x=NULL,type_for_x,
                            given.y=F,input.y=NULL,parameters_for_y=NULL,model_for_y=NULL,initial.y=NULL,lag_x_on_y=NULL,interaction_pairs=NULL,n,printFlag=F){
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
  }
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
}

periods = c(1100)

# parameters for c
parameters_for_c=list(baseline=extend(5,periods),
                      c=extend(0.25,periods),
                      sd=extend(sqrt(1),periods))
model_for_c=c("baseline","c")
# parameters for x
parameters_for_x=list(baseline=extend(5,periods),
                      x=extend(0.5,periods),
                      sd=extend(sqrt(1),periods))
model_for_x=c("baseline","x")
# parameters for y
lag_x_on_y=2;treatment_effects=matrix(c(-1.5,-0.5),ncol=lag_x_on_y,byrow=T)
model_for_y=c("baseline","c","x","y")
type_for_x="continuous"
parameters_for_y=list(baseline=extend(40,periods),
                      y=extend(0.5,periods),
                      x=extend_into_matrix(treatment_effects,periods),
                      c=extend(-1,periods),
                      sd=extend(sqrt(1),periods))
colnames(parameters_for_y$x)=paste("t",0:(ncol(parameters_for_y$x)-1),sep="-")

sim_parameters=list(given.c=F,input.c=NULL,parameters_for_c=parameters_for_c,
                    model_for_c=model_for_c,initial.c=NULL,
                    given.x=F,input.x=NULL,parameters_for_x=parameters_for_x,
                    model_for_x=model_for_x,initial.x=NULL,type_for_x=type_for_x,
                    given.y=F,input.y=NULL,parameters_for_y=parameters_for_y,
                    model_for_y=model_for_y,initial.y=NULL,
                    lag_x_on_y=lag_x_on_y,
                    interaction_pairs=NULL,n=sum(periods))

simulation = sim_y_x_single_c(given.c = sim_parameters$given.c,
                              input.c = sim_parameters$input.c,
                              parameters_for_c = sim_parameters$parameters_for_c,
                              model_for_c = sim_parameters$model_for_c,
                              initial.c = sim_parameters$initial.c,
                              given.x = sim_parameters$given.x,
                              input.x = sim_parameters$input.x,
                              parameters_for_x = sim_parameters$parameters_for_x,
                              model_for_x = sim_parameters$model_for_x,
                              initial.x = sim_parameters$initial.x,
                              type_for_x = sim_parameters$type_for_x,
                              given.y = sim_parameters$given.y,
                              input.y = sim_parameters$input.y,
                              parameters_for_y = sim_parameters$parameters_for_y,
                              model_for_y = sim_parameters$model_for_y,
                              initial.y = sim_parameters$initial.y,
                              lag_x_on_y = sim_parameters$lag_x_on_y,
                              interaction_pairs = sim_parameters$interaction_pairs,
                              n = sim_parameters$n,
                              printFlag = T)
simulation$y=simulation$y[-(1:100)]
simulation$x=simulation$x[-(1:100)]
simulation$c=simulation$c[-(1:100)]
data_stationary = data.frame(Date=Sys.Date()-days(length(simulation$y):1),y=simulation$y,x=simulation$x,c=simulation$c)
usethis::use_data(data_stationary, overwrite = TRUE)
