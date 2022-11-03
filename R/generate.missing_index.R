#' A function for generating
#'
#' @param type three missing type MCAR,MAR,MNAR
#' @param param must be a list and \code{param$data} must be a data.frame
#' @param n the length of original data
#' @param printFlag print intermediate results and plots, default is T
#'
#' @return A list of result, containing the missing index
#' @export
#'
#' @examples \dontrun{generate.missing_index(type = "MNAR", n=length(data$y),param = list(data = data.frame(y=data$y,x=data$x), MNAR.type = "increasing", coeff = c(0.5,-1,0), MNAR.drawplot = c(TRUE, "y")))$missing_index}
#'
#'
generate.missing_index=function(type,param,n,printFlag=T){
  # The function is mainly written by Xinru Wang, edited by Xiaoxuan Cai from 03/17/21 to 03/20/21
  # The general idea is from the page 63 of book <<Flexible Imputation of Missing Data>> by Stef van Buuren

  # This function generate missing index giving three machenisms: MCAR, MAR, MNAR
  # type: [MCAR,MAR,MNAR]
  # n: the length of original data
  # param: must be a list -> param$data must be a data.frame
  # (1) MCAR: <p>
  #     Example: generate.missing_index(type = "MCAR", n=length(data$y), param = list(p=0.5))
  # (2) MAR: (generate missing depending on observed covariates)
  #          <data> of covariates;
  #          <MAR.type> (choose in ["increasing","tail","middle])
  #          <coeff> (#coeff=#covariates+1 for "increasing;
  #                   #coeff=#covariates+1 for "tail";
  #                   #coeff=#covariates+1 for "middle";)
  #         <MMAR.drawplot> = c(T/F,<name of the covariate>)
  #     Example: generate.missing_index(type = "MAR", n=length(data$y),
  #                                     param =list(data = data.frame(y_1=data$y_1),
  #                                                 MAR.type = "increasing",
  #                                                 coeff = c(0.2,-8),
  #                                                 MAR.drawplot = c(TRUE, "y_1"))
  #                                      )
  # (3) MNAR: generate missing depending on unobserved values of y and c/x
  #           <data> unmeasured y or c/x or both
  #           <MNAR.type> (choose in ["increasing","tail","middle])
  #           <coeff> (#coeff=#covariates+1 for "increasing;
  #                  #coeff=#covariates+1 for "tail";
  #                  #coeff=#covariates+1 for "middle";)
  #           <MNAR.drawplot> = c(T/F,<name of the covariate>)
  #     Example: generate.missing_index(type = "MNAR", n=length(data$y),
  #                                     param = ist(data = data.frame(y=data$y,x=data$x),
  #                                                 MNAR.type = "increasing",
  #                                                 coeff = c(0.5,-1,0),
  #                                                 MNAR.drawplot = c(TRUE, "y"))
  #                                     )

  # Pre-checking
  # check parameters: check param$type choose %in%  [MCAR,MAR,MNAR]
  #                   check param$data is data frame
  #                   check the length of coeff agrees with the chooseing missing machenism
  # type=="MCAR", length(list)==1  note: Xinru change to length(list)==1 (p)
  # type=="MAR":  length(list)==4 (data,coeff,MAR.type,MAR.drawplot)
  # type=="MNAR": length(list)==4 (data,coeff,MNAR.type,MNAR.drawplot)

  logistic <- function(x) {exp(x)/(1+exp(x)) }

  ## check type
  if(!type %in% c("MCAR","MAR","MNAR")){
    stop("The missing type shoud choose from ''MCAR'',''MAR'',''MNAR''!\n")
  }
  if(!is.list(param)){
    stop("The parameter type is wrong!\n")
  }

  ## Check the length of parameter is 1 for MCAR
  if(type=="MCAR"){
    if(length(param)!=1)
    { stop("The number of parameters for MCAR shoudl be 1!")}
    # create missingness indicator
    p=param$p
    cat("The missing probability is:",p,"\n")
    missing=(2:n)[as.numeric(runif(n-1)<p)==1]
    if(T){cat(paste(paste("The missing rate for MCAR is", length(missing)/n, sep = " "),"\n",sep=""))}
  }

  if(type=="MAR"){
    # requirement: param$data must be a data frame
    #              param$MAR.type in c("increasing","tail","middle")
    #              param$coeff
    #              param$MAR.drawplot = T/F
    if(length(param)!=4){
      stop("We require 4 parameters for MAR.")
    }
    if(!all(names(param) == c("data","MAR.type","coeff","MAR.drawplot"))){
      stop("The input parameter must be [data],[MAR.type],[coeff],[MAR.drawplot]")
    }
    if(!is.data.frame(param$data) ){
      stop("The data in param$data must be a data frame.")
    }
    if(nrow(param$data)!=n){
      stop("The length of given dataset for MAR is wrong.\n")
    }
    if(!param$MAR.type %in% c("increasing","tail","middle")){
      stop("The MAR.type must be increasing/tail/middle.")
    }

    coeff <- param$coeff
    dat <- param$data %>% mutate(intercept=1)

    # model is logit(p) = coeff[last]*intercept + X\beta
    # the number of coeff is one more than the number of variables
    if(param$MAR.type=="increasing"){
      if(ncol(param$data) != (length(coeff)-1)){
        stop("The number of coefficient for MAR.type 'increasing' is wrong!")
      }
      if(ncol(dat)!=length(coeff)){
        stop("The dat is not correct.")
      }

      ## missing probability
      ## for one variable logistic(coeff[1]*y_1+coeff[2])
      p <- logistic(as.matrix(dat) %*% coeff)
      missing=(2:n)[rbinom((n-1),size=1,prob=(1-p[-1]))==0]
      if(printFlag){
        cat("The missing indicators for MAR 'increasing' are generated!\n")
      }
    }

    # model: logit(p)= coeff[1] + | coeff[2]y_1 (+coeffi[3] c) + coeff[last] |
    #        the number in coeffi is 2 more than the number of variables
    if(param$MAR.type=="tail"){
      if(ncol(param$data) != (length(coeff)-2)){
        stop("The number of coefficients for MAR.type 'tail' is wrong.")
      }
      if(ncol(dat)!=length(coeff)-1){
        stop("The dat is not correct.")
      }

      ## missing probability
      ## for one variable logistic(coeff[1]+abs(coeff[2]*y_1+coeff[3]))
      p <- logistic(coeff[1]+abs(as.matrix(dat) %*% coeff[-1]))
      missing=(2:n)[rbinom((n-1),size=1,prob=(1-p[-1]))==0]
      if(printFlag){
        cat("The missing indicators for MAR 'tail' are generated!\n")
      }
    }

    # model: logit(p)= coeff[1] - | coeff[2]y_1 (+coeffi[3] c) + coeff[last] |
    #        the number in coeffi is 2 more than the number of variables
    if(param$MAR.type=="middle"){
      if(ncol(param$data) != (length(coeff)-2)){
        stop("The number of coefficients for MAR.type 'middle' is wrong.")
      }
      if(ncol(dat)!=length(coeff)-1){
        stop("The dat is not correct.")
      }
      ## missing probability
      ## for one variable logistic(coeff[1]+abs(coeff[2]*y_1+coeff[3]))
      p <- logistic(coeff[1]-abs(as.matrix(dat) %*% coeff[-1]))
      missing=(2:n)[rbinom((n-1),size=1,prob=(1-p[-1]))==0]
      if(printFlag){
        cat("The missing indicators for MAR 'middle' are generated!\n")
      }
    }

    if(printFlag){
      if(param$MAR.drawplot[1]){
        # hist(p)
        for(j in 2:length(param$MAR.drawplot)){
          data_plot=data.frame(x=dat[-1,param$MAR.drawplot[j]],p=p[-1])
          g1 = ggplot(data_plot, aes(x=x, y=p))+
            geom_line()+
            theme_bw()+
            ylab("Missing rate")+
            theme(plot.margin = unit(c(1,1,1,1), "cm"),
                  axis.title.x = element_text(size=15,vjust = -1),
                  axis.title.y = element_text(size=15,vjust=4),
                  axis.text = element_text(size=13))+
            xlab(param$MAR.drawplot[j])
          print(g1)
        }
      }
    }
    cat(paste(paste("The missing rate for MAR is", length(missing)/n, sep = " "),"\n",sep=""))
  }

  if(type=="MNAR"){
    # requirement: param$data must be a data frame
    #              param$MAR.type in c("increasing","tail","middle")
    #              param$coeff
    #              param$MAR.drawplot = T/F
    if(length(param)!=4){
      stop("We require 4 parameters for MAR.")
    }
    if(!all(names(param) == c("data","MNAR.type","coeff","MNAR.drawplot"))){
      stop("The input parameter must be [data],[MNAR.type],[coeff],[MNAR.drawplot]")
    }
    if(!is.data.frame(param$data) ){
      stop("The data in param$data must be a data frame.")
    }
    if(nrow(param$data)!=n){
      stop("The length of given dataset for MAR is wrong.\n")
    }
    if(!param$MNAR.type %in% c("increasing","tail","middle")){
      stop("The MAR.type must be increasing/tail/middle.")
    }

    coeff <- param$coeff
    dat <- param$data %>% mutate(intercept=1)

    if(param$MNAR.type=="increasing"){
      if(ncol(param$data) != (length(coeff)-1)){
        stop("The number of coefficient for MNAR.type 'increasing' is wrong!")
      }
      if(ncol(dat)!=length(coeff)){
        stop("The dat is not correct.")
      }

      ## missing probability
      ## for one variable logistic(coeff[1]*y_1+coeff[2])
      p <- logistic(as.matrix(dat) %*% coeff)
      missing=(2:n)[rbinom((n-1),size=1,prob=(1-p[-1]))==0]
      if(printFlag){cat("The missing indicators for MNAR 'increasing' are generated!\n")}
    }

    if(param$MNAR.type=="tail"){
      if(ncol(param$data) != (length(coeff)-2)){
        stop("The number of coefficients for MNAR.type 'tail' is wrong.")
      }
      if(ncol(dat)!=length(coeff)-1){
        stop("The dat is not correct.")
      }

      ## missing probability
      ## for one variable logistic(coeff[1]+abs(coeff[2]*y_1+coeff[3]))
      p <- logistic(coeff[1]+abs(as.matrix(dat) %*% coeff[-1]))
      missing=(2:n)[rbinom((n-1),size=1,prob=(1-p[-1]))==0]
      if(printFlag){cat("The missing indicators for MNAR 'tail' are generated!\n")}
    }

    if(param$MNAR.type=="middle"){
      if(ncol(param$data) != (length(coeff)-2)){
        stop("The number of coefficients for MNAR.type 'middle' is wrong.")
      }
      if(ncol(dat)!=length(coeff)-1){
        stop("The dat is not correct.")
      }
      ## missing probability
      ## for one variable logistic(coeff[1]+abs(coeff[2]*y_1+coeff[3]))
      p <- logistic(coeff[1]-abs(as.matrix(dat) %*% coeff[-1]))
      missing=(2:n)[rbinom((n-1),size=1,prob=(1-p[-1]))==0]
      if(printFlag){cat("The missing indicators for MNAR 'middle' are generated!\n")}
    }

    if(printFlag){
      if(param$MNAR.drawplot[1])
      {
        # hist(p)
        for(m in 2:length(param$MNAR.drawplot)){
          data_plot=data.frame(x=dat[-1,param$MNAR.drawplot[m]],p=p[-1])
          g1 = ggplot(data_plot, aes(x=x, y=p))+
            geom_line()+
            theme_bw()+
            ylab("Missing rates")+
            theme(plot.margin = unit(c(1,1,1,1), "cm"),
                  axis.title.x = element_text(size=15,vjust = -1),
                  axis.title.y = element_text(size=15,vjust=4),
                  axis.text = element_text(size=13))+
            xlab(param$MNAR.drawplot[m])
          print(g1)
        }
      }
    }
    cat(paste(paste("The missing rate for MNAR is", length(missing)/n, sep = " "),"\n",sep=""))
  }

  return(list(missing_index=missing,missing_rates=p))
}
