# correct version: 2021/11/09 11:48 pm

#####################################################
##          small function in preparation         ###
#####################################################
# [LOCF] and [LOCF.default] directly copied from R package "DescTools", avoiding downloading the whole R package "DescTools"
LOCF <- function(x) UseMethod("LOCF")
LOCF.default <- function(x) {
  l <- !is.na(x)
  rep(c(NA, x[l]), diff(c(1L, which(l), length(x) + 1L)))
  
}
# [nocb] directly copied from a R package (forgot which R package)
nocb <- function(x, rule=0.75) {
  for (row in 1:nrow(x)) {
    if (sum(is.na(x[row,]))/length(x[row,]) <= (1-rule)) {
      v <- !is.na(x[row,])
      x[row,] <- c(NA, x[row,][v])[cumsum(v)+as.numeric(!v)+1]
    }
  }
  return(x)
}
# [na_kalman] directly copied from R package "imputeTS" directly, avoiding downloading the whole R package "imputeTS"
na_kalman <- function(x, model = "StructTS", smooth = TRUE, nit = -1, maxgap = Inf) {
  data <- x
  
  
  #----------------------------------------------------------
  # Mulivariate Input
  # The next 20 lines are just for checking and handling multivariate input.
  #----------------------------------------------------------
  
  # Check if the input is multivariate
  if (!is.null(dim(data)[2]) && dim(data)[2] > 1) {
    # Go through columns and impute them by calling this function with univariate input
    for (i in 1:dim(data)[2]) {
      if (!anyNA(data[, i])) {
        next
      }
      # if imputing a column does not work - mostly because it is not numeric - the column is left unchanged
      tryCatch(data[, i] <- na_kalman(data[, i], model, smooth, nit, maxgap,...), error = function(cond) {
        warning(paste("imputeTS: No imputation performed for column", i, "because of this", cond), call. = FALSE)
      })
    }
    return(data)
  }
  
  
  #----------------------------------------------------------
  # Univariate Input
  # All relveant imputation / pre- postprocessing  code is within this part
  #----------------------------------------------------------
  
  else {
    missindx <- is.na(data)
    
    ##
    ## 1. Input Check and Transformation
    ##
    
    
    # 1.1 Check if NAs are present
    if (!anyNA(data)) {
      return(data)
    }
    
    # 1.2 special handling data types
    if (any(class(data) == "tbl")) {
      data <- as.vector(as.data.frame(data)[, 1])
    }
    
    # 1.3 Check for algorithm specific minimum amount of non-NA values
    if (sum(!missindx) < 3) {
      stop("Input data needs at least 3 non-NA data point for applying na_kalman")
    }
    
    
    # 1.4 Checks and corrections for wrong data dimension
    
    # Check if input dimensionality is not as expected
    if (!is.null(dim(data)[2]) && !dim(data)[2] == 1) {
      stop("Wrong input type for parameter x")
    }
    
    # Altering multivariate objects with 1 column (which are essentially
    # univariate) to be dim = NULL
    if (!is.null(dim(data)[2])) {
      data <- data[, 1]
    }
    
    # 1.5 Check if input is numeric
    if (!is.numeric(data)) {
      stop("Input x is not numeric")
    }
    
    # 1.6 Check if type of parameter smooth is correct
    if (!is.logical(smooth)) {
      stop("Parameter smooth must be of type logical ( TRUE / FALSE)")
    }
    
    # 1.7 Transformation to numeric as 'int' can't be given to KalmanRun
    data[1:length(data)] <- as.numeric(data)
    
    ##
    ## End Input Check and Transformation
    ##
    
    
    ##
    ## 2. Imputation Code
    ##
    
    # 2.1 Selection of state space model
    
    # State space representation of a arima model
    if (model[1] == "auto.arima") {
      mod <- forecast::auto.arima(data)$model
    }else if (model[1] == "StructTS") {
      # Fallback, because for StructTS first value is not allowed to be NA
      if (is.na(data[1])) {
        data[1] <- nocb(matrix(data,nrow=1), rule=0.2)[1]
      }
      mod <- stats::StructTS(data)$model0
    }else {
      mod <- model
      if (length(mod) < 7) {
        stop("Parameter model has either to be \"StructTS\"/\"auto.arima\" or a user supplied model in
            form of a list with at least components T, Z, h , V, a, P, Pn specified")
      }
      
      if (is.null(mod$Z)) {
        stop("Something is wrong with the user supplied model. Either choose \"auto.arima\" or \"StructTS\"
             or supply a state space model with at least components T, Z, h , V, a, P, Pn as specified
             under Details on help page for KalmanLike")
      }
    }
    
    
    # 2.2 Selection if KalmanSmooth or KalmanRun
    
    if (smooth == TRUE) {
      kal <- stats::KalmanSmooth(data, mod, nit)
      erg <- kal$smooth # for kalmanSmooth
    }else {
      kal <- stats::KalmanRun(data, mod, nit)
      erg <- kal$states # for kalmanrun
    }
    
    # Check if everything is right with the model
    if (dim(erg)[2] != length(mod$Z)) {
      stop("Error with number of components $Z")
    }
    
    # 2.3 Getting Results
    
    # Out of all components in $states or$smooth only the ones
    # which have 1 or -1 in $Z are in the model
    # Therefore matrix multiplication is done
    karima <- erg[missindx, , drop = FALSE] %*% as.matrix(mod$Z)
    
    # Add imputations to the initial dataset
    data[missindx] <- karima
    
    ##
    ## End Imputation Code
    ##
    
    
    ##
    ## 3. Post Processing
    ##
    
    # 3.1 Check for Maxgap option
    
    # If maxgap = Inf then do nothing and when maxgap is lower than 0
    if (is.finite(maxgap) && maxgap >= 0) {
      
      # Get logical vector of the time series via is.na() and then get the
      # run-length encoding of it. The run-length encoding describes how long
      # the runs of FALSE and TRUE are
      rlencoding <- rle(is.na(x))
      
      # Runs smaller than maxgap (which shall still be imputed) are set FALSE
      rlencoding$values[rlencoding$lengths <= maxgap] <- FALSE
      
      # The original vector is being reconstructed by reverse.rls, only now the
      # longer runs are replaced now in the logical vector derived from is.na()
      # in the beginning all former NAs that are > maxgap are also FALSE
      en <- inverse.rle(rlencoding)
      
      # Set all positions in the imputed series with gaps > maxgap to NA
      # (info from en vector)
      data[en == TRUE] <- NA
    }
    
    ##
    ## End Post Processing
    ##
    
    
    ##
    ## 4. Final Output Formatting
    ##
    
    # Give back the object originally supplied to the function
    # (necessary for multivariate input with only 1 column)
    if (!is.null(dim(x)[2])) {
      x[, 1] <- data
      return(x)
    }
    
    ##
    ## End Final Output Formatting
    ##
    
    return(data)
  }
}
# [na_interpolation] directly copied from R package "imputeTS" directly, avoiding downloading the whole R package "imputeTS"
na_interpolation <- function(x, option = "linear", maxgap = Inf) {
  data <- x
  
  #----------------------------------------------------------
  # Mulivariate Input
  # The next 20 lines are just for checking and handling multivariate input.
  #----------------------------------------------------------
  
  # Check if the input is multivariate
  if (!is.null(dim(data)[2]) && dim(data)[2] > 1) {
    # Go through columns and impute them by calling this function with univariate input
    for (i in 1:dim(data)[2]) {
      if (!anyNA(data[, i])) {
        next
      }
      # if imputing a column does not work - mostly because it is not numeric - the column is left unchanged
      tryCatch(data[, i] <- na_interpolation(data[, i], option, maxgap), error = function(cond) {
        warning(paste("imputeTS: No imputation performed for column", i, "because of this", cond), call. = FALSE)
      })
    }
    return(data)
  }
  
  
  #----------------------------------------------------------
  # Univariate Input
  # All relveant imputation / pre- postprocessing  code is within this part
  #----------------------------------------------------------
  
  else {
    missindx <- is.na(data)
    
    ##
    ## 1. Input Check and Transformation
    ##
    
    
    # 1.1 Check if NAs are present
    if (!anyNA(data)) {
      return(data)
    }
    
    # 1.2 special handling data types
    if (any(class(data) == "tbl")) {
      data <- as.vector(as.data.frame(data)[, 1])
    }
    
    # 1.3 Check for algorithm specific minimum amount of non-NA values
    if (sum(!missindx) < 2) {
      stop("Input data needs at least 2 non-NA data point for applying na_interpolation")
    }
    
    # 1.4 Checks and corrections for wrong data dimension
    
    # Check if input dimensionality is not as expected
    if (!is.null(dim(data)[2]) && !dim(data)[2] == 1) {
      stop("Wrong input type for parameter x")
    }
    
    # Altering multivariate objects with 1 column (which are essentially
    # univariate) to be dim = NULL
    if (!is.null(dim(data)[2])) {
      data <- data[, 1]
    }
    
    # 1.5 Check if input is numeric
    if (!is.numeric(data)) {
      stop("Input x is not numeric")
    }
    
    ##
    ## End Input Check
    ##
    
    
    ##
    ## 2. Imputation Code
    ##
    
    n <- length(data)
    
    allindx <- 1:n
    indx <- allindx[!missindx]
    
    data_vec <- as.vector(data)
    
    if (option == "linear") {
      interp <- stats::approx(indx, data_vec[indx], 1:n, rule = 2)$y
    }
    else if (option == "spline") {
      interp <- stats::spline(indx, data_vec[indx], n = n)$y
    }
    else {
      stop("Wrong parameter 'option' given. Value must be either 'linear', or 'spline'.")
    }
    
    # Merge interpolated values back into original time series
    data[missindx] <- interp[missindx]
    
    ##
    ## End Imputation Code
    ##
    
    
    ##
    ## 3. Post Processing
    ##
    
    # 3.1 Check for Maxgap option
    
    # If maxgap = Inf then do nothing and when maxgap is lower than 0
    if (is.finite(maxgap) && maxgap >= 0) {
      
      # Get logical vector of the time series via is.na() and then get the
      # run-length encoding of it. The run-length encoding describes how long
      # the runs of FALSE and TRUE are
      rlencoding <- rle(is.na(x))
      
      # Runs smaller than maxgap (which shall still be imputed) are set FALSE
      rlencoding$values[rlencoding$lengths <= maxgap] <- FALSE
      
      # The original vector is being reconstructed by reverse.rls, only now the
      # longer runs are replaced now in the logical vector derived from is.na()
      # in the beginning all former NAs that are > maxgap are also FALSE
      en <- inverse.rle(rlencoding)
      
      # Set all positions in the imputed series with gaps > maxgap to NA
      # (info from en vector)
      data[en == TRUE] <- NA
    }
    
    ##
    ## End Post Processing
    ##
    
    
    ##
    ## 4. Final Output Formatting
    ##
    
    # Give back the object originally supplied to the function
    # (necessary for multivariate input with only 1 column)
    if (!is.null(dim(x)[2])) {
      x[, 1] <- data
      return(x)
    }
    
    ##
    ## End Final Output Formatting
    ##
    
    return(data)
  }
}
# [merge.closepoints] merge points departing no bigger than <band>
merge.closepoints=function(points,band){
  points=unique(sort(points))
  to_be_merged=which(as.numeric(diff(points)<band)==1)
  points[to_be_merged]=points[to_be_merged+1]=floor((points[to_be_merged]+points[to_be_merged+1])/2)
  points=unique(sort(points))
  if(all((diff(points)<band)==F)){
    return(points)
  }else{
    return(merge.closepoints(points,band))
  }
}

#####################################################
##        generate and insert missingness         ###
#####################################################
# generate missing indexes of a given type, length and relavent parameters (previously named as missing_index)
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
#  insert missingness into certain variables, given the missing_index
insert.missingness=function(data,variables,missing_index){
  result=data
  result[missing_index,variables]=NA
  return(result)
}

