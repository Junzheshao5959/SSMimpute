#' Main function, state space model multiple imputation with unanimous changepoints
#' @param data_ss_ori contains all information, and only selected variables in formula_var enters the state space model
#' @param formula_var select variables from \code{data_ss} into the statespace model
#' @param ss_param_temp a list of parameters to set up state-space model
#' \itemize{
#' \item{\code{m0}: }{initial values for states}
#' \item{\code{C0}: }{initial values for variance of states}
#' \item{\code{inits}: }{initial values for the estimating of all NA terms, via maximizing likelihood}
#' \item{\code{AR1_coeffi}: }{variables, whose coefficient is a AR(1) process; if none, then is NULL}
#' \item{\code{rw_coeffi}: }{variables, whose coefficient is a random walk process;if none, then is NULL}
#' \item{\code{w_cp_param}: }{variables, whose coefficients are periodic fixed (may shift to other levels over time, but fixed within periods), it contains a list of parameters for each variable whose coefficient level shifts to different values. Details for variable: \code{variable} : the name of the variable (Required); \code{segments} : how many segments of constant coefficient (Required); \code{changepoints} :the corresponding changepoints for the separated segments, it can either be specified by the user or automatically inferred, \code{fixed_cpts} only exist when \code{changepoints} exists }
#' \item{\code{v_cp_param}: }{information about periodic observational variance V (may decrease or increase over time, but fixed within periods), it contains a list of parameters for each variable whose coefficient level shifts to different values. Details for variable: \code{segments} : how many segments of constant coefficient (Required); \code{changepoints} : the corresponding changepoints for the separated segments, it can either be specified by the user or automatically inferred, \code{fixed_cpts} only exist when \code{changepoints} exists}
#' }
#' @param max_iteration control for the convergence of changepoints, a positive integer
#' @param cpt_learning_param a list of variable for change point learning
#' \itemize{
#' \item{\code{cpt_method}: }{either "mean" or "meanvar"}
#' \item{\code{burnin}: }{a positive number in (0,1)}
#' \item{\code{mergeband}: }{a positive integer}
#' \item{\code{AR1_coeffi}: }{variables, whose coefficient is a AR(1) process; if none, then is NULL}
#' \item{\code{rw_coeffi}: }{variables, whose coefficient is a random walk process;if none, then is NULL}
#' \item{\code{w_cp_param}: }{variables, whose coefficients are periodic fixed (may shift to other levels over time, but fixed within periods)}
#' \item{\code{v_cp_param}: }{information about periodic observational variance V (may decrease or increase over time, but fixed within periods)}
#' }
#'
#' @param initial_imputation_option for the first iteration of imputing missing y, choose StructTS or others, and can't be "ignore"
#' @param estimate_convergence_cri critic value for convergence check, default 0.01
#' @param lik_convergence_cri critic value for convergence check, default 0.01
#' @param stepsize_for_newpart stepsize specified, default 1/3
#' @param max_iteration max iteration, default 100
#' @param cpt_initial_guess_option option for initially learning cpts in preparation period
#' @param dlm_option choose between smooth or filter
#' @param m number of draws for multiple imputation
#' @param seed random seed
#' @param printFlag whether we need to print the Flag plots.
#'
#' @return A list of result
#' @export
#' @importFrom stats rnorm
#' @importFrom utils head
#' @importFrom MASS mvrnorm
#' @importFrom crayon red
#' @importFrom crayon blue
#' @import changepoint dlm tidyverse imputeTS
SSMimpute_unanimous_cpts=function(data_ss_ori,formula_var,ss_param_temp,
                                      initial_imputation_option="StructTS",estimate_convergence_cri=0.01,lik_convergence_cri=0.01,stepsize_for_newpart=1/3,max_iteration=100,
                                      cpt_learning_param=list(cpt_method="mean",burnin=1/10,mergeband=20,convergence_cri=15),
                                      cpt_initial_guess_option="ignore", # only this function needs this. no imputation has complete data, no use for StrutST; likelihood no fit model
                                      dlm_option="smooth",
                                      m=5,seed=1,printFlag=T){
  # !! Please read this paragraph first !!
  # Limitation 1: this current version of function has only tested on models, no more complicated than y ~ y_1 + x + x_1 + c
  #               for a continuous y, a binary/continuous x, and a continuous c,
  #               other scenarios, including (i) interactions and (ii) binary y and c, are not allowed and haven't been tested.
  # Limitation 2: this current version of function only allow for simutaniously learning of the shared changepoints for w_cpt_param.
  #               Thus, it only works for (i) a initial collection of changepoints are given for a collection of variables,
  #                                           those changpoints are fixed, and no more updating
  #                                       (ii) a initial collection of changepoints are given for a collection of variables,
  #                                            but those changepoints are not for sure, and need updating
  #                                       (iii) no initial guess is given in terms of changepoints, but only how many
  #                                             the changepoints learnt need to be the same across all variables, whose coefficient shifts levels
  #                Working flowchart for [changepoints]
  #                (i) if changepoints are given in <changepoints> and <fixed_cpts>=T for all variabels in [w_cp_param]
  #                     -> 1) no initailization
  #                        2) no estimation for the location of changepoints
  #                        3) no iteration until convergence for changepoints
  #                     we allow changepoints differ for each variables, when they are all given in advance in this case
  #                (ii) if changepoints are given in <changepoints> and <fixed_cpts>=F for all variabels in [w_cp_param]
  #                     -> 1) no initailization
  #                        2) estimate the location of changepoints for those variables with fixed_cpts=F
  #                        3) update changepoints via iteration, until convergence achieved for those variables with fixed_cpts=F
  #                     changepoints should be same for variables for with fixed_cpts=F
  #                (iii) if changepoints are not given and thus <fixed_cpts> do not exist in [w_cp_param]
  #                    -> 1) initialization
  #                       2) require estimation for the location of changepoints
  #                       3) require iteration until convergence for changepoints
  #                     we require changepoints should be same across all variables with shifed levels
  #                Otherwise, all other cases are not permitted for the moment, including
  #                    -> [not allowed] 1) changepoints are partially given
  #                    -> [not allowed] 2) <fixed_cpts> are partially T
  # Limitation 3: currently the initiation of changepoints parts miss the update of data with interaction between y_1 and x
  # Limitation 4[to be changed]: simple and no iterative "StructTS" imputation for learning cpts in V
  #                              estimated_cpts and iter in the return result only reflects those for W, not for V

  # Explanation for parameters:
  # <data_ss_ori>: contains all information, and only selected variables in formula_var enters the statespace model
  # <formula_var>: select variables from <data_ss_ori> into the statespace model
  # <ss_param_temp>:
  #             <m0>: initial values for states
  #             <C0>: initial values for variance of states
  #             <inits>: initial values for the estimating of all NA terms, via maximizing likelihood
  #             <AR1_coeffi>: variables, whose coefficient is a AR(1) process;
  #                           if none, then is NULL
  #             <rw_coeffi>: variables, whose coefficient is a random walk process;
  #                          if none, then is NULL
  #             <w_cp_param>: variables, whose coefficients are periodic fixed (may shift to other levels over time, but fixed within periods)
  #                          [structure] a list of lists, containing <variable>, <segments>, <changepoints>, <fixed_cpts> for each variable whose coefficient level shifts to different values
  #                                      - <variable> the name of the variable [must exist]
  #                                      - <segments> how many segments of constant coefficient [must exist]
  #                                      - <changepoints> the corresponding changepoints for the separated segments
  #                                        [note]: <changepoints> can be learnt from <segments>, or directly given
  #                                      - <fixed_cpts>: only exist when <changepoints> exists
  #                          if none, then is NULL
  #                          [Requirement] <changepoints> + <fixed_cpts> exist either for all variables or for none  -> [may be futher altered]
  #                          [Requirement] when <changepoints> is estimated, all changepoints should be the same across different variables  -> [may be futher altered]
  #                          [Requirement] when <changepoints> is given, changepoints for different variables may differ
  #             <v_cp_param>: information about periodic observational variance V (may decrease or increase over time, but fixed within periods)
  #                          [structure] only one list containing <segments>, <changepoints>, and <fixed_cpts>
  #                                      - <segments>: how many segments of constant coefficient [must exist]
  #                                      - <changepoints>: the corresponding changepoints for the separated segments
  #                                        [note]: <changepoints> can be learnt from <segments>, or directly given
  #                                      - <fixed_cpts>=T: only exist when <changepoints> exists
  #                          if none, then is NULL
  # <initial_imputation_option>: for the first iteration of imputing missing y, choose StructTS or others, and can't be "ignore"
  # <lik_convergence_cri>:
  # <max_iteration>: control for the convergence of EM algorithm
  # <cpt_learning_param>:  <cpt_method> either "mean" or "meanvar"
  #                        <burnin> a positive number in (0,1)
  #                        <mergeband> a positive integer
  #                        <convergence_cri> a positive integer
  #                        Caution: all used for learning changepoints in coefficients, not for changepoints in observational variance
  # <cpt_initial_guess_option>: option for initially learning cpts in preparation period
  # <dlm_option>: choose between smooth or filter
  # <m>: number of draws for multiple imputation

  # Example 1: (fixed effects for all coefficients)
  #            ss_param=list(inits=c(log(1)),m0=c(40,0.5,-1.5,-0.5,1),C0=diag(rep(10^3),5),
  #                          AR1_coeffi=NULL,rw_coeffi=NULL,v_cp_param=NULL,w_cp_param=NULL,max_iteration=50)
  # Example 2: (AR(1) for C's coefficient and fixed for all other coefficients)
  #            ss_param=list(inits=c(0,-4.60517,0), # 0 for autocorrelation, log(0.01) for W_i, last 0 for V
  #                          m0=c(40,0.5,-1.5,-0.5,1,0.5),C0=diag(rep(10^3),6),
  #                          AR1_coeffi="c",rw_coeffi=NULL,v_cp_param=NULL,w_cp_param=NULL,max_iteration=50)
  # Example 3: (random walk for C's coefficient and fixed for all other coefficients)
  #            ss_param=list(inits=c(0,-4.60517),m0=c(40,0.5,-0.5,-1.5,1.5),C0=diag(rep(10^3),5),
  #                          AR1_coeffi=NULL,rw_coeffi="c",v_cp_param=NULL,w_cp_param=NULL,max_iteration=50)
  # Example 4: (AR(1) for both X and X_1, and fixed for all other coefficients)
  #            ss_param=list(inits=c(0,0,-4.60517,-4.60517,0), # first and second 0 for 2 autocorrelations, 3rd and 4th for W_i of X and X_1, and last 0 for V
  #                          m0=c(40,0.5,-1.5,-0.5,1,-0.75,-0.25),C0=diag(rep(10^3),7),
  #                          AR1_coeffi=c("x","x_1"),rw_coeffi=NULL,v_cp_param=NULL,w_cp_param=NULL,max_iteration=50)
  # Example 5: (periodic V of 3 periods)
  #            ss_param=list(inits=c(log(1),log(20),log(1)),m0=c(40,0.5,-1.5,-0.5,1),C0=diag(rep(10^3),5),
  #                          AR1_coeffi=NULL,rw_coeffi=NULL,
  #                          v_cp_param=list(segments=3),w_cp_param=NULL,max_iteration=50)
  # Example 5.1: (periodic V of 3 known periods)
  #              ss_param=list(inits=log(c(1,20,1)),m0=c(40,0.5,-1.5,-0.5,1),C0=diag(rep(10^3),5),
  #                            AR1_coeffi=NULL,rw_coeffi=NULL,
  #                            v_cp_param=list(segments=3,changepoints=c(400,700),fixed_cpts=T),w_cp_param=NULL,max_iteration=max_iteration)
  # Example 6: (periodic coefficient for X of 3 periods)
  #            ss_param=list(inits=c(log(1)),m0=c(40,0.5,-1.5,-0.5,1),C0=diag(rep(10^3),5),
  #                          AR1_coeffi=NULL,rw_coeffi=NULL,v_cp_param=NULL,
  #                          w_cp_param=list(list(variable="x",segments=3)),max_iteration=50)
  # Example 6.1: (periodic coefficient for X of 3 known periods)
  #              ss_param=list(inits=c(log(1)),m0=c(40,0.5,-1.5,-0.5,1),C0=diag(rep(10^3),5),
  #                       AR1_coeffi=NULL,rw_coeffi=NULL,
  #                       v_cp_param=NULL,
  #                       w_cp_param=list(list(variable="x","segments"=3,changepoints=c(400,700),fixed_cpts=T)),
  #                       max_iteration=max_iteration)
  # Example 6.2: (periodic coefficient for X and X_1 of 3 uncertain known periods)
  #              ss_param=list(inits=c(log(1)),m0=c(40,0.5,-1.5,-0.5,1),C0=diag(rep(10^3),5),
  #                            AR1_coeffi=NULL,rw_coeffi=NULL,
  #                            v_cp_param=NULL,
  #                            w_cp_param=list(list(variable="x","segments"=3,changepoints=c(400,700),fixed_cpts=F),
  #                                            list(variable="x_1","segments"=3,changepoints=c(400,700),fixed_cpts=F)),
  #                            max_iteration=max_iteration)
  # Example 7: (period v of 3 periods and periodic coefficient for X of 3 periods)
  #            ss_param=list(inits=c(log(1),log(20),log(1)),m0=c(40,0.5,-1.5,-0.5,1),C0=diag(rep(10^3),5),
  #                          AR1_coeffi=NULL,rw_coeffi=NULL,
  #                          v_cp_param=list(segments=3),
  #                          w_cp_param=list(list(variable="x",segments=3)),max_iteration=50)

  # ------------------- Preparation ---------------------- #
  # check: 1. all variables specified in the model exist the data_ss_ori set
  #        2. ss_param is a list
  #           1) m0 is either NULL(use the default in dlm), or given (contains no NA, all numeric, same length as states)
  #           2) C0 is either NULL(use the default in dlm), or given (contains no NA, all numeric, same dim as n_states*n_states)
  #           3) the variables, whose coefficiens is AR(1) process, are char string, no NA, and in both model and data_ss_ori set
  #           4) the variables, whose coefficiens is randeom walk process, are char string, no NA, and in both model and data_ss_ori set
  #           5) check if changepoints of V need to be estimated
  #              -> as no imputation is applied for y -> no change on data_ss_ori -> no iteration needed
  #              In short, changepoints of V is only estimated once, in the initialization part
  #              Otherwise, if changepoints are given, no action is needed.
  #           6) check cpt_learning_param under w_param
  #           7) check if changepoints of W for certain variables need to be estimated
  #           8) check if the length of inits is correct
  #         3. initial_imputation_option: only "StructTS"
  #         4. lik_convergence_cri: a positive number
  #         5. max_iteration is positive
  #         6. cpt_learning_param is a list
  #            1) <cpt_method>: "mean" or "meanvar"
  #            2) <burnin>: a positive number in (0,1)
  #            3) <mergeband>: a positive integer
  #            4) <convergence_cri>: a positive integer
  #         7. <cpt_initial_guess_option>: "StructTS" or "ignore"
  #         8. dlm_option: smooth or filter
  #         9. m: a positive integer

  if(!all(formula_var %in% colnames(data_ss_ori))){
    stop("The variables specified in the formula are not contained in the data_ss_ori set.")
  }
  if(!is.list(ss_param_temp)){
    stop("ss_param must be a list.")
  }else{
    m0=ss_param_temp$m0                    # initial for states -> checked
    C0=ss_param_temp$C0                    # initial for variance of states -> checked
    inits=ss_param_temp$inits              # initial for NA terms in mle
    AR1_coeffi=ss_param_temp$AR1_coeffi    # add 1 additional state and 2 additional NA (autocorrelation and W) -> checked
    rw_coeffi=ss_param_temp$rw_coeffi      # add 1 additional NA (W) -> checked
    w_cp_param=ss_param_temp$w_cp_param    # change points to be estimated -> change W_t to be 10
    v_cp_param=ss_param_temp$v_cp_param    # change points estimated once in the begining, and add #(segments-1) NA for V -> checked

    # n_states includes: 1 for intercept,
    #                    n_variable for all selected variable,
    #                    length(AR1_coeffi) for additional baselines for the AR(1) process
    #              Note: random walk only change W, and do not add states
    #                    shift in beta only change W to 10, do not add anything
    #                    varying V only change segments of V_t, do not add states
    n_states=1+length(formula_var)+length(AR1_coeffi)
    if(!is.null(m0)){ # m0 is either NULL(use the default in dlm), or given (contains no NA, all numeric, same length as states)
      if(!is.numeric(m0)){stop("Error in m0: m0 is not all numeric.")}
      if(!all(!is.na(m0))){stop("Error in m0: m0 contains NA.")}
      if(length(m0)!=n_states){stop("Error in m0: m0 has different length as that of the states.")}
    }
    if(!is.null(C0)){ # C0 is either NULL(use the default in dlm), or given (contains no NA, all numeric, same dim as n_states*n_states)
      if(!is.numeric(C0)){stop("Error in C0: C0 is not numeric.")}
      if(!all(!is.na(C0))){stop("Error in C0: C0 contains NA.")}
      if(!all(dim(C0)==n_states)){stop("Error in C0: C0 has different dimention than n_states*n_states.")}
    }
    if(!is.null(AR1_coeffi)){ # variables are of type 'char' string, no NA, and in both model and data_ss set
      if(!is.character(AR1_coeffi)){stop("Error in AR1_coeffi: the specified variables in AR1_coeffi are not all char string.")}
      if(!all(AR1_coeffi %in% formula_var)){stop("Error in AR1_coeffi: the selected time-varying coefficient has no corresponding variable in the regression formula.")}
      if(!all(AR1_coeffi %in% colnames(data_ss_ori))){stop("Error in AR1_coeffi: the selected time-varying coefficient has no corresponding variable in the new data_ss_ori set.")}
      if(!all(!is.na(AR1_coeffi))){stop("Error in AR1_coeffi: AR1_coeffi contains NA.")}
    }
    if(!is.null(rw_coeffi)){ # variables are char string, no NA, and in both model and data_ss set
      if(!is.character(rw_coeffi)){stop("Error in rw_coeffi: the specified variables in rw_coeffi are not all char string")}
      if(!all(rw_coeffi %in% c("intercept",formula_var))){stop("Error in rw_coeffi: the selected time-varying coefficient has no corresponding variable in the regression formula.")}
      if(!all(rw_coeffi %in% c("intercept",colnames(data_ss_ori)))){stop("Error in rw_coeffi: the selected time-varying coefficient has no corresponding variable in the data_ss_ori set.")}
      if(!all(!is.na(rw_coeffi))){stop("Error in rw_coeffi: rw_coeffi contains NA.")}
    }
    if(!is.null(v_cp_param)){
      v_cp_param_learncps=v_cp_param
      # check 1) only allow [segments], [changepoints], or [fixed_cpts]
      #       2) #elements is 1-3.
      #       3-4) [segments] must exist, and must be a positive number
      #       5) when either [changepoints] or [fixed_cpts] are detected: (for the case when cpts are given)
      #           5.1) both [changepoints] and [fixed_cpts] exist.
      #           5.2) changepoints are all positive numbers, and within the range of timeline
      #          when [changepoints] doesn't exist, (for the case when cpts aren't given)
      #           as NA is not allowed in the program of finding change points,
      #           if there are NA in the original outcome -> do one time imputation "StructTS"
      #           -> no iteration is needed, without iterative imputation
      #       6) [segments]-1 must equals to #changepoints
      if(!all(names(v_cp_param_learncps) %in% c("segments","changepoints","fixed_cpts"))){stop("Error in v_cp_param: only [segments], [changepoints], [fixed_cpts] are allowed.")}
      if(!(length(v_cp_param_learncps)>0 & length(v_cp_param_learncps)<4)){stop("Error in v_cp_param: length is [1,3].")}
      if(!("segments" %in% names(v_cp_param_learncps))){stop("Error in v_cp_param: <segments> must exist for v_cp_param.")}
      if(!(is.numeric(v_cp_param_learncps$segments) &
           length(v_cp_param_learncps$segments)==1 &
           v_cp_param_learncps$segments>0)){stop("Error in v_cp_param: <segments> is not a positive number.")}
      if(!all(!names(v_cp_param_learncps) %in% c("changepoints","fixed_cpts"))){ # if detected one of changepoints or fixed_cpts
        if(!all(c("changepoints","fixed_cpts") %in% names(v_cp_param_learncps))){stop("Error in v_cp_param: [changepoints] and [fixed_cpts] must exist at the same time.")}
        if(!is.numeric(v_cp_param_learncps$changepoints)){stop("Error in v_cp_param: [changepoints] is not all numeric")}
        if(!all(v_cp_param_learncps$changepoints > 1 & v_cp_param_learncps$changepoints < nrow(data_ss_ori))){stop("Error in v_cp_param: [changepoints] is not range of [1,nrow(data_ss_ori)-1].")}
        if(printFlag){
          if(v_cp_param_learncps$fixed_cpts){
            cat(red("Note: ''changepoints'' of V is given and fixed.\n"))
          }else{
            cat(red("Note: ''changepoints'' of V is given a initial guess, but need to be updated according to the fitting result.\n"))
          }
        }
      }else{
        # when initial guess of changepoints of V is unkwown
        ######### initialize changepoints for V (part I of initialization)  #########
        if(sum(is.na(data_ss_ori$y))!=0){
          # for the ignore case
          cpt_v_temp=cpt.var(na_kalman(data_ss_ori$y,model="StructTS"),penalty="Manual",Q=v_cp_param_learncps$segments-1, method="BinSeg")
        }else{
          # for the full and cc case
          cpt_v_temp=cpt.var(data_ss_ori$y,penalty="Manual",Q=v_cp_param_learncps$segments-1, method="BinSeg")
        }
        v_cp_param[["changepoints"]]=cpts(cpt_v_temp)
        v_cp_param[["fixed_cpts"]]=T
        if(printFlag){
          invisible(plot(cpt_v_temp))
          cat("Estimated ''changepoints'' for V are",cpts(cpt_v_temp),".\n")
        }
      }
      if(length(v_cp_param[["changepoints"]]) != (v_cp_param[["segments"]]-1)){
        stop("For varying V, the number <changepoints> doesn't agree with those indicated by <segments>")
      }
    }
    if(!is.null(w_cp_param)){
      # only when w_cp_param (shifted levels of coefficients) exists and cpts for W need to be estimated (either not given, or given but not fixed)
      #  ->  we need [cpt_learning_param]
      if(!is.null(cpt_learning_param) & is.list(cpt_learning_param)){
        if(all(c("cpt_method","burnin","mergeband","convergence_cri") %in% names(cpt_learning_param)) & length(unique(names(cpt_learning_param)))==4){
          if(length(cpt_learning_param$cpt_method) == 1 & is.character(cpt_learning_param$cpt_method) & !all(cpt_learning_param$cpt_method %in% c("mean","meanvar"))){stop("Error in cpt_learning_param$cpt_method: should be mean or meanvar.")}
          if(!(length(cpt_learning_param$burnin)==1 & is.numeric(cpt_learning_param$burnin) & cpt_learning_param$burnin>0 & cpt_learning_param$burnin<1)){stop("Error in cpt_learning_param$burnin: should be number in (0,1).")}
          if(!(length(cpt_learning_param$mergeband)==1 & is.numeric(cpt_learning_param$mergeband) & cpt_learning_param$mergeband>0)){stop("Error in cpt_learning_param$mergeband: should be a positive number.")}
          if(!(length(cpt_learning_param$convergence_cri)==1 & is.numeric(cpt_learning_param$convergence_cri) & cpt_learning_param$convergence_cri>0)){stop("Error in cpt_learning_param$mergeband: should be a positive number.")}
          cpt_method=cpt_learning_param$cpt_method
          burnin=cpt_learning_param$burnin
          mergeband=cpt_learning_param$mergeband
          cpts_convergence_cri=cpt_learning_param$convergence_cri
        }else{
          stop("Error in cpt_learning_param: this is a list of 4 elements of [cpt_method],[burnin],[mergeband], and [convergence_cri].")
        }
      }else{
        stop("Error in cpt_learning_param: cpt_learning_param (list) is needed when some variable has shifted levles of coefficient.")
      }
      # only when w_cp_param (shifted levels of coefficients), we need [cpt_learning_param]
      if(!is.null(cpt_initial_guess_option)){
        if(!(length(cpt_initial_guess_option)==1 & is.character(cpt_initial_guess_option) & cpt_initial_guess_option %in% c("StructTS","ignore"))){
          stop("Error in cpt_initial_guess_option: only [StructTS] and [ignore] are allowed!")
        }
      }else{
        stop("Error in cpt_initial_guess_option: missing!")
      }
    }
    if(!is.null(w_cp_param)){
      w_cp_param_learncps=w_cp_param
      # Check: 1) The variables with changepoints must be in data_ss_ori and in formula
      #           the order of variables aligns with the order in the formula
      #        2) for each variable with changepoints,
      #           check 2.1) [variable], [segments], [changepoints], [fixed_cpts] may exist
      #                 2.2) #elements for each sublist is 1-4.
      #                 2.3) [variable] and [segments] must exist
      #                      [variable] must be a char
      #                      [segments] must be a positive number
      #                      [segments] for all variables must be the same
      #                 2.4) when [changepoints] exists -> [fixed_cpts] must exist and be logical
      #                 2.5) when [changepoints] doesn't exist, it need to be estimated only once from the data_ss_ori
      #                      -> no iteration is needed, without iterative imputation
      #                 2.6) [variable], [segments], [changepoints] must exist for all variables after possibly estimation
      #            Note: # changepoints need to strictly equal to [segments]-1 eventually, but not necessary in the intialization
      # check 2.1)-2.3)
      for(aa in 1:length(w_cp_param_learncps)){
        if(!all(names(w_cp_param_learncps[[aa]]) %in% c("variable","segments","changepoints","fixed_cpts"))){stop("Error in w_cp_param: only [variable], [segments], [changepoints], [fixed_cpts] for each variable are allowed.")}
        if(!(length(w_cp_param_learncps[[aa]])>0 & length(w_cp_param_learncps[[aa]])<5)){stop("Error in w_cp_param: the length of each list is wrong.")}
        if(!all(c("variable","segments") %in% names(w_cp_param_learncps[[aa]]))){stop("Error in w_cp_param: <variable> and <segments> must exist for each level")}
        if(!(is.character(w_cp_param_learncps[[aa]]$variable) &
             length(w_cp_param_learncps[[aa]]$variable)==1)){stop("Error in w_cp_param: <variable> is not a char.")}
        if(!(is.numeric(w_cp_param_learncps[[aa]]$segments) &
             length(w_cp_param_learncps[[aa]]$segments)==1 &
             w_cp_param_learncps[[aa]]$segments>0)){stop("Error in w_cp_param: <segments> is not a positive integer.")}
      }
      if(length(unique(unlist(lapply(1:length(w_cp_param_learncps),function(aa1){w_cp_param_learncps[[aa1]][["segments"]]}))))>1){
        stop("Error in w_cp_param: the segments value for different variables need to be the same")
      }
      # check 1)
      w_cp_param_variables=unlist(lapply(1:length(w_cp_param_learncps),function(bb){w_cp_param_learncps[[bb]][["variable"]]}))
      if(!all(w_cp_param_variables %in% formula_var)){stop("The selected level-shifted coefficient has no corresponding variable in the regression formula.")}
      if(!all(w_cp_param_variables %in% colnames(data_ss_ori))){stop("The selected level-shifted coefficient has no corresponding variable in the data_ss_ori set.")}
      if(length(w_cp_param_variables)>1){
        w_cp_param_variables_order=c()
        for(bb in w_cp_param_variables){
          w_cp_param_variables_order=c(w_cp_param_variables_order,which(bb == c("intercept",formula_var)))
        }
        if(!all(diff(w_cp_param_variables_order)>0)){stop("Error in w_cp_param: the variable order should be the same as in formula.")}
      }
      # check 2.4)-2.5):  if <changepoints> exist and "fixed_cpts"=T for all variables -> no need to learn changespoints and no need for iteration
      #                   if <changepoints> exist and "fixed_cpts"=F for all variables -> no need to learn changespoints and need further iteration
      #                   if <changepoints> do not exist -> "fixed_cpts" shouldn't exist either
      #                                      -> generate changepoints and iterate
      #                   Caution: For now, we require changepoints+fixed_cpts to exist for all variables or for none
      #                            when changepoints+fixed_cpts exist, fixed_cpts =T for all or =F for all
      #                   Caution 2: For initialization, we do not require # segments=1 to be equal to the length of the changepoints
      if(all(unlist(lapply(1:length(w_cp_param_learncps), function(cc1){ "changepoints" %in% names(w_cp_param_learncps[[cc1]]) & "fixed_cpts" %in% names(w_cp_param_learncps[[cc1]])})))){
        # if [changepoints] and [fixed_cpts] both exist for each variable
        if(!all(unlist(lapply(1:length(w_cp_param_learncps),function(cc3){is.numeric(w_cp_param_learncps[[cc3]]$changepoints) &
            all(w_cp_param_learncps[[cc3]]$changepoints>1 & w_cp_param_learncps[[cc3]]$changepoints < nrow(data_ss_ori))})))){
          stop("Error in w_cp_param: [changepoints] for all variables are not all numeric and within the right range.")
        }
        if(length(unique(unlist(lapply(1:length(w_cp_param_learncps),function(aa2){w_cp_param_learncps[[aa2]][["changepoints"]]}))))!=w_cp_param_learncps[[1]]$segments-1){
          stop("Error in w_cp_param: given [changepoints] may not be same for all variables")
        }
        if(all(unlist(lapply(1:length(w_cp_param_learncps), function(cc4){w_cp_param_learncps[[cc4]]$fixed_cpts})))){
          cat(red("Note: ''changepoints'' of W is given and fixed.\n"))
        }else if(all(unlist(lapply(1:length(w_cp_param_learncps), function(cc4){!w_cp_param_learncps[[cc4]]$fixed_cpts})))){
          cat(red("Note: ''changepoints'' of W is given a initial guess, but need to be updated according to the fitting result.\n"))
        }else{
          stop("Error in w_cp_param: [fixed_cpts] need to be all T or all F when changepoints exist.")
        }
      }else if(all(unlist(lapply(1:length(w_cp_param_learncps), function(cc2){ (!"changepoints" %in% names(w_cp_param_learncps[[cc2]])) & (!"fixed_cpts" %in% names(w_cp_param_learncps[[cc2]]))})))){
        # if [changepoints] and [fixed_cpts] both do not exist for each varible
        if(printFlag){cat("The changepoints of variables, whose coefficients shift over time, are not given.\n")}
        ######### initialization of changepoints of W (part II of initialization) ##########
        # if the changepoints don't exist, we ask the corresponding coeffi be random walk, and then learn changepoints from the results
        ss_param_learncps=ss_param_temp
        # 1) remove the w_cp_param part
        ss_param_learncps$w_cp_param=NULL
        # 2) adjust the random walk part
        # w_cp_param_variables=unlist(lapply(1:length(w_cp_param_learncps),function(ee){w_cp_param_learncps[[ee]][["variable"]]}))
        ss_param_learncps$rw_coeffi=unique(c(ss_param_learncps$rw_coeffi,w_cp_param_variables))
        # 3) generate new inits
        v_part_length=ifelse(is.null(ss_param_learncps$v_cp_param),1,ss_param_learncps$v_cp_param$segments)
        v_part=ss_param_learncps$inits[(length(ss_param_learncps$inits)-v_part_length+1):length(ss_param_learncps$inits)]
        ar1_part_length=length(ss_param_learncps$AR1_coeffi)
        if(ar1_part_length>0){
          ar1_part=ss_param_learncps$inits[1:(2*ar1_part_length)]
        }else{
          ar1_part=NULL
        }
        ss_param_learncps$inits=c(ar1_part, rep(0,length(ss_param_learncps$rw_coeffi)),v_part)
        # [Unique point] (begin)
        # create a initial copy of data_ss_ori as data_ss_temp_init_cpts_w
        data_ss_temp_init_cpts_w=data_ss_ori
        if(cpt_initial_guess_option=="StructTS"){ # update data method 1 (impute using really advance univariate time series imputation method)
          missing_part=which(is.na(data_ss_temp_init_cpts_w$y))[which(is.na(data_ss_temp_init_cpts_w$y))<nrow(data_ss_temp_init_cpts_w)] # if the last y is missing, ignore it
          data_ss_temp_init_cpts_w$y_1[missing_part+1]=na_kalman(data_ss_temp_init_cpts_w$y,model="StructTS")[missing_part]
        }else if(cpt_initial_guess_option=="ignore"){ # update data method 2 (ignore missingness)
          ignore_part=which(rowSums(is.na(data_ss_temp_init_cpts_w))!=0)
          for(l in formula_var){
            data_ss_temp_init_cpts_w[is.na(data_ss_temp_init_cpts_w[,l]),l]=mean(data_ss_temp_init_cpts_w[,l],na.rm=T)
          }
          data_ss_temp_init_cpts_w$y[ignore_part]=NA
        }
        # fit dlm model
        out_filter_init_cpts_w=run.SSM_unanimous_cpts(data_ss=data_ss_temp_init_cpts_w,formula_var=formula_var,ss_param_temp=ss_param_learncps,
                                                      max_iteration=max_iteration,cpt_learning_param=cpt_learning_param,
                                                      dlm_option=dlm_option,printFlag=printFlag)$out_filter
        out_smooth_init_cpts_w=dlmSmooth(out_filter_init_cpts_w)
        # [Unique point] (end)
        # get the changepoints -> and update it into the real ss_param
        w_cps_new=c()
        start_pt=floor(nrow(out_filter_init_cpts_w$m)*burnin)
        for(ee in 1:length(w_cp_param)){
          if(dlm_option=="filter"){
            temp=out_filter_init_cpts_w$m[start_pt:nrow(out_filter_init_cpts_w$m),which(colnames(out_filter_init_cpts_w$mod$GG)==w_cp_param[[ee]]$variable)]
          }else if(dlm_option=="smooth"){
            temp=out_smooth_init_cpts_w$s[start_pt:nrow(out_smooth_init_cpts_w$s),which(colnames(out_filter_init_cpts_w$mod$GG)==w_cp_param[[ee]]$variable)]
          }
          if("meanvar" %in% cpt_method){
            cpt_temp=cpt.meanvar(temp,penalty="Manual",Q=ss_param_temp$w_cp_param[[ee]]$segments-1,method="BinSeg") # changepoints learning method cannot be chosen
          }
          if("mean" %in% cpt_method){
            cpt_temp=cpt.mean(temp,penalty="Manual",Q=ss_param_temp$w_cp_param[[ee]]$segments-1,method="BinSeg") # changepoints learning method cannot be chosen
          }
          if(printFlag){plot(cpt_temp)}
          w_cps_new=c(w_cps_new,cpts(cpt_temp)+start_pt-1)
        }
        w_cps_new=merge_closepoints(points=w_cps_new,band=mergeband)
        if(printFlag){cat(red("The inital guess of changepoints are:",w_cps_new,"\n"))}
        for(ff in 1:length(w_cp_param)){
          w_cp_param[[ff]][["changepoints"]]=w_cps_new
          w_cp_param[[ff]][["fixed_cpts"]]=F
        }
      }else{
        stop("Error in w_cp_param: [changepoints] and [fixed_cpts] must both exist or not exist for all variables")
      }
      # final check 2.6): <variable>, <segments> and <changpoints> all exist for each variable
      for(gg in 1:length(w_cp_param)){
        if(!all(c("variable","segments","changepoints","fixed_cpts") %in% names(w_cp_param[[gg]]))){
          stop("<variable>, <segments>, <changepoints>, and <fixed_cpts> must exist for each level shift coefficient.")
        }
      }
    }
    if(is.null(inits)){
      stop("inits is required before estimating NA in dlm model.")
    }else{
      if(!all(!is.na(inits))){stop("Error in inits: NA in inits.")}
      if(!is.numeric(inits)){stop("Error in inits: inits is not numeric")}
      if(length(inits)!=(length(AR1_coeffi)*2+length(rw_coeffi)+ifelse(is.null(v_cp_param),1,v_cp_param$segments))){stop("Error in inits: Error in the length of inits.")}
    }
  }
  if(!(length(initial_imputation_option)==1 & is.character(initial_imputation_option) & initial_imputation_option %in% c("StructTS"))){stop("Error in initial_imputation_option: choose between ignore or StructTS.")}
  if(!(length(lik_convergence_cri) & is.numeric(lik_convergence_cri) & lik_convergence_cri>0)){stop("Error in lik_convergence_cri.")}
  if(!(length(max_iteration)==1 & is.numeric(max_iteration) & max_iteration>0 & max_iteration%%1==0)){stop("max_iteration is not a positive integer!")}
  if(!(length(dlm_option)==1 & is.character(dlm_option) & dlm_option %in% c("smooth","filter"))){stop("Error in dlm_option: choose between smooth or filter.")}
  if(!(length(m)==1 & is.numeric(m) & m>1 & m%%1==0)){stop("m in MC is not a positive integer!")}

  # functions to build up statespace model, which certains terms to be estimates as NA
  # the items in u are: [1:length(AR1_coeffi)] terms -> logit(autocorrelation) for all AR(1) coefficients
  #                                                     default value is 0=exp(0.5)/(1+exp(0.5))
  #                     [length(AR1_coeffi)+1 : 2*length(AR1_coeffi)] -> log(variance for Wj) for all AR(1) coefficients
  #                     [2*length(AR1_coeffi)+1 : 2*length(AR1_coeffi)+length(rw_coeffi)] -> log(variance for Wj) for random walk coefficients
  #                     [2*length(AR1_coeffi)+length(rw_coeffi)+1:2*length(AR1_coeffi)+length(rw_coeffi)+length(v_changepoint)] -> (possibly sereis of) log(variance for V)
  buildCamp_iterative=function(u){
    model=dlmModReg(data_temp[,formula_var],dV=exp(u[2*length(AR1_coeffi)+length(rw_coeffi)+1]))
    FF=matrix(c(model$FF,rep(0,length(AR1_coeffi))),nrow=1)
    if(length(AR1_coeffi)>0){
      states_name=c("intercept",formula_var,paste("rho.",AR1_coeffi,sep=""))
    }else{
      states_name=c("intercept",formula_var)
    }
    colnames(FF)=states_name
    GG=diag(rep(1,length(states_name)))
    colnames(GG)=rownames(GG)=states_name
    W=diag(rep(0,length(states_name)))
    colnames(W)=rownames(W)=states_name
    if(length(AR1_coeffi)>0){
      for(l in 1:length(AR1_coeffi)){
        variable=AR1_coeffi[l]
        GG[variable,grep(paste("^rho.",variable,"$",sep=""),colnames(GG))]=1
        GG[variable,variable]=exp(u[l])/(1+exp(u[l]))
        W[variable,variable]=exp(u[l+length(AR1_coeffi)])
      }
    }
    if(length(rw_coeffi)>0){
      for(k in 1:length(rw_coeffi)){
        variable=rw_coeffi[k]
        W[variable,variable]=exp(u[2*length(AR1_coeffi)+k])
      }
    }
    JFF=matrix(c(model$JFF,rep(0,ncol(FF)-ncol(model$JFF))),nrow=1)
    model$FF=FF
    model$JFF=JFF
    model$GG=GG
    model$W=W
    model$C0=diag(rep(1e07,ncol(model$FF)))
    model$m0=rep(0,ncol(model$FF))
    if(!is.null(m0) & all(!is.na(m0))){
      model$m0=m0
    }
    if(!is.null(C0) & all(!is.na(C0))){
      model$C0=C0
    }
    # adjust model for time-varying V
    if(!is.null(v_cp_param)){
      JV(model)=ncol(model$X)+1
      v_changepoint=c(1,v_cp_param$changepoints,nrow(model$X))
      X_JV=rep(NA,nrow(model$X))
      for(jv in 1:(length(sort(v_changepoint))-1)){
        X_JV[v_changepoint[jv]:v_changepoint[jv+1]]=exp(u[2*length(AR1_coeffi)+length(rw_coeffi)+jv])
      }
      model$X=cbind(model$X,X_JV)
    }
    # adjust model for sudden shift in coefficients
    if(!is.null(w_cp_param)){
      JW=model$W
      JW[JW!=0]=0
      for(jw in 1:length(w_cp_param_variables)){
        JW[w_cp_param_variables[jw],w_cp_param_variables[jw]]=ncol(model$X)+jw
      }
      w_changepoint=matrix(0,nrow=nrow(model$X),ncol=length(w_cp_param_variables))
      colnames(w_changepoint)=w_cp_param_variables
      for(jw in 1:length(w_cp_param_variables)){
        jw_cps_temp=unique(unlist(lapply(w_cp_param[[jw]]$changepoints,function(jwpoint){(jwpoint-10):(jwpoint+10)})))
        jw_cps_temp=jw_cps_temp[jw_cps_temp>0 & jw_cps_temp<=nrow(model$X)]
        w_changepoint[jw_cps_temp,w_cp_param_variables[jw]]=10
      }
      X(model)=cbind(model$X,w_changepoint)
      JW(model)=JW
    }
    dlm(model)
    return(model)
  }
  get.filter_result=function(out_filter,formula_var,AR1_coeffi,rw_coeffi,w_cp_param){
    # organize result
    result_var_names=c("(Intercept)",formula_var)

    if(is.null(w_cp_param)){
      Estimate=Std.Error=rep(NA,length(formula_var)+1)
      w_cp_param_variables=NULL
      result_var_names_new=result_var_names
    }else{
      Estimate=Std.Error=rep(NA,length(formula_var)+1+length(unlist(lapply(1:length(w_cp_param),function(i){w_cp_param[[i]][["changepoints"]]}))))
      result_var_names_new=result_var_names
      for(rr in 1:length(w_cp_param)){
        pointer_id=grep(paste("^",w_cp_param[[rr]]$variable,"$",sep=""),result_var_names)
        temp_id=grep(paste("^",w_cp_param[[rr]]$variable,"$",sep=""),result_var_names_new)
        result_var_names_new=c(result_var_names_new[1:(temp_id-1)],
                               paste(w_cp_param[[rr]]$variable,"(period",1:(length(w_cp_param[[rr]]$changepoints)+1),")",sep=""))
        if((pointer_id+1)<=length(result_var_names)){
          result_var_names_new=c(result_var_names_new,result_var_names[(pointer_id+1):length(result_var_names)])
        }
      }
    }
    last_id=nrow(out_filter$m) # one more than nrow(data_ss_ori)
    # for fixed coefficient -> pick the last timepoint estimation
    # for AR(1) coefficient -> pick the average of 1/6 to 6/6 timepoints
    # for random walk -> pick the last timepoint estimation, since any point doesn't make sense
    pos=1
    for(a in 1:(length(formula_var)+1)){
      if(result_var_names[a] %in% AR1_coeffi){
        Estimate[pos]=mean(out_filter$m[round(last_id/6,0):last_id,a]) # if it's an AR(1) process, return the mean of 1/6 to 6/6
        pos=pos+1
      }else if(result_var_names[a] %in% rw_coeffi){
        Estimate[pos]=NA
        pos=pos+1
      }else if(result_var_names[a] %in% w_cp_param_variables){
        changes_id=c(w_cp_param[[which(result_var_names[a]==w_cp_param_variables)]][["changepoints"]]-10,last_id)
        Estimate[pos:(pos+length(changes_id)-1)]=out_filter$m[changes_id,a]
        Std.Error[pos:(pos+length(changes_id)-1)]=unlist(lapply(changes_id,function(temp_id){sqrt(diag(dlmSvd2var(out_filter$U.C[[temp_id]],out_filter$D.C[temp_id,])))[a]}))
        pos=pos+length(changes_id)
      }else{
        Estimate[pos]=out_filter$m[last_id,a]
        Std.Error[pos]=sqrt(diag(dlmSvd2var(out_filter$U.C[[last_id]],out_filter$D.C[last_id,])))[a]
        pos=pos+1
      }
    }
    result=as.data.frame(cbind(Estimate=Estimate,Std.Error=Std.Error))
    rownames(result)=result_var_names_new
    return(result)
  }
  get.smooth_result=function(out_filter,formula_var,AR1_coeffi,rw_coeffi,w_cp_param){
    out_smooth=dlmSmooth(out_filter)
    # organize result
    result_var_names=c("(Intercept)",formula_var)

    if(is.null(w_cp_param)){
      Estimate=Std.Error=rep(NA,length(formula_var)+1)
      w_cp_param_variables=NULL
      result_var_names_new=result_var_names
    }else{
      Estimate=Std.Error=rep(NA,length(formula_var)+1+length(unlist(lapply(1:length(w_cp_param),function(i){w_cp_param[[i]][["changepoints"]]}))))
      result_var_names_new=result_var_names
      for(rr in 1:length(w_cp_param)){
        pointer_id=grep(paste("^",w_cp_param[[rr]]$variable,"$",sep=""),result_var_names)
        temp_id=grep(paste("^",w_cp_param[[rr]]$variable,"$",sep=""),result_var_names_new)
        result_var_names_new=c(result_var_names_new[1:(temp_id-1)],
                               paste(w_cp_param[[rr]]$variable,"(period",1:(length(w_cp_param[[rr]]$changepoints)+1),")",sep=""))
        if((pointer_id+1)<=length(result_var_names)){
          result_var_names_new=c(result_var_names_new,result_var_names[(pointer_id+1):length(result_var_names)])
        }
      }
    }
    last_id=nrow(out_smooth$s) # one more than nrow(data_ss_ori)
    # for fixed coefficient -> pick the last timepoint estimation
    # for AR(1) coefficient -> pick the average of 1/6 to 6/6 timepoints
    # for random walk -> pick the last timepoint estimation, since any point doesn't make sense
    pos=1
    for(a in 1:(length(formula_var)+1)){
      if(result_var_names[a] %in% AR1_coeffi){
        Estimate[pos]=mean(out_smooth$s[round(last_id/6,0):last_id,a]) # if it's an AR(1) process, return the mean of 1/6 to 6/6
        pos=pos+1
      }else if(result_var_names[a] %in% rw_coeffi){
        Estimate[pos]=NA
        pos=pos+1
      }else if(result_var_names[a] %in% w_cp_param_variables){
        changes_id=c(w_cp_param[[which(result_var_names[a]==w_cp_param_variables)]][["changepoints"]]-10,last_id)
        Estimate[pos:(pos+length(changes_id)-1)]=out_smooth$s[changes_id,a]
        Std.Error[pos:(pos+length(changes_id)-1)]=unlist(lapply(changes_id,function(temp_id){sqrt(diag(dlmSvd2var(out_filter$U.C[[temp_id]],out_filter$D.C[temp_id,])))[a]}))
        pos=pos+length(changes_id)
      }else{
        Estimate[pos]=out_smooth$s[last_id,a]
        Std.Error[pos]=sqrt(diag(dlmSvd2var(out_filter$U.C[[last_id]],out_filter$D.C[last_id,])))[a]
        pos=pos+1
      }
    }
    result=as.data.frame(cbind(Estimate=Estimate,Std.Error=Std.Error))
    rownames(result)=result_var_names_new
    return(result)
  }
  missing_y_sample=function(data_temp,out_filter_temp,dlm_option="smooth",stepsize_for_newpart=1){
    out_smooth_temp=dlmSmooth(out_filter_temp)
    Fpart_temp=cbind(1,data_temp[,formula_var])
    if(!all(out_filter_temp$f-rowSums(Fpart_temp*out_filter_temp$a[,which(out_filter_temp$mod$FF[1,]==1)]) < 1E-10)){
      stop("The time-varying F matrix has not been updated correctly.")
    }
    if(dlm_option=="filter"){
      # check dimentions before multiplied together
      if(!all(dim(Fpart_temp)==dim(out_filter_temp$m[-1,which(out_filter_temp$mod$FF[1,]==1)]))){stop("Error in F and estimated coefficients dimentions")}
      init_y_temp=rowSums(Fpart_temp*out_filter_temp$m[-1,which(out_filter_temp$mod$FF[1,]==1)])[missing_part]
    }else if(dlm_option=="smooth"){
      # check dimentions before multiplied together
      if(!all(dim(Fpart_temp)==dim(out_smooth_temp$s[-1,which(out_filter_temp$mod$FF[1,]==1)]))){stop("Error in F and estimated coefficients dimentions")}
      init_y_temp=rowSums(Fpart_temp*out_smooth_temp$s[-1,which(out_filter_temp$mod$FF[1,]==1)])[missing_part]
    }
    init_y_new=init_y_old*(1-stepsize_for_newpart)+init_y_temp*stepsize_for_newpart
    return(init_y_new)
  }

  ######### initialization of data (part III of initialization) ##########
  # preparation part
  data_temp=data_ss_ori
  missing_part=which(is.na(data_temp$y))[which(is.na(data_temp$y))<nrow(data_temp)] # if the last y is missing, ignore it
  if(initial_imputation_option=="StructTS"){
    init_y_old=na_kalman(data_temp$y,model="StructTS")[missing_part]
  }else{stop("initial_imputation_option only accept StructTS for now." )}
  # fill in initial guess of y into y_1
  if(all(is.na(data_temp$y_1[missing_part+1]))){
    data_temp$y_1[missing_part+1]=init_y_old # [Initialization complete]
    if(length(grep("_and_",formula_var))==1){
      interaction_var=grep("_and_",formula_var,value=T)
      interaction_pair=unlist(strsplit(formula_var[grep("_and_",formula_var)],"_and_"))
      if(all(is.na(data_temp[missing_part+1,interaction_var]))){
        data_temp[missing_part+1,interaction_var]=data_temp[missing_part+1,interaction_pair[1]]*data_temp[missing_part+1,interaction_pair[2]]
      }else{
        stop("Missing value positions mismatch for interaction y_1_and_x.")
      }
    }
  }else{stop("Error in positions of NA in y_1.")}

  # Imputation Iteration
  # preparation -> I) SSM part -> II) Update part -> III) Tracking part -> IV)  Checking convergence part -> final output
  # ------------------- Imputation Iteration (begin)---------------------- #
  estimate_old=estimate_new=data.frame()
  lik_old=lik_new=c()
  iter=1
  convergence=cpt.convergence=F

  while(convergence==F){
    # SSM model fitting part
    # 1) estimate unknown variance parameter
    # 2) get estimations for time-varying coefficients, filter version or smooth version
    outMLE_temp=dlmMLE(data_temp$y,parm=inits,buildCamp_iterative)
    out_filter_temp=dlmFilter(data_temp$y,buildCamp_iterative(outMLE_temp$par))
    out_smooth_temp=dlmSmooth(out_filter_temp)

    # Update part
    # 1) based on previous round's F_t and current updated beta_t, update guess of Y_{t-1} into current round's F_t
    init_y_new=missing_y_sample(data_temp,out_filter_temp,dlm_option=dlm_option,stepsize_for_newpart=stepsize_for_newpart)
    data_temp$y_1[missing_part+1]=init_y_new  # update F_t using new y_1 info
    if(length(grep("_and_",formula_var))==1){ # update interaction y_1 and x
      data_temp[missing_part+1,interaction_var]=data_temp[missing_part+1,interaction_pair[1]]*data_temp[missing_part+1,interaction_pair[2]]
    }
    # 2) (optional) update changepoints in coefficients, given the current updated beta_t
    if(!is.null(w_cp_param)){
      # if fixed_cpts = False for all variables
      if(all(unlist(lapply(1:length(w_cp_param), function(cc5){w_cp_param[[cc5]]$fixed_cpts==F})))){
        cat(red("The changepoints are not all fixed, and need to be updated."))
        start_pt=floor(nrow(out_filter_temp$m)*burnin)
        w_cps_new=c()
        for(hh in 1:length(w_cp_param)){
          if(dlm_option=="filter"){
            temp=out_filter_temp$m[start_pt:nrow(out_filter_temp$m),which(colnames(out_filter_temp$mod$GG)==w_cp_param[[hh]]$variable)]
          }else if(dlm_option=="smooth"){
            temp=out_smooth_temp$s[start_pt:nrow(out_smooth_temp$s),which(colnames(out_filter_temp$mod$GG)==w_cp_param[[hh]]$variable)]
          }
          if("meanvar" %in% cpt_method){
            cpt_temp=cpt.meanvar(temp,penalty="Manual",Q=ss_param_temp$w_cp_param[[hh]]$segments-1,method="BinSeg")
            w_cps_new=c(w_cps_new,cpts(cpt_temp)+start_pt-1)
            if(printFlag){
              cat(red("The changepoints are:",w_cps_new,"\n"))
              plot(cpt_temp)
            }
          }
          if("mean" %in% cpt_method){
            cpt_temp=cpt.mean(temp,penalty="Manual",Q=ss_param_temp$w_cp_param[[hh]]$segments-1,method="BinSeg")
            w_cps_new=c(w_cps_new,cpts(cpt_temp)+start_pt-1)
            if(printFlag){
              cat(red("The changepoints are:"),w_cps_new,"\n")
              plot(cpt_temp)
            }
          }
        }
        w_cps_new=merge_closepoints(w_cps_new,band=mergeband)
        cat(red("The changepoints are:"),w_cps_new,"; ")
        if(length(w_cps_new)==w_cp_param[[1]]$segments-1){
          if(all(abs(w_cp_param[[1]][["changepoints"]]-w_cps_new)<cpts_convergence_cri)){
            cpt.convergence=T
          }else{
            cpt.convergence=F
          }
        }else{
          cpt.convergence=F
        }
        cat(red("changepoint convergence is"),cpt.convergence,"\n")
        for(ii in 1:length(w_cp_param)){
          w_cp_param[[ii]][["changepoints"]]=w_cps_new
        }
      }else if(all(unlist(lapply(1:length(w_cp_param), function(cc5){w_cp_param[[cc5]]$fixed_cpts==T})))){
        cpt.convergence=T
        if(printFlag){cat("All changepoints are fixed, and no need to update after fitting.\n")}
      }else{stop("Error in w_cp_param after fitting: fixed_cpts are not all T or all F.")}
      # Caution: the code here only test the lenth of the last list in the w_cp_param, which is problematic
    }

    # Tracking part
    # 1) likelihood
    # 2) estimates of time-varying coefficients beta_t, filter version or smooth version
    # 3) (optional) cpts in time-varying coefficients (finished in the update part simultaneously)
    lik_new=outMLE_temp$value
    estimate_filter_new=get.filter_result(out_filter=out_filter_temp,formula_var=formula_var,AR1_coeffi=AR1_coeffi,rw_coeffi=rw_coeffi,w_cp_param=w_cp_param)
    estimate_smooth_new=get.smooth_result(out_filter=out_filter_temp,formula_var=formula_var,AR1_coeffi=AR1_coeffi,rw_coeffi=rw_coeffi,w_cp_param=w_cp_param)
    if(dlm_option=="filter"){estimate_new=estimate_filter_new}else if(dlm_option=="smooth"){estimate_new=estimate_smooth_new}
    # print out results after this iteration
    if(T){
      cat("i=",iter,"\n",sep="")
      cat("likelihood:",lik_new,"\n")
      cat("Estimate:",estimate_new[,1],"\n")
      cat("Std.Error:", estimate_new[,2],"\n")
    }

    # Checking convergence part
    # if it's converged, stop here; otherwise, continue to next iteraction
    if(iter>1 & iter <= max_iteration){
      if(is.null(w_cp_param)){
        if(sum(abs(estimate_new[,1]-estimate_old[,1]),na.rm = T)< estimate_convergence_cri & abs(lik_new-lik_old) < lik_convergence_cri){
          convergence=T
        }
      }else{
        if(cpt.convergence==T){
          if(sum(abs(estimate_new[,1]-estimate_old[,1]),na.rm = T)< estimate_convergence_cri & abs(lik_new-lik_old) < lik_convergence_cri){
            convergence=T
          }
        }
      }
    }else if(iter > max_iteration){
      cat(red("Don't converge before max number of iteration."))
      return(NULL)
    }

    # save result of current round, for the convergence checking with the next round
    init_y_old=init_y_new
    estimate_old=estimate_new
    lik_old=lik_new
    iter=iter+1
  }

  # Save final result
  init_y_final=init_y_new # one iteration more than out_filter and out_smooth
  Fpart=cbind(1,data_temp[,formula_var]) # one iteration more than out_filter and out_smooth
  out_smooth=out_smooth_temp
  out_filter=out_filter_temp
  # ------------------- Imputation Iteration (end)---------------------- #

  # result1: without multiple imputation
  if(dlm_option=="filter"){
    result_convergence=get.filter_result(out_filter=out_filter,formula_var=formula_var,AR1_coeffi=AR1_coeffi,rw_coeffi=rw_coeffi,w_cp_param=w_cp_param)
  }else if(dlm_option=="smooth"){
    result_convergence=get.smooth_result(out_filter=out_filter,formula_var=formula_var,AR1_coeffi=AR1_coeffi,rw_coeffi=rw_coeffi,w_cp_param=w_cp_param)
  }
# Rubin's Rule
  pool.SSMimpute_mp=function(result){
    # calculate mean of coefficients
    name=rownames(result[[as.character(1)]])
    coefficients=matrix(unlist(lapply(1:m,function(i){result[[as.character(i)]][,1]})),byrow=T,ncol=length(name))
    colnames(coefficients)=name
    a=colMeans(coefficients)
    # calculate adjusted variance
    variances=matrix(unlist(lapply(1:m,function(i){(result[[as.character(i)]][,2])^2})),byrow=T,ncol=length(name))
    U_bar=colMeans(variances)
    B=0
    for(j in 1:m){
      B=B+(result[[as.character(j)]][,1]-a)%*%t(result[[as.character(j)]][,1]-a)
    }
    B_adjusted=diag(B)
    coefficients_sd=sqrt((1+1/m)*B_adjusted+U_bar)
    return(cbind(Estimate=a,Std.Error=coefficients_sd))
  }
  # result 2: multiple imputation for all variables
  if(!is.null(seed)){set.seed(seed);cat(red("seed for multiple imputation after SSMimpute is:",seed),"\n")}
  missing_y_1_mp=matrix(NA,ncol=m,nrow=length(init_y_final))
  mp_results=list()
  for(j1 in 1:m){
    if(dlm_option=="filter"){
      imputed_states=matrix(unlist(lapply(1:nrow(Fpart),function(j2){mvrnorm(n=1,mu=out_filter$m[j2+1,which(out_filter$mod$FF==1)],Sigma=dlmSvd2var(out_filter$U.C[[j2+1]],out_filter$D.C[j2+1,])[which(out_filter$mod$FF==1),which(out_filter$mod$FF==1)])})),
                            nrow=nrow(Fpart),ncol=ncol(Fpart),byrow=T)
    }else if(dlm_option=="smooth"){
      imputed_states=matrix(unlist(lapply(1:nrow(Fpart),function(j2){mvrnorm(n=1,mu=out_smooth$s[j2+1,which(out_filter$mod$FF==1)],Sigma=dlmSvd2var(out_smooth$U.S[[j2+1]],out_smooth$D.S[j2+1,])[which(out_filter$mod$FF==1),which(out_filter$mod$FF==1)])})),
                            nrow=nrow(Fpart),ncol=ncol(Fpart),byrow=T)
    }
    missing_y_1_mp[,j1]=rowSums(Fpart*imputed_states)[missing_part]
  }
  if(!is.null(w_cp_param)){
    if(w_cp_param[[1]]$fixed_cpts==T){
      w_cps_new=c()
      cat(red("Changepoints created to be skipped when they are given and fixed"))
      for(j3 in 1:length(w_cp_param)){
        w_cps_new=c(w_cps_new,w_cp_param[[j3]]$changepoints)
      }
    }
    cpts_region=unique(unlist(lapply(w_cps_new,function(kk){(kk-10):(kk+10)})))
    non_impute_region=missing_part[missing_part %in% cpts_region]
    data_temp_deletecp_y=data_temp$y
    data_temp_deletecp_y[non_impute_region+1]=NA
    for(j in 1:m){
      data_temp$y_1[missing_part+1]=missing_y_1_mp[,j]
      if(length(grep("_and_",formula_var))==1){
        data_temp[missing_part+1,interaction_var]=data_temp[missing_part+1,interaction_pair[1]]*data_temp[missing_part+1,interaction_pair[2]]
      }
      outMLE_temp=dlmMLE(data_temp_deletecp_y,parm=inits,buildCamp_iterative)  # M step
      out_filter_temp=dlmFilter(data_temp_deletecp_y,buildCamp_iterative(outMLE_temp$par))
      if(dlm_option=="filter"){
        mp_results[[as.character(j)]]=get.filter_result(out_filter=out_filter_temp,formula_var=formula_var,AR1_coeffi=AR1_coeffi,rw_coeffi=rw_coeffi,w_cp_param=w_cp_param)
      }else if(dlm_option=="smooth"){
        mp_results[[as.character(j)]]=get.smooth_result(out_filter=out_filter_temp,formula_var=formula_var,AR1_coeffi=AR1_coeffi,rw_coeffi=rw_coeffi,w_cp_param=w_cp_param)
      }
    }
  }else{
    for(j in 1:m){
      data_temp$y_1[missing_part+1]=missing_y_1_mp[,j]
      if(length(grep("_and_",formula_var))==1){
        data_temp[missing_part+1,interaction_var]=data_temp[missing_part+1,interaction_pair[1]]*data_temp[missing_part+1,interaction_pair[2]]
      }
      outMLE_temp=dlmMLE(data_temp$y,parm=inits,buildCamp_iterative)  # M step
      out_filter_temp=dlmFilter(data_temp$y,buildCamp_iterative(outMLE_temp$par))
      if(dlm_option=="filter"){
        mp_results[[as.character(j)]]=get.filter_result(out_filter=out_filter_temp,formula_var=formula_var,AR1_coeffi=AR1_coeffi,rw_coeffi=rw_coeffi,w_cp_param=w_cp_param)
      }else if(dlm_option=="smooth"){
        mp_results[[as.character(j)]]=get.smooth_result(out_filter=out_filter_temp,formula_var=formula_var,AR1_coeffi=AR1_coeffi,rw_coeffi=rw_coeffi,w_cp_param=w_cp_param)
      }
    }
  }
  print(mp_results)
  result_convergence_mp=as.data.frame(pool.SSMimpute_mp(mp_results))

  # result 3: multiple imputation with added V noise for all variables
  if(!is.null(seed)){set.seed(seed);cat(red("seed for multiple imputation after SSMimpute is:",seed),"\n")}
  missing_y_1_mp_addV=matrix(NA,ncol=m,nrow=length(init_y_final))
  mp_addV_results=list()
  for(j0 in 1:m){
    if(dlm_option=="filter"){
      imputed_states_addV=matrix(unlist(lapply(1:nrow(Fpart),function(j1){mvrnorm(n=1,mu=out_filter$m[j1+1,which(out_filter$mod$FF==1)],Sigma=dlmSvd2var(out_filter$U.C[[j1+1]],out_filter$D.C[j1+1,])[which(out_filter$mod$FF==1),which(out_filter$mod$FF==1)])})),
                                 nrow=nrow(Fpart),ncol=ncol(Fpart),byrow=T)
    }else if(dlm_option=="smooth"){
      imputed_states_addV=matrix(unlist(lapply(1:nrow(Fpart),function(j1){mvrnorm(n=1,mu=out_smooth$s[j1+1,which(out_filter$mod$FF==1)],Sigma=dlmSvd2var(out_smooth$U.S[[j1+1]],out_smooth$D.S[j1+1,])[which(out_filter$mod$FF==1),which(out_filter$mod$FF==1)])})),
                                 nrow=nrow(Fpart),ncol=ncol(Fpart),byrow=T)
    }
    missing_y_1_mp_addV[,j0]=rowSums(Fpart*imputed_states_addV)[missing_part]+rnorm(n=length(missing_part),mean=0,sd=sqrt(out_filter$mod$V))
  }
  if(!is.null(w_cp_param)){
    if(w_cp_param[[1]]$fixed_cpts==T){
      w_cps_new=c()
      cat(red("Changepoints created to be skipped when they are given and fixed"))
      for(j2 in 1:length(w_cp_param)){
        w_cps_new=c(w_cps_new,w_cp_param[[j2]]$changepoints)
      }
    }
    cpts_region=unique(unlist(lapply(w_cps_new,function(kk){(kk-10):(kk+10)})))
    non_impute_region=missing_part[missing_part %in% cpts_region]
    data_temp_deletecp_y=data_temp$y
    data_temp_deletecp_y[non_impute_region+1]=NA
    for(j in 1:m){
      data_temp$y_1[missing_part+1]=missing_y_1_mp_addV[,j]
      if(length(grep("_and_",formula_var))==1){
        data_temp[missing_part+1,interaction_var]=data_temp[missing_part+1,interaction_pair[1]]*data_temp[missing_part+1,interaction_pair[2]]
      }
      outMLE_temp=dlmMLE(data_temp_deletecp_y,parm=inits,buildCamp_iterative)  # M step
      out_filter_temp=dlmFilter(data_temp_deletecp_y,buildCamp_iterative(outMLE_temp$par))
      if(dlm_option=="filter"){
        mp_addV_results[[as.character(j)]]=get.filter_result(out_filter=out_filter_temp,formula_var=formula_var,AR1_coeffi=AR1_coeffi,rw_coeffi=rw_coeffi,w_cp_param=w_cp_param)
      }else if(dlm_option=="smooth"){
        mp_addV_results[[as.character(j)]]=get.smooth_result(out_filter=out_filter_temp,formula_var=formula_var,AR1_coeffi=AR1_coeffi,rw_coeffi=rw_coeffi,w_cp_param=w_cp_param)
      }
    }
  }else{
    for(j in 1:m){
      data_temp$y_1[missing_part+1]=missing_y_1_mp_addV[,j]
      if(length(grep("_and_",formula_var))==1){
        data_temp[missing_part+1,interaction_var]=data_temp[missing_part+1,interaction_pair[1]]*data_temp[missing_part+1,interaction_pair[2]]
      }
      outMLE_temp=dlmMLE(data_temp$y,parm=inits,buildCamp_iterative)  # M step
      out_filter_temp=dlmFilter(data_temp$y,buildCamp_iterative(outMLE_temp$par))
      if(dlm_option=="filter"){
        mp_addV_results[[as.character(j)]]=get.filter_result(out_filter=out_filter_temp,formula_var=formula_var,AR1_coeffi=AR1_coeffi,rw_coeffi=rw_coeffi,w_cp_param=w_cp_param)
      }else if(dlm_option=="smooth"){
        mp_addV_results[[as.character(j)]]=get.smooth_result(out_filter=out_filter_temp,formula_var=formula_var,AR1_coeffi=AR1_coeffi,rw_coeffi=rw_coeffi,w_cp_param=w_cp_param)
      }
    }
  }
  print(mp_addV_results)
  result_convergence_mp_addV=as.data.frame(pool.SSMimpute_mp(mp_addV_results))

  if(printFlag){
    cat("The converged result is:\n")
    print(result_convergence)
    cat("The converged result after multiple imputation is:\n")
    print(result_convergence_mp)
    cat("The converged result after multiple imputation with added V noise is:\n")
    print(result_convergence_mp_addV)
  }
  # multiple imputation for AR(1) and rw coefficient (end)

  # organize result
  if(!is.null(w_cp_param) & all(unlist(lapply(1:length(w_cp_param), function(cc5){w_cp_param[[cc5]]$fixed_cpts==F})))){
    estimated_cpts=w_cps_new
  }else{
    estimated_cpts=NULL
  }
  return(list(result_convergence=result_convergence,
              result_convergence_mp=result_convergence_mp,
              result_convergence_mp_addV=result_convergence_mp_addV,
              estimated_cpts=estimated_cpts,
              out_filter=out_filter,
              iter=iter,
              data_temp = data_temp, y_final = init_y_final))
}
