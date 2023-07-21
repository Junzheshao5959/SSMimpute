# SSM: statespace model without imputation and without missing data in covariates
# first wrote on  2023.05.18
# This function combines previous function of "run.SSM_unanimous_cpts" and "run.SSM_partial_converged_cpts"
# Explanation: 0. no missing data in response and explanatory variables, and thus no imputation needed
#                 1) no need for <lik_convergence_cri>: convergence only depends on cpt, if exists
#                 2) no missing imputation: <cpt_guess_imputation_option> and <dlm_imputation_option> no long need
#                 3) no need of <m>: no imputation, no m for multiple imputation.
#              1. different from "run.SSM_unanimous_cpts" (previous version), when un-given change points are learnt from cpt.meanvar
#                 1) it's learn by the given cpt_method, which choose from "mean" or "meanvar"
#                 2) change points for different variables can be different, only very similarly ones are forced to be the same
#              2. different from "run.SSM_unanimous_cpts" and "run.SSM_partial_converged_cpts",
#                 1) formula_var is removed in the parameter list and replaced with formula
#                 2) formula_var is not given, but learned from "formula". 
#                    formula_var is also renamed as explanatory_var except two places in "run.SSM_unanimous_cpts"
#                 4) when learned change points are either forced to be same or learned separately (05/16/2023)
#                    -> it's chosen by cpt_merge_option for "unanimous" or "separate"
run.SSM=function(data_ss,formula,ss_param_temp,max_iteration=100,
                 cpt_learning_param=list(cpt_method="mean",burnin=1/10,mergeband=20,convergence_cri=10),
                 cpt_merge_option = "unanimous",
                 dlm_option="smooth",
                 dlm_cpt_learning_option="smooth",bandwidth=20, cpt_V=5, printFlag=T){
  {
    # !! Please read this paragraph first !!
    # Limitation 1: this current version of function has only tested on models that are no more complicated than y ~ y_1 + y_2 + ... + y_1*x + x + x_1 + c
    #               for a continuous y, a binary/continuous x, and a continuous c,
    #               other scenarios, including (i) interactions except for y_1 with x/c and (ii) binary y and c, are not allowed and haven't been tested.
    # [Improvement  on limitation 2]
    # Limitation 2: allow choosing "unanimous" change points or "separate" change points "cpt_merge_option"
    #               for "unanimous" option, change points are same across variables, except when they are all fixed in the very first beginning
    #               for "separate" option, different variable can have different number of change points at different locations, but change points closer 
    #                                      to each other across variables are forced to merge as the same
    #               Set-up works for  (i)  an initial collection of change points are given for a collection of variables,
    #                                      since those change points are fixed, and no more updating
    #                                (ii)  an initial collection of change points are given for a collection of variables,
    #                                      but those change points are not for sure, and need to be updated
    #                                (iii) no initial guess is given for change points, but only how many
    #                                      ["unanimous" option] the change points are learnt simultaneously and forced to be the same across all variables that are periodic stable
    #                                      ["separate" option] the change points are learnt separately for each variables, but forces to be the same if certain change points are very close.
    #                Working flowchart for [changepoints]
    #                (i) if change points are given in <changepoints> and <fixed_cpts>=T for all variables in [w_cp_param] 
    #                     -> 1) no initialization
    #                        2) no estimation for the location of change points
    #                        3) no iteration until convergence for change points 
    #                     we allow change points differ for each variables, as they are all given in advance in this case
    #                (ii) if change points are given in <changepoints> and <fixed_cpts>=F for all variables in [w_cp_param] 
    #                     -> 1) no initialization
    #                        2) estimate the location of change points for those variables with fixed_cpts=F
    #                        3) update change points via iteration, until convergence achieved for those variables with fixed_cpts=F
    #                     ["unanimous"] change points should be same for variables with fixed_cpts=F
    #                     ["separate"]  change points could be different for variables with fixed_cpts=F
    #                (iii) if change points are not given and thus <fixed_cpts> do not exist in [w_cp_param] 
    #                    -> 1) initialization
    #                       2) require estimation for the location of change points
    #                       3) require iteration until convergence for change points 
    #                     ["unanimous"] change points should be same for all periodic variables 
    #                     ["separate"]  change points could be different for all periodic variables 
    #                     changepoints could be different for variables for with fixed_cpts=F
    #                Otherwise, all other cases are not permitted for the moment, including
    #                    -> [not allowed] 1) changelings are partially given, and partially learnt
    #                    -> [not allowed] 2) <fixed_cpts> are partially T
    # Limitation 4: [to be changed]: simple and no iterative "StructTS" imputation for learning cpts in V
    # Limitation 5: estimated_cpts and iter in the return result only reflects those for W, not for V
    #
    # Explanation: 0. no missing data in response and explanatory variables, and thus no imputation needed
    #                 1) no need for <lik_convergence_cri>: convergence only depends on cpt, if exists
    #                 2) no missing imputation: <cpt_guess_imputation_option> and <dlm_imputation_option> no long need
    #                 3) no need of <m>: no imputation, no m for multiple imputation.
    #              1. different from "run.SSM_unanimous_cpts" (previous version), when un-given change points are learnt from cpt.meanvar
    #                 1) it's learn by the given cpt_method, which choose from "mean" or "meanvar"
    #                 2) change points for different variables can be different, only very similarly ones are forced to be the same
    #              2. different from "run.SSM_unanimous_cpts" and "run.SSM_partial_converged_cpts",
    #                 1) formula_var is removed in the parameter list and replaced with formula
    #                 2) formula_var is not given, but learned from "formula". 
    #                    formula_var is also renamed as explanatory_var except two places in "run.SSM_unanimous_cpts"
    #                 4) when learned change points are either forced to be same or learned separately (05/16/2023)
    #                    -> it's chosen by cpt_merge_option for "unanimous" or "separate"
    #
    # Explanation for parameters:
    # <data_ss>: contains all information, and only selected variables in formula_var enters the statespace model
    # <formula>
    # <ss_param>: 
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
    # <max_iteration>: control for the convergence of changepoints, a positive integer
    # <cpt_learning_param>:  <cpt_method> either "mean" or "meanvar"
    #                        <burnin> a positive number in (0,1)
    #                        <mergeband> a positive integer
    #                        <convergence_cri> a positive integer
    #                        Caution: all used for learning changepoints in coefficients, not for changepoints in observational variance
    # <cpt_merge_option>: "unanimous"(default) or "separate"
    # <dlm_option>: choose between "smooth" or "filter"
    # <dlm_cpt_learning_option>: choose smooth or filter results for initial cpt learning
    # <bandwidth>: bandwidth with increase V to allow quick adjust to new level
    # <cpt_V>
  }
  # extract response and explanatory variables from the formula
  # Question: how to specify interaction?
  response_var=unlist(strsplit(formula,"~"))[1]
  explanatory_var=unlist(strsplit(unlist(strsplit(formula,"~"))[2],"+",fixed=T))
  n_var=length(explanatory_var)+1
  
  # function to merge similar change points
  select_closest=function(ps_collection,old_ps){
    new_ps=rep(NA,length(old_ps))
    for(ps in 1:length(old_ps)){
      new_ps[ps]=ps_collection[which(abs(old_ps[ps]-ps_collection)==min(abs(old_ps[ps]-ps_collection)))]
    }
    return(new_ps)
  }
  
  # ------------------- Preparation ---------------------- #
  # check: 1. all variables specified in the model exist the data_ss set
  #        2. ss_param is a list
  #           1) m0 is either NULL(use the default in dlm), or given (contains no NA, all numeric, same length as states)
  #           2) C0 is either NULL(use the default in dlm), or given (contains no NA, all numeric, same dim as n_states*n_states)
  #           3) the variables, whose coefficients is AR(1) process, are char string, no NA, and in both model and data_ss set
  #           4) the variables, whose coefficients is random walk process, is char string and in both model and data_ss set
  #           5) check if change points of V need to be estimated 
  #              -> as no imputation is applied for y -> no change on data_ss -> no iteration needed
  #              In short, change points of V is only estimated once, in the initialization part
  #              Otherwise, if change points are given, no action is needed.
  #           6) check cpt_learning_param under w_param
  #           7) check if change points of W for certain variables need to be estimated
  #           8) check if the length of inits is correct 
  #         3. <max_iteration>: a positive integer
  if(!all(c(response_var,explanatory_var) %in% colnames(data_ss))){
    stop("The variables specified in the formula are not contained in the data_ss set.")
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
    v_cp_param=ss_param_temp$v_cp_param    # change points estimated once in the beginning, and add #(segments-1) NA for V -> checked
    
    # n_states includes: 1 for intercept, 
    #                    n_variable for all selected variable, 
    #                    length(AR1_coeffi) for additional baselines for the AR(1) process
    #              Note: random walk only change W, and do not add states
    #                    shift in beta only change W to 10, do not add anything
    #                    varying V only change segments of V_t, do not add states
    n_states=1+length(explanatory_var)+length(AR1_coeffi) 
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
      if(!all(!is.na(AR1_coeffi))){stop("Error in AR1_coeffi: AR1_coeffi contains NA.")}
      if(!all(AR1_coeffi %in% explanatory_var)){stop("Error in AR1_coeffi: the selected time-varying coefficient has no corresponding variable in the regression formula.")}
      if(!all(AR1_coeffi %in% colnames(data_ss))){stop("Error in AR1_coeffi: the selected time-varying coefficient has no corresponding variable in the data_ss set.")}
    }
    if(!is.null(rw_coeffi)){ # variables are char string, no NA, and in both model and data_ss set
      if(!is.character(rw_coeffi)){stop("Error in rw_coeffi: the specified variables in rw_coeffi are not all char string")}
      if(!all(!is.na(rw_coeffi))){stop("Error in rw_coeffi: rw_coeffi contains NA.")}
      if(!all(rw_coeffi %in% c("intercept",explanatory_var))){stop("Error in rw_coeffi: the selected time-varying coefficient has no corresponding variable in the regression formula.")}
      if(!all(rw_coeffi %in% c("intercept",colnames(data_ss)))){stop("Error in rw_coeffi: the selected time-varying coefficient has no corresponding variable in the data_ss set.")}
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
        if(!all(v_cp_param_learncps$changepoints > 1 & v_cp_param_learncps$changepoints < nrow(data_ss))){stop("Error in v_cp_param: [changepoints] is not range of [1,nrow(data_ss)-1].")}
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
        if(sum(is.na(data_ss[,response_var]))!=0){
          # for the ignore case
          cpt_v_temp=cpt.var(na_kalman(data_ss[,response_var],model="StructTS"),penalty="Manual",Q=v_cp_param_learncps$segments-1, method="BinSeg")
        }else{
          # for the full and cc case
          cpt_v_temp=cpt.var(data_ss[,response_var],penalty="Manual",Q=v_cp_param_learncps$segments-1, method="BinSeg")
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
          if(!all(cpt_learning_param$cpt_method %in% c("mean","meanvar"))){stop("Error in cpt_learning_param$cpt_method: should be mean or meanvar.")}
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
    }
    if(!is.null(w_cp_param)){
      w_cp_param_learncps=w_cp_param
      # Check: 1) The variables with changepoints must be in data_ss and in formula
      #           the order of variables aligns with the order in the formula
      #        2) for each variable with changepoints, 
      #           check 2.1) [variable], [segments], [changepoints], [fixed_cpts] may exist
      #                 2.2) #elements for each sublist is 1-4.
      #                 2.3) [variable] and [segments] must exist
      #                      [variable] must be a char
      #                      [segments] must be a positive number
      #                 2.4) when [changepoints] exists -> [fixed_cpts] must exist and be logical
      #                 2.5) when [changepoints] doesn't exist, it need to be estimated only once from the data_ss 
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
      if(cpt_merge_option == "unanimous"){
        if(length(unique(unlist(lapply(1:length(w_cp_param_learncps),function(aa1){w_cp_param_learncps[[aa1]][["segments"]]}))))>1){
          stop("Error in w_cp_param: the segments value for different variables need to be the same")
        }
      }
      # check 1)
      w_cp_param_variables=unlist(lapply(1:length(w_cp_param_learncps),function(bb){w_cp_param_learncps[[bb]][["variable"]]}))
      if(!all(w_cp_param_variables %in% explanatory_var)){stop("The selected level-shifted coefficient has no corresponding variable in the regression formula.")}
      if(!all(w_cp_param_variables %in% colnames(data_ss))){stop("The selected level-shifted coefficient has no corresponding variable in the data_ss set.")}
      if(length(w_cp_param_variables)>1){
        w_cp_param_variables_order=c()
        for(bb in w_cp_param_variables){
          w_cp_param_variables_order=c(w_cp_param_variables_order,which(bb == c("intercept",explanatory_var)))
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
            all(w_cp_param_learncps[[cc3]]$changepoints>1 & w_cp_param_learncps[[cc3]]$changepoints < nrow(data_ss))})))){
          stop("Error in w_cp_param: [changepoints] for all variables are not all numeric and within the right range.")
        }
        if(cpt_merge_option == "unanimous"){
          if(length(unique(unlist(lapply(1:length(w_cp_param_learncps),function(aa2){w_cp_param_learncps[[aa2]][["changepoints"]]}))))!=w_cp_param_learncps[[1]]$segments-1){
            stop("Error in w_cp_param: given [changepoints] may not be same for all variables")
          }
        }
        if(all(unlist(lapply(1:length(w_cp_param_learncps), function(cc4){w_cp_param_learncps[[cc4]]$fixed_cpts})))){
          cat(red("Note: ''changepoints'' of W is given and fixed.\n"))
        }else if(all(unlist(lapply(1:length(w_cp_param_learncps), function(cc4){!w_cp_param_learncps[[cc4]]$fixed_cpts})))){
          cat(red("Note: ''changepoints'' of W is given a initial guess, but need to be updated according to the fitting result.\n"))
        }else{
          stop("Error in w_cp_param: [fixed_cpts] need to be all T or all F when changepoints exist.")
        }
      }else if(all(unlist(lapply(1:length(w_cp_param_learncps), function(cc2){ (!"changepoints" %in% names(w_cp_param_learncps[[cc2]])) & (!"fixed_cpts" %in% names(w_cp_param_learncps[[cc2]]))})))){
        # if [changepoints] and [fixed_cpts] both do not exist for each variable
        if(printFlag){cat("The changepoints of variables, whose coefficients shift over time, are not given.\n")}
        ######### initialization of changepoints of W (part II of initialization) ##########        
        # if the changepoints don't exist, we ask the corresponding coeffi be randome walk, and then learn changepoints from the results
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
        # fit dlm model
        # [Unique point] (begin)
        data_ss_temp_init_cpts_w=data_ss
        out_filter_init_cpts_w=run.SSM(data_ss=data_ss_temp_init_cpts_w,formula=formula,ss_param_temp=ss_param_learncps,max_iteration=max_iteration,
                                       cpt_learning_param=cpt_learning_param,
                                       cpt_merge_option=cpt_merge_option,
                                       dlm_option=dlm_option,
                                       dlm_cpt_learning_option=dlm_cpt_learning_option,bandwidth=bandwidth,cpt_V=cpt_V,printFlag=printFlag)$out_filter
        out_smooth_init_cpts_w=dlmSmooth(out_filter_init_cpts_w)
        # [Unique point] (end)
        # get the change points
        w_cps_new=c()
        w_cp_param_new=list()
        start_pt=floor(nrow(out_filter_init_cpts_w$m)*burnin)
        for(ee in 1:length(w_cp_param)){
          if(dlm_cpt_learning_option=="filter"){
            temp=out_filter_init_cpts_w$m[start_pt:nrow(out_filter_init_cpts_w$m),which(colnames(out_filter_init_cpts_w$mod$GG)==w_cp_param[[ee]]$variable)]
          }else if(dlm_cpt_learning_option=="smooth"){
            temp=out_smooth_init_cpts_w$s[start_pt:nrow(out_smooth_init_cpts_w$s),which(colnames(out_filter_init_cpts_w$mod$GG)==w_cp_param[[ee]]$variable)]
          }
          if("meanvar" %in% cpt_method){
            cpt_temp=cpt.meanvar(temp,penalty="Manual",Q=ss_param_temp$w_cp_param[[ee]]$segments-1,method="BinSeg") 
            w_cps_new=c(w_cps_new,cpts(cpt_temp)+start_pt-1) 
            w_cp_param_new[[ee]]=cpts(cpt_temp)+start_pt-1
            if(printFlag){
              plot(cpt_temp,ylab=w_cp_param[[ee]]$variable)
            }
          }
          if("mean" %in% cpt_method){
            cpt_temp=cpt.mean(temp,penalty="Manual",Q=ss_param_temp$w_cp_param[[ee]]$segments-1,method="BinSeg")
            w_cps_new=c(w_cps_new,cpts(cpt_temp)+start_pt-1)
            w_cp_param_new[[ee]]=cpts(cpt_temp)+start_pt-1
            if(printFlag){
              plot(cpt_temp,ylab=w_cp_param[[ee]]$variable)
            }
          }
        }
        w_cps_new=merge.closepoints(points=w_cps_new,band=mergeband)
        for(ff in 1:length(w_cp_param)){
          if(cpt_merge_option == "unanimous"){
            w_cp_param[[ff]][["changepoints"]]=w_cps_new
          }else if(cpt_merge_option == "separate"){
            w_cp_param[[ff]][["changepoints"]]=select_closest(w_cps_new,w_cp_param_new[[ff]])
          }
          w_cp_param[[ff]][["fixed_cpts"]]=F
        }
        if(printFlag){
          cat(red("The inital guess of the collection of changepoints for all variables are:"),w_cps_new,"\n")
          cat(red("The initial guess of changepoints for each variables are:\n"))
          names(w_cp_param_new)=w_cp_param_variables
          print(w_cp_param_new)
          cat(red("The updated parameters for changepoints are:\n"))
          print(w_cp_param)
        }
      }else{
        stop("Error in w_cp_param: [changepoints] and [fixed_cpts] must both exist or not exist for all variables")
      }
      # final check 2.6): <variable>, <segments> and <changpoints> all exist for each variable
      for(gg in 1:length(w_cp_param)){
        if(!all(c("variable","segments","changepoints","fixed_cpts") %in% names(w_cp_param[[gg]]))){
          stop("<variable>, <segments>, <changepoints>, and <fixed_cpts> must exist for each level shift coefficient.")
        }
        if(length(w_cp_param[[gg]]$changepoints)!=w_cp_param[[gg]]$segments-1){stop("length of changepoints must be segments-1.")}
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
  if(!(length(max_iteration)==1 & is.numeric(max_iteration) & max_iteration>0 & max_iteration%%1==0)){stop("max_iteration is not a positive integer!")}
  if(!(length(dlm_cpt_learning_option)==1 & is.character(dlm_cpt_learning_option) & dlm_cpt_learning_option %in% c("smooth","filter"))){stop("Error in dlm_cpt_learning_option: choose between smooth or filter.")}
  if(!(dlm_option %in% c("smooth","filter") & is.character(dlm_option) & length(dlm_option)==1)){stop("dlm_option can either be filter or smooth")}
  # functions to build up statespace model, which certains terms to be estimates as NA
  # the items in u are: [1:length(AR1_coeffi)] terms -> logit(autocorrelation) for all AR(1) coefficients
  #                                                     default value is 0=exp(0.5)/(1+exp(0.5))
  #                     [length(AR1_coeffi)+1 : 2*length(AR1_coeffi)] -> log(variance for Wj) for all AR(1) coefficients
  #                     [2*length(AR1_coeffi)+1 : 2*length(AR1_coeffi)+length(rw_coeffi)] -> log(variance for Wj) for random walk coefficients
  #                     [2*length(AR1_coeffi)+length(rw_coeffi)+1 : 2*length(AR1_coeffi)+length(rw_coeffi)+length(v_changepoint)] -> (possibly sereis of) log(variance for V)
  #                     [2*length(AR1_coeffi)+1 : 3*length(AR1_coeffi)] -> log(variance for Wj) for random walk coefficients
  buildCamp=function(u){
    model=dlmModReg(data_ss[,explanatory_var],dV=exp(u[2*length(AR1_coeffi)+length(rw_coeffi)+1]))
    FF=matrix(c(model$FF,rep(0,length(AR1_coeffi))),nrow=1)
    if(length(AR1_coeffi)>0){
      states_name=c("intercept",explanatory_var,paste("rho.",AR1_coeffi,sep=""))
    }else{
      states_name=c("intercept",explanatory_var)
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
        jw_cps_temp=unique(unlist(lapply(w_cp_param[[jw]]$changepoints,function(jwpoint){(jwpoint-bandwidth):(jwpoint+bandwidth)})))
        jw_cps_temp=jw_cps_temp[jw_cps_temp>0 & jw_cps_temp<=nrow(model$X)]
        w_changepoint[jw_cps_temp,w_cp_param_variables[jw]]=cpt_V
      }
      X(model)=cbind(model$X,w_changepoint)
      JW(model)=JW
    }
    dlm(model)
    return(model)
  }
  get.filter_result=function(out_filter,explanatory_var,AR1_coeffi,rw_coeffi,w_cp_param){
    # organize result
    result_var_names=c("(Intercept)",explanatory_var)
    
    if(is.null(w_cp_param)){
      Estimate=Std.Error=rep(NA,length(explanatory_var)+1)
      w_cp_param_variables=NULL
      result_var_names_new=result_var_names
    }else{
      Estimate=Std.Error=rep(NA,length(explanatory_var)+1+length(unlist(lapply(1:length(w_cp_param),function(i){w_cp_param[[i]][["changepoints"]]}))))
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
    last_id=nrow(out_filter$m) # one more than nrow(data_ss_new)
    # for fixed coefficient -> pick the last timepoint estimation
    # for AR(1) coefficient -> pick the average of 1/6 to 6/6 timepoints
    # for random walk -> pick the last timepoint estimation, since any point doesn't make sense
    pos=1
    for(a in 1:(length(explanatory_var)+1)){
      if(result_var_names[a] %in% AR1_coeffi){
        Estimate[pos]=mean(out_filter$m[round(last_id/6,0):last_id,a]) # if it's an AR(1) process, return the mean of 1/6 to 6/6
        pos=pos+1
      }else if(result_var_names[a] %in% rw_coeffi){
        Estimate[pos]=NA
        pos=pos+1
      }else if(result_var_names[a] %in% w_cp_param_variables){
        changes_id=c(w_cp_param[[which(result_var_names[a]==w_cp_param_variables)]][["changepoints"]]-bandwidth,last_id)
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
  get.smooth_result=function(out_filter,explanatory_var,AR1_coeffi,rw_coeffi,w_cp_param){
    out_smooth=dlmSmooth(out_filter)
    # organize result
    result_var_names=c("(Intercept)",explanatory_var)
    
    if(is.null(w_cp_param)){
      Estimate=Std.Error=rep(NA,length(explanatory_var)+1)
      w_cp_param_variables=NULL
      result_var_names_new=result_var_names
    }else{
      Estimate=Std.Error=rep(NA,length(explanatory_var)+1+length(unlist(lapply(1:length(w_cp_param),function(i){w_cp_param[[i]][["changepoints"]]}))))
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
    last_id=nrow(out_smooth$s) # one more than nrow(data_ss)
    # for fixed coefficient -> pick the last timepoint estimation
    # for AR(1) coefficient -> pick the average of 1/6 to 6/6 timepoints
    # for random walk -> pick the last timepoint estimation, since any point doesn't make sense
    pos=1
    for(a in 1:(length(explanatory_var)+1)){
      if(result_var_names[a] %in% AR1_coeffi){
        Estimate[pos]=mean(out_smooth$s[round(last_id/6,0):last_id,a]) # if it's an AR(1) process, return the mean of 1/6 to 6/6
        pos=pos+1
      }else if(result_var_names[a] %in% rw_coeffi){
        Estimate[pos]=NA
        pos=pos+1
      }else if(result_var_names[a] %in% w_cp_param_variables){
        changes_id=c(w_cp_param[[which(result_var_names[a]==w_cp_param_variables)]][["changepoints"]]-bandwidth,last_id)
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
  
  # ------------------- main part (begin)---------------------- #  
  # 1. given data without missing values in the covariates -> fit state space model
  # 2. if cpts is not given, alternate between i)  learn cpts
  #                                            ii) refit state sapce model
  #    until convergence
  
  # 1. find MLE for init, build up spacestate model, and get filter states
  if(printFlag){
    cat("==============================================================================\n")
    cat("Build statespace model with unknown variance and auto-correlation terms estimated.\n")
    cat("\n The initial m0, C0, and inits given by the users is:\n m0: ",m0,
        "\n C0: ",C0,
        "\n inits: ",inits, 
        ", or after transformation is: ")
    if(length(AR1_coeffi)>0){
      cat(c(exp(inits[1:length(AR1_coeffi)])/(1+exp(inits[1:length(AR1_coeffi)])),exp(inits[(length(AR1_coeffi)+1):length(inits)])))
    }else{
      cat(exp(inits[(length(AR1_coeffi)+1):length(inits)]))
    }
    cat("\n")
  }
  outMLE=dlmMLE(data_ss[,response_var],parm=inits,buildCamp) 
  out_filter=dlmFilter(data_ss[,response_var],buildCamp(outMLE$par))
  out_smooth=dlmSmooth(out_filter)
  
  iter=1
  # added for change points(begin)
  if(!is.null(w_cp_param)){
    # if fixed_cpts = False for all variables
    if(all(unlist(lapply(1:length(w_cp_param), function(cc5){w_cp_param[[cc5]]$fixed_cpts==F})))){
      cat(red("The changepoints are all not fixed, and need to be updated."))
      out_filter_temp=out_filter
      out_smooth_temp=out_smooth
      start_pt=floor(nrow(out_filter_temp$m)*burnin)
      convergence=F
      while(convergence==F){
        w_cps_new=c()
        w_cp_param_new=list()
        for(hh in 1:length(w_cp_param)){
          if(dlm_cpt_learning_option=="filter"){
            temp=out_filter_temp$m[start_pt:nrow(out_filter_temp$m),which(colnames(out_filter_temp$mod$GG)==w_cp_param[[hh]]$variable)]
          }else if(dlm_cpt_learning_option=="smooth"){
            temp=out_smooth_temp$s[start_pt:nrow(out_smooth_temp$s),which(colnames(out_filter_temp$mod$GG)==w_cp_param[[hh]]$variable)]
          }
          if("meanvar" %in% cpt_method){
            cpt_temp=cpt.meanvar(temp,penalty="Manual",Q=ss_param_temp$w_cp_param[[hh]]$segments-1,method="BinSeg")
            w_cps_new=c(w_cps_new,cpts(cpt_temp)+start_pt-1)
            w_cp_param_new[[hh]]=cpts(cpt_temp)+start_pt-1
            if(printFlag){
              cat(red("The changepoints are:",w_cps_new,"\n"))
              plot(cpt_temp)
            }
          }
          if("mean" %in% cpt_method){
            cpt_temp=cpt.mean(temp,penalty="Manual",Q=ss_param_temp$w_cp_param[[hh]]$segments-1,method="BinSeg")
            w_cps_new=c(w_cps_new,cpts(cpt_temp)+start_pt-1)
            w_cp_param_new[[hh]]=cpts(cpt_temp)+start_pt-1
            if(printFlag){
              cat(red("The changepoints are:"),w_cps_new,"\n")
              plot(cpt_temp)
            }
          }
        }
        w_cps_new=merge.closepoints(w_cps_new,band=mergeband)
        cat(red("The changepoints collection are:"),w_cps_new,";\n")
        for(ii in 1:length(w_cp_param)){
          w_cp_param_new[[ii]]=select_closest(w_cps_new,w_cp_param_new[[ii]])
        }
        print(w_cp_param_new)
        if(cpt_merge_option=="unanimous"){
          if(length(w_cps_new)==w_cp_param[[1]]$segments-1){
            if(all(abs(w_cp_param[[1]][["changepoints"]]-w_cps_new)<cpts_convergence_cri)){
              convergence=T
            }else{
              convergence=F
            }
          }else{
            convergence=F
          }
        }else if(cpt_merge_option=="separate"){
          cpt.convergence_collection=rep(F,length(w_cp_param))
          for(ii in 1:length(w_cp_param)){
            if(length(w_cp_param_new[[ii]])==w_cp_param[[ii]]$segments-1){
              if(all(abs(w_cp_param_new[[ii]]-w_cp_param[[ii]]$changepoints)<cpts_convergence_cri)){
                cpt.convergence_collection[ii]=T
              }
            }
          }
          if(all(cpt.convergence_collection==T)){
            convergence=T
          }else{
            convergence=F
          }
        }
        cat(red("changepoint convergence is"),convergence,"\n")
        for(ii in 1:length(w_cp_param)){
          if(cpt_merge_option=="unanimous"){
            w_cp_param[[ii]][["changepoints"]]=w_cps_new
          }else if(cpt_merge_option=="separate"){
            w_cp_param[[ii]][["changepoints"]]=w_cp_param_new[[ii]]
          }
        }
        
        outMLE_temp=dlmMLE(data_ss[,response_var],parm=inits,buildCamp) 
        out_filter_temp=dlmFilter(data_ss[,response_var],buildCamp(outMLE_temp$par))
        
        iter=iter+1
        # check if it never converges
        if(convergence==F & iter > max_iteration){
          cat(red("Error in convergence."))
          return(NULL)
        }
      }
      out_filter=out_filter_temp
    }else if(all(unlist(lapply(1:length(w_cp_param), function(cc5){w_cp_param[[cc5]]$fixed_cpts==T})))){
      if(printFlag){cat("All changepoints are fixed, and no need to update after fitting.\n")}
    }else{stop("Error in w_cp_param after fitting: fixed_cpts are not all T or all F.")}
  }
  # added for changepoints (end)
  
  if(dlm_option=="smooth"){
    result=get.smooth_result(out_filter=out_filter,explanatory_var=explanatory_var,AR1_coeffi=AR1_coeffi,rw_coeffi=rw_coeffi,w_cp_param=w_cp_param)
  }else if(dlm_option=="filter"){
    result=get.filter_result(out_filter=out_filter,explanatory_var=explanatory_var,AR1_coeffi=AR1_coeffi,rw_coeffi=rw_coeffi,w_cp_param=w_cp_param)
  }
  if(printFlag){print(result)}
  
  # organize result of change points
  if(!is.null(w_cp_param) & all(unlist(lapply(1:length(w_cp_param), function(cc5){w_cp_param[[cc5]]$fixed_cpts==F})))){
    if(cpt_merge_option=="unanimous"){
      estimated_cpts=w_cps_new
    }else if(cpt_merge_option=="separate"){
      names(w_cp_param_new)=w_cp_param_variables
      estimated_cpts=w_cp_param_new
    }
  }else{
    estimated_cpts=NULL
  }
  # output result
  return(list(result=result,estimated_cpts=estimated_cpts,out_filter=out_filter,out_smooth=out_smooth,iter=iter))
} 