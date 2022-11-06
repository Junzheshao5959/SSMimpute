#' Exploration function as initial guess
#'
#' @param res result subject, generated from \code{shrinkTVP} function
#' @param plotFlag whether we need to plot the result, default value is T
#'
#' @return Tables including non-stationary coefficient, possible change points, fitted ARIMA model
#' @export
#'
#' @examples
explore_SSM = function(res = res, plotFlag = T){
  if (plotFlag){
    plot(res)
    plot(res, pars = "theta_sr")
  }
  S1 = summary(res)
  #S1$summaries$beta_mean
  #cat("variables with non-zero theta are:","\n")
  theta_res = S1$summaries$theta_sr %>% as.data.frame() %>%  rownames_to_column( var = "variance") %>% as.tibble() %>% filter(HPD1 > 1e-3)
  vec_res = unlist(strsplit(theta_res$variance,  "_sr_"))[2*seq(1, length(theta_res$variance))]
  theta_res$variance = vec_res
  beta_res = S1$summaries$beta_mean %>% as.data.frame() %>%  rownames_to_column(var = "coef") %>% as.tibble()
  beta_res = beta_res %>% mutate(coef = unlist(strsplit(beta_res$coef,  "_mean_"))[2*seq(1, length(beta_res$coef))])
  print(beta_res)
  readline(prompt="Press [enter] to see the non-stationary coefficient: ")
  print(theta_res)
  id_list = match(vec_res,beta_res$coef)
  #summary(res)
  #summary(res,showprior = F)
  beta = res$beta


  cat("indicating non-stationary process for: ", vec_res,"; \n")
  readline(prompt="Press [enter] to further checking for change points: ")
  for ( i in id_list){
    b_1 = as.data.frame(beta[[i]])
    B_1 = b_1 %>% summarise_if(is.numeric, mean) %>% t() %>% ts()
    plot(bcp(B_1))
  }
  readline(prompt="Press [enter] to further checking fitted ARIMA model")
  for ( i in id_list){
    b_1 = as.data.frame(beta[[i]])
    B_1 = b_1 %>% summarise_if(is.numeric, mean) %>% t() %>% ts()

    res_arima = forecast::auto.arima(B_1)
    print(summary(res_arima))
  }
}
