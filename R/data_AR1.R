#' Sample data for a complex structure, with  AR(1) \eqn{\beta_{C,t}}
#'
#'
#' @format A data frame with 1000 rows and 4 variables, the data generation process follows the following equation
#' \deqn{Y_t = 40+0.5Y_{t-1}-\beta_{X,t}X_t-0.5X_{t-1}-\beta_{C}C_t+v_t, v_t \sim \text{i.i.d }\mathcal{N}(0,1)\\ \text{where } \beta_{C,t} =0.5 \beta_{C,t-1}-1+w_{\beta,t}, w_{\beta,t}\sim \mathcal{N}(0,0.1)}
#' \describe{
#'   \item{Date}{A simulated date}
#'   \item{y}{outcome, with data generating process as above}
#'   \item{c}{continous covarite}
#'   \item{x}{continous treatment}
#'
#' }
"data_AR1"
