#' Sample data for a complex structure, with periodic stable \eqn{\beta_{X,t}} and random walk \eqn{\beta_{C,t}}
#'
#'
#' @format A data frame with 1000 rows and 4 variables, the data generation process follows the following equation
#' \deqn{Y_t = 40+0.5Y_{t-1}-\beta_{X,t}X_t-0.5X_{t-1}-\beta_{C}C_t+v_t, v_t \sim \text{i.i.d }\mathcal{N}(0,1)\\ \text{where } \beta_{X,t}= -I(t\in[0,400])-2.5I(t\in(400,700]-I(t\in(700,1000])), \\\beta_{C,t} =\beta_{C,t-1}+1.5+w_{\beta,t}, w_{\beta,t}\sim \mathcal{N}(0,0.05)}
#' \describe{
#'   \item{Date}{A simulated date}
#'   \item{y}{outcome, with data generating process as above}
#'   \item{y_1}{one step laggged out come}
#'   \item{c}{continous covarite}
#'   \item{x}{continous treatment}
#'   \item{x_1}{one step lagged continous treatment}
#'
#' }
"data_complex"
