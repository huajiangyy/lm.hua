#' A custom function that calculates multiple linear regression
#'
#' It returns a table of beta coefficients, standard error, t value, and p value
#' as variables, and predictors as observations. It also prints out confidence
#' intervals and descriptive univariate statistics.
#'
#' @param predictor the predictor variables, can consist of multiple variables.
#' It has to be in nxm matrix format
#'
#' @param outcome the outcome variables. It has to be in nx1 matrix format
#'
#' @return returns a table of beta, SE, t value, p value of the linear regression
#'
#' A sample health dataset is included and loaded automatically when loading the
#' package
#'
#' @examples
#' Y = matrix(c(5.6, 7.9, 10.8), ncol=1)
#' X = matrix(c(1, 2, 3), ncol=1)
#' linear_regression(X, Y)
#'
#' Y = matrix(c(1714,1664,1760,1685,1693,1670,1764,1764,1792,1850,1735,1775),
#' ncol=1)
#' X = matrix(c(2.4,2.52,2.54,2.74,2.83,2.91,3,3,3.01,3.01,3.02,3.07), ncol=1)
#' linear_regression(X, Y)
#'
#' Y = mydata$death_rate
#' X = cbind(mydata$doctor_num, mydata$hos_num, mydata$capita_nuuanl_income,
#' mydata$population_den)
#' linear_regression(X, Y)
#'
#' @export
#' linear_regression(predictor, outcome)

linear_regression = function(predictor, outcome){
  #get the size of outcome
  n = nrow(predictor)

  #prepare predictor matrix
  temp = rep(1,n)
  predictor = cbind(temp,predictor)

  #get the number of predictor
  p = ncol(predictor)

  #calculate transpose of predictor matrix
  predictor_tran = t(predictor)

  #calculate beta vector
  beta_hat = (solve(predictor_tran %*% predictor)) %*% predictor_tran %*% outcome

  #calculate hat matrix
  hat_mat = predictor %*% (solve(predictor_tran %*% predictor)) %*% predictor_tran

  #calculate y_hat
  outcome_hat = predictor %*% beta_hat

  #calculate residual
  e_hat = outcome - outcome_hat

  #calculate sigma^2
  sigma_sq = t(e_hat) %*% e_hat / (n - p)

  #Variance of betahat
  var_beta_hat = diag( solve(t(predictor)%*%predictor))*c(sigma_sq)
  se_beta_hat = sqrt(var_beta_hat)

  #Inference: t statistics and p val for H0: beta = 0
  t_statistic = c(beta_hat/se_beta_hat)
  P_value = c(2* (1-pt(q=abs(t_statistic), df=n-p)))

  output = cbind(Estimate=c(beta_hat),
                 Std_Err=se_beta_hat,
                 t_statistic=t_statistic,
                 p_value=P_value)

  beta_name = rep("Intercept", p)
  for (i in c(2:p)){
    temp_name = paste0("X", i-1)
    beta_name[i] = temp_name
  }
  rownames(output) = beta_name

  #get outcome table
  outcome_table = cbind(outcome, outcome_hat, e_hat)
  colnames(outcome_table) = c('Yi', 'Yi_hat', 'ei')

  #calculate SSE
  sse = t(e_hat) %*% e_hat

  #calculate SSY
  ident_mat = diag(n)
  one_mat = matrix(1/n, nrow = n, ncol = n)
  ssy = t(outcome) %*% (ident_mat - one_mat) %*% outcome

  #calculate SSR
  ssr = t(outcome) %*% (hat_mat - one_mat) %*% outcome

  #calculate adjusted R-squared
  R_squared = ssr/ssy

  #calculate F statistics
  Fstat = (ssr/(p-1))/(sse/(n-p))

  #calculate confidence interval, using Rcpp
  library(Rcpp)
  src = "
  #include <Rcpp.h>
  #include <vector>
  using namespace Rcpp;
  //' this function calcualte the confidence interval
  //'
  //' @param beta_hat beta hat matrix
  //' @param se_beta_hat standar error matrix
  //' @return confidence interval

  // [[Rcpp::export]]
  NumericMatrix calculate_CI(NumericVector beta_hat, NumericVector se_beta_hat){
    int p = beta_hat.size();
    NumericMatrix CI(p,2);
    // Initialize a vector of vector or 2D vector of size 5X4 with value 10
    //vector < vector<int> > CI (p, vector<int>(2, -1) );
    for (int j =0 ; j < p; j++){
      CI[j] = beta_hat[j] - 1.96*se_beta_hat[j];
    }
    for (int k =0 ; k < p; k++){
      CI[p+k] = beta_hat[k] + 1.96*se_beta_hat[k];
    }
    return(CI);
  }"
  #sourceCpp("src/code.cpp")
  sourceCpp(code = src)
  CI = calculate_CI(beta_hat, se_beta_hat)
  rownames(CI) = beta_name
  colnames(CI) = c('2.5 %', '97.5 %')

  #calculate descriptive univariate statistics
  uni_stat_predictor = matrix(NA, nrow = p-1, ncol = 4)
  for (i in c(1:p-1)){
    uni_stat_predictor[i, 1] = n
    uni_stat_predictor[i, 2] = sum(predictor[, i+1])/n
    uni_stat_predictor[i, 3] = sd(predictor[, i+1])
    uni_stat_predictor[i, 4] = min(predictor[, i+1])
  }
  uni_stat_outcome = matrix(NA, nrow = 1, ncol = 4)
  uni_stat_outcome[1, 1] = n
  uni_stat_outcome[1, 2] = sum(outcome)/n
  uni_stat_outcome[1, 3] = sd(outcome)
  uni_stat_outcome[1, 4] = min(outcome)
  uni_stat = rbind(uni_stat_predictor, uni_stat_outcome)
  row_name = c()
  for (i in c(1:p-1)){
    temp_name = paste0("X", i)
    row_name[i] = temp_name
  }
  row_name = append(row_name, "Outcome")
  rownames(uni_stat) = row_name
  colnames(uni_stat) = c('N', 'Mean', 'StdDev', 'Minimum')

  print(uni_stat)
  print(CI)
  print(paste0("R-squred: ", R_squared))
  a = paste0("F-statistic: ", Fstat)
  b = paste0(" on ", (p-1))
  c = paste0(" and ", (n-p))
  d = paste0(" DF ")
  print(paste0(a,b,c,d))

  return(output)
}
