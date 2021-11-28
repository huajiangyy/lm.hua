#' A custom function that calculates multiple linear regress
#'
#' It returns a table of beta coefficient, standard error, t value, and p value
#' as variables, and predictors as observations,
#' @param predictor the predictor variables
#' @param outcome the outcome variables
#' @return returns beta, SE, t value, p value of the linear regression
#'
#' @examples
#' Y = c(5.6, 7.9, 10.8)
#' X = c(1, 2, 3)
#' linear_regression(X, Y)
#'
#' @export
#' linear_regression(predictor1, predictor2, ..., outcome)

linear_regression = function(predictor, outcome){
  #get the size of outcome
  n = length(predictor)

  #prepare predictor matrix
  temp = rep(1,n)
  predictor = cbind(temp,predictor)

  #get the number of predictor
  p = ncol(predictor)

  #calculate transpose of predictor matrix
  predictor_tran = t(predictor)
  #predictor_tran = rbind(temp,predictor_tran)

  #calculate beta vector
  beta_hat = (solve(predictor_tran %*% predictor)) %*% predictor_tran %*% outcome

  #calculate hat matrix
  hat_mat = predictor %*% (solve(predictor_tran %*% predictor)) %*% predictor_tran

  #calculate y_hat
  #outcome_hat = hat_mat %*% outcome
  outcome_hat = predictor %*% beta_hat

  #calculate residual
  e_hat = outcome - outcome_hat

  #calculate sigma^2
  #sigma_sq = sse/degree_freedom
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
  #print(output)

  #get outcome table
  outcome_table = cbind(outcome, outcome_hat, e_hat)
  colnames(outcome_table) = c('Yi', 'Yi_hat', 'ei')
  #print(outcome_table)

  #calculate SSE
  sse = t(e_hat) %*% e_hat
  #(sse)

  #calculate SSY
  ident_mat = diag(n)
  one_mat = matrix(1/n, nrow = n, ncol = n)
  ssy = t(outcome) %*% (ident_mat - one_mat) %*% outcome
  #print(ssy)

  #calculate SSR
  ssr = t(outcome) %*% (hat_mat - one_mat) %*% outcome
  #print(ssr)

  #calculate adjusted R-squared
  R_squared = ssr/ssy

  #calculate F statistics
  Fstat = (ssr/(p-1))/(sse/(n-p))

  #calculate confidence interval
  CI = matrix(NA, nrow = p, ncol = 2)
  for (i in c(1:p)){
    CI[i, 1] = beta_hat[i] - 1.96*se_beta_hat[i]
    CI[i, 2] = beta_hat[i] + 1.96*se_beta_hat[i]
  }
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
  uni_stat_outcome[i, 1] = n
  uni_stat_outcome[i, 2] = sum(outcome)/n
  uni_stat_outcome[i, 3] = sd(outcome)
  uni_stat_outcome[i, 4] = min(outcome)
  uni_stat = rbind(uni_stat_predictor, uni_stat_outcome)
  row_name = ""
  for (i in c(1:p-1)){
    temp_name = paste0("X", i)
    row_name[i] = temp_name
  }
  row_name = rbind(row_name, "Outcome")
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

#package.skeleton("mypackage", list = c("linear_regression"))
#this is a new function

