#' A custom function that calculates multiple linear regress
#'
#' It returns a table of beta coefficient, standard error, t value, and p value
#' as variables, and predictors as observations,
#' @param predictor the predictor variables
#' @param target the target variables
#' @return returns beta, SE, t value, p value of the linear regression
#'
#' @examples
#' Y = c(5.6, 7.9, 10.8)
#' X = c(1, 2, 3)
#' linear_regression(X, Y)
#'
#' @export
#' linear_regression(predictor1, predictor2, ..., target)

linear_regression = function(predictor, target){
  #get the size of target
  sample_size = length(target)

  #prepare predictor matrix
  temp = rep(1,sample_size)
  predictor = cbind(temp,predictor)

  #get the number of predictor
  n_predictor = ncol(predictor)

  #calculate transpose of predictor matrix
  predictor_tran = t(predictor)
  #predictor_tran = rbind(temp,predictor_tran)

  #calculate beta vector
  beta = (solve(predictor_tran %*% predictor)) %*% predictor_tran %*% target
  beta_name = rep(NA, n_predictor)
  for (i in c(1:n_predictor)){
    temp_name = paste0("beta", i-1)
    beta_name[i] = temp_name
  }
  rownames(beta) = beta_name
  colnames(beta) = c("Coefficient")

  #calcuate SEs
  se = rep(NA, n_predictor)
  for (i in c(1:n_predictor)){
   predictor_i = predictor[,i]
   predictor_i_mean = rep(sum(predictor_i)/sample_size,sample_size)
   diff_matrix = predictor_i-predictor_i_mean
   sum = 0
   for (j in c(1:sample_size)){
     sum = sum + (abs(diff_matrix[i]))^2
   }
   se[i] =  sum/sample_size/sqrt(sample_size)
  }
  beta = cbind(beta,se)

  #calculate t values
  t_value = rep(NA, n_predictor)
  for (i in c(1:n_predictor)){
    t_value[i] =  beta[,1][i]/beta[,2][i]
  }
  beta = cbind(beta,t_value)

  #calculate degree of freedom
  degree_freedom = sample_size - n_predictor

  #calculate p values
  p_value = rep(NA, n_predictor)
  for (i in c(1:n_predictor)){
    beta[,1][i]/se
    p_value[i] =  2*pt(abs(beta[,3][i]), degree_freedom, lower.tail = F)
  }
  beta = cbind(beta,p_value)

  #calculate hat matrix
  hat_mat = predictor %*% (solve(predictor_tran %*% predictor)) %*% predictor_tran

  #calculate y_hat
  target_hat = hat_mat %*% target

  #calculate residual
  e = target - target_hat

  #get target table
  target_table = cbind(target, target_hat, e)
  colnames(target_table) = c('yi', 'yi_hat', 'ei')

  #calculate SSres
  ss_res = t(e) %*% e

  #calculate sigma^2
  sigma_sq = ss_res/degree_freedom

  #print and return
  #print(beta)
  return(beta)
}

#package.skeleton("mypackage", list = c("linear_regression"))
#this is a new function

