#include <Rcpp.h>
#include <vector>
using namespace Rcpp;
// this is a cpp file to calculate CI matrix

// [[Rcpp::export]]
NumericMatrix calculate_CI(NumericVector beta_hat, NumericVector se_beta_hat){
  int p = beta_hat.size();
  NumericMatrix CI(p,2);
  // Initialize a vector of vector or 2D vector of size 5X4 with value 10
  //vector < vector<int> > CI (p, vector<int>(2, -1) );
  for (int i =0 ; i < p; i++){
    CI[i] = beta_hat[i] - 1.96*se_beta_hat[i];
    //CI[i] = beta_hat[i] + 1.96*se_beta_hat[i];
  }
  for (int j =0 ; j < p; j++){
    CI[p+j] = beta_hat[j] + 1.96*se_beta_hat[j];
  }
  return(CI);
}
