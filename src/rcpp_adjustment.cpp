#include <string>
#include <iostream>
#include <Rcpp.h>
#include <stdio.h>
#include <math.h>    

// [[Rcpp::plugins(cpp11)]]
using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
double chooseC(double n, double k) {
  return Rf_choose(n, k);
}

/*define one function */
// [[Rcpp::export]]
NumericVector adjust_Rcpp_min(NumericMatrix mat, double prob) {  
NumericVector out;
out=mat(_,1);
NumericVector corrected_S(out.size());
for(int j=0; j<out.size(); j++){ 
corrected_S[j]=mat(j,0)-(log(Rf_choose( min(mat(j,1),mat(j,2)),mat(j,0))) /-log(prob) );
}
 return corrected_S;
}



