#include <string>
#include <iostream>
#include <Rcpp.h>
#include <stdio.h>
#include <math.h>    

// [[Rcpp::plugins(cpp11)]]
using namespace std;
using namespace Rcpp;


/*define one function */
// [[Rcpp::export]]
NumericVector adjust_R_Rcpp_short_binomial(NumericMatrix mat, double prob) {  
NumericVector out;
out=mat(_,1);
NumericVector corrected_pv(out.size());
for(int j=0; j<out.size(); j++){ 
corrected_pv[j]=( Rf_choose( (min(mat(j,1),mat(j,2))),mat(j,0)) * pow (prob,mat(j,0)) * pow((1-prob),min(mat(j,1),mat(j,2))-mat(j,0))  );
}
 return corrected_pv;
}
