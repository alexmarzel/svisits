\name{adjust_R_Rcpp_short_binomial}
\alias{adjust_R_Rcpp_short_binomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probabilites of sharing a given number of visits
}
\description{
This is the C++ function that is called by \code{Bonferroni_m}.
}
\usage{
adjust_R_Rcpp_short_binomial(mat, prob)
}
\arguments{
  \item{mat}{
matrix with 3 columns (unadjusted number of shared visits, number of visits of the first member, number of visits of the second member)
}
  \item{prob}{
probability of sharing a single visit, default 1/75
}
}
\details{
written in C++
}
\value{
a vector of probabilities to share the given number of visits for each pair 
}
\examples{
visits <- data.frame(Freq = c(1,1,1,2,3),
                     N_visits.x = c(46,39,10,40,55),
                     N_visits.y = c(68,68,68,24,24))
mat <- as.matrix(visits)
Bonferroni_mat <- adjust_R_Rcpp_short_binomial(mat, prob=1/75)
Bonferroni_mat
}
