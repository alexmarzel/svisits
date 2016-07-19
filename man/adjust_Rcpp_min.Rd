\name{adjust_Rcpp_min}
\alias{adjust_Rcpp_min}

\title{
Adjustment function
 }
\description{
This is the C++ adjustment function that is called by adjust_visits. It adjusts the numbers of shared visits for each pair.
}
\usage{
adjust_Rcpp_min(mat, prob)
}
\arguments{
  \item{mat}{
matrix with 3 columns (unadjusted number of shared visits, number of visits of the first member, number of visits of the second member) to adjust
}
  \item{prob}{
probability of sharing a single visit, default 1/75
}
}
\details{
written in C++
}
\value{
a vector of adjusted number of shared visits
}

\examples{
visits <- data.frame(Freq = c(1,1,1,2,3),
                     N_visits.x = c(46,39,10,40,55),
                     N_visits.y = c(68,68,68,24,24))
mat <- as.matrix(visits)
adj_mat <- adjust_Rcpp_min(mat, prob=1/75)
adj_mat
}