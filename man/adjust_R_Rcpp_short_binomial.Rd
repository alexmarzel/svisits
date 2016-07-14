\name{adjust_R_Rcpp_short_binomial}
\alias{adjust_R_Rcpp_short_binomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probabilites of sharing a given number of visits
}
\description{
this is the C++ function that is called by Bonferroni_m
}
\usage{
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
pairs as found by the function get_observed_pairs
}
  \item{prob}{
probability of sharing a single visit, default 1/75
}
}
\details{
written in C++
}
\value{
returns a matrix with a probability to share the given number of visits for each pair 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

# load data
data("simulated_data")
db_dates<-prepare_db(your_database = simulated_data,ids_column = "subject",dates_column = "sim_dates")

# first get unadjusted pairs
unadjusted_observed_pairs<-get_observed_pairs(db_dates)
# prepare ids
ids<- data.table(ids=as.character(names(table(db_dates$subject))), 
                 N_visits=as.character(as.numeric(table(db_dates$subject)))) 
setkey(ids, "ids")

# run
Bonferroni_m_output<-Bonferroni_m(unadjusted_observed_pairs,prob = 1/75,alpha =0.01)
length(Bonferroni_m_output[,1])
}
