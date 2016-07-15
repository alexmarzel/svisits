\name{adjust_Rcpp_min}
\alias{adjust_Rcpp_min}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
adjustment function
 }
\description{
This is the C++ adjustment function that is called by adjust_visits
}
\usage{
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
matrix to adjust
}
  \item{prob}{
probability of sharing a single visit, default 1/75
}
}
\details{
written in C++
}
\value{
returns the adjusted matrix, the column "Corrected" is the adjusted number
of shared visits per pair
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
ids<- data.table(ids=as.character(names(table(db_dates$subject))), 
                 N_visits=as.character(as.numeric(table(db_dates$subject)))) 
setkey(ids, "ids") 

# now adjust the pairs
adjusted_observed<-adjust_visits(unadjusted_pairs=unadjusted_observed_pairs,ids=ids,prob = 1/75)
head(adjusted_observed)
}