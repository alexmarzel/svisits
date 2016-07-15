\name{adjust_visits}
\alias{adjust_visits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adjust the number of shared visits }
\description{
Adjusts the total number of shared visits per pair
}
\usage{
adjust_visits(unadjusted_pairs,ids=ids, prob = 1/75)
}
\arguments{
  \item{unadjusted_pairs}{
The pairs
}
 \item{ids}{
patient ids}
  \item{prob}{
probability of sharing a single visit, depends on the frequency of visits. Default 1/75
~visit every three months.
}
}
\details{
}
\value{
returns data.frame with the column "Corrected" representing the adjsuted number 
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
adjusted_observed<-adjust_visits(unadjusted_pairs=unadjusted_observed_pairs,ids=ids, prob = 1/75)
head(adjusted_observed)
}
