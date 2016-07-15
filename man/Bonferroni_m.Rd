\name{Bonferroni_m}
\alias{Bonferroni_m}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Detecting the pairs using Bonferroni correction 
 }
\description{
much faster than the shuffling
}
\usage{
Bonferroni_m(unadjusted_pairs, ids=ids,prob = 1/75, alpha = 0.01, only_significant = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{unadjusted_pairs}{
pairs as found by the function get_observed_pairs
}
 \item{ids}{
patient ids}
  \item{prob}{
probability of sharing a single visit, default 1/75
}
  \item{alpha}{
type 1 error
}
  \item{only_significant}{
return only significant pairs above the threshold (default) or all the pairs}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
Bonferroni_m_output<-Bonferroni_m(unadjusted_observed_pairs,ids=ids,prob = 1/75,alpha =0.01)
length(Bonferroni_m_output[,1])
}
