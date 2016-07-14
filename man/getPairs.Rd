\name{getPairs}
\alias{getPairs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pairs that collided at a certain date}
\description{
Pairs that shared at least a single visit at a certain date}
\usage{
getPairs(ids)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ids}{
%%     ~~Describe \code{ids} here~~
}
}
\details{
this function is called by get_observed_pairs which goes date by date for each unique date
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

# get unadjusted pairs
unadjusted_observed_pairs<-get_observed_pairs(db_dates)

}
