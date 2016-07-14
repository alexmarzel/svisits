\name{get_observed_pairs}
\alias{get_observed_pairs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find all the pairs that shared at least a single visit }
\description{
Find all the pairs that shared at least a single visit.
}
\usage{
get_observed_pairs(tested_db)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tested_db}{
visits in a long format, each row a visit date, see data("simulated_data") for an example
}
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

# get unadjusted pairs
unadjusted_observed_pairs<-get_observed_pairs(db_dates)

}

