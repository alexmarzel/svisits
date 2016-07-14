\name{prepare_db}
\alias{prepare_db}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare the database
}
\description{
Prepare the database for finding pairs.
}
\usage{
prepare_db(your_database = simulated_data,ids_column = "subject",dates_column = "sim_dates")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{your_database}{
the data.frame to be analyzed
}
\item{ids_column}{
column with the person identifiers 
}
\item{dates_column}{
column with the visit dates
}
}
\details{
the dates are ordered in an increasing order within each subject
}
\value{
returns the database with two columns "subject" and "dates"
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

data("simulated_data")
db_dates<-prepare_db(your_database = simulated_data,ids_column = "subject",dates_column = "sim_dates")
}
