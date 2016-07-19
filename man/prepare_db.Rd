\name{prepare_db}
\alias{prepare_db}

\title{
Prepare the database
}
\description{
Transform the database to the appropriate format for finding pairs.
}
\usage{
prepare_db(your_database, ids_column, dates_column)
}
\arguments{
  \item{your_database}{
the data frame to be analyzed
}
\item{ids_column}{
character, column name with the person identifiers 
}
\item{dates_column}{
character, column name with the visit dates
}
}
\details{
The dates are ordered in an increasing order within each subject.
}
\value{
a data frame with two columns:
\describe{
  \item{\code{"subject"}}{person identifiers}
  \item{\code{"dates"}}{dates of the clinic visits}
}
}
\examples{
data("simulated_data")
db_dates<-prepare_db(your_database = simulated_data,
                     ids_column = "subject",
                     dates_column = "sim_dates")
}
