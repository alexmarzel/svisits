\name{get_observed_pairs}
\alias{get_observed_pairs}
\title{
Pairs with at least one shared visit}
\description{
FInd all the pairs that shared at least a single visit.
}
\usage{
get_observed_pairs(tested_db)
}
\arguments{
  \item{tested_db}{
visits in long format with each row representing a visit date. See \code{\link{simulated_data}} for an example. It should contain columns \code{"subject"} with (numeric) patient identifiers and \code{"dates"} with the dates of clinic visits.
}
}
\details{
}
\value{
a table object containing the number of shared visits for all pairs with at least one shared visit 
}
\examples{
# a short example (with artificial data)
dates_long <- data.frame(subject=c(12,15,2,202,31,15, 202), 
                         dates=as.factor(c(rep("2001-03-16",5),rep("2011-06-20",2))))
get_observed_pairs(dates_long)

# example with the simulated dates based on the real distribution                         
# load data
data("simulated_data")
db_dates <- prepare_db(your_database = simulated_data,
                       ids_column = "subject",
                       dates_column = "sim_dates")

# get unadjusted pairs
unadjusted_observed_pairs <- get_observed_pairs(db_dates)
}

