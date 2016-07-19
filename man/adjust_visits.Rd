\name{adjust_visits}
\alias{adjust_visits}

\title{
Adjust the number of shared visits}
\description{
The function adjusts the total number of shared visits per pair.
}
\usage{
adjust_visits(unadjusted_pairs, ids, prob = 1/75)
}
\arguments{
  \item{unadjusted_pairs}{a \code{"table"} object containing the pairs and the number of shared visits
}
  \item{ids}{a \code{"data.table"} object containing columns "ids" and "N_visits". The "ids" column represents patient identifiers and is the reference (see \code{setkey} from package \code{data.table} for more details). The "N_visits" contains the number of visits for each patient.}
  \item{prob}{probability of sharing a single visit, depending on the frequency of visits. Default \code{1/75}~visit every three months.
}
}
\details{
}
\value{
a data frame with the columns:
\describe{
  \item{\code{"allPairs"}}{pair identifiers}
  \item{\code{"Freq"}}{unadjusted number of shared visits}
  \item{\code{"id_1"}}{first member of the pair's id}
  \item{\code{"id_2"}}{second member of the pair's id}
  \item{\code{"N_visits.x"}}{first member of the pair's number of visits}
  \item{\code{"N_visits.y"}}{second member of the pair's number of visits}
  \item{\code{"Corrected"}}{the adjusted number of shared visits per pair}
}
}
\examples{
# load data
data("simulated_data")
db_dates <- prepare_db(your_database = simulated_data,
                       ids_column = "subject",
                       dates_column = "sim_dates")

# first get unadjusted pairs
unadjusted_observed_pairs <- get_observed_pairs(db_dates)
ids <- data.table(ids = as.character(names(table(db_dates$subject))), 
                  N_visits = as.character(as.numeric(table(db_dates$subject)))) 
setkey(ids, "ids") 

# now adjust the pairs
adjusted_observed <- adjust_visits(unadjusted_pairs = unadjusted_observed_pairs,
                                   ids = ids, 
                                   prob = 1/75)
head(adjusted_observed)
}
