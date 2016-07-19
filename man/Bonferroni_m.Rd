\name{Bonferroni_m}
\alias{Bonferroni_m}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Detecting the pairs using Bonferroni correction 
 }
\description{
This method is much faster than the shuffling method, but less precise.
}
\usage{
Bonferroni_m(unadjusted_pairs, 
             ids=ids, 
             prob = 1/75, 
             alpha = 0.01, 
             only_significant = TRUE)
}
\arguments{
  \item{unadjusted_pairs}{
a \code{"table"} object containing the pairs and the number of shared visits. See \code{\link{get_observed_pairs}}.
}
 \item{ids}{a \code{"data.table"} object with columns "ids" and "N_visits". The "ids" column represents patient identifiers and is the reference (see \code{setkey} from package \code{data.table} for more details). The "N_visits" contains the number of visits for each patient.}
  \item{prob}{
probability of sharing a single visit. Default \code{1/75}.
}
  \item{alpha}{
type I error
}
  \item{only_significant}{
returns only significant pairs above the threshold (default) or all the pairs}
}
\details{
}
\value{a data frame with
\describe{
  \item{\code{"allPairs"}}{pair identifier}
  \item{\code{"Freq"}}{the unadjusted number of shared visits}
  \item{\code{"id_1"}}{identifier of the first pair member}
  \item{\code{"id_2"}}{identifier of the second pair member}
  \item{\code{"N_visits.x"}}{total number of visits of the first pair member}
  \item{\code{"N_visits.y"}}{total number of visits of the second pair member}
  \item{\code{"Prob_for_Bonferr"}}{probability of sharing the given number of shared visits}
  \item{\code{"BP"}}{row number}
  \item{\code{"ltp"}}{the log-transformed probablity of sharing the given number shared visits}
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
# prepare ids
ids <- data.table(ids = as.character(names(table(db_dates$subject))), 
                  N_visits = as.character(as.numeric(table(db_dates$subject)))) 
setkey(ids, "ids")

# run
Bonferroni_m_output <- Bonferroni_m(unadjusted_observed_pairs,
                                    ids = ids, prob = 1/75, alpha = 0.01)
# number of significant pairs
nrow(Bonferroni_m_output)
}
