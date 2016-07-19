\name{Shuffling_simulation}
\alias{Shuffling_simulation}

\title{
Shuffling simulation to determine shared visits false-positive threshold}
\description{
Shuffling simulation to determine shared visits false-positive threshold. It performs only one run. Please use \code{"replicate"} for several (see Example).
}
\usage{
Shuffling_simulation(db_dates, ids, adjusted_observed, max_threshold=NULL)
}

\arguments{
  \item{db_dates}{
a data frame with visits to be shuffled (in long format) which column "fupdatePeriod" containing the time period which each clinic visit belongs to.
}
 \item{ids}{a \code{"data.table"} object with columns "ids" and "N_visits". The "ids" column represents patient identifiers and it is the reference (see \code{setkey} from package \code{data.table} for more details). The "N_visits" contains the number of visits for each patient.}
  \item{adjusted_observed}{
a data frame with adjusted observed number of visits.
}
  \item{max_threshold}{
  a maximum false-positive threshold to be inspected. If it is \code{NULL}, the 99.995\%-quantile of the distribution of the unadjusted number of shared visits is used.
  }
}
\details{
Due to the multiple comparisons problem and strongly simplifying assumptions about the uniform distribution of visits, the corrected number of shared visits is not an optimal test statistic. Instead, the visits are first shuffled within each period (such that the original distribution of the number of visits per individual is preserved) and the number of randomly collided shared visits per pair is counted and penalized using "adjust_visits". In other words, the observed visit dates from a given period (by default a quater) are randomly re-assigned between the patients that attended the clinic during this period.

Define the number of desired shuffling simulations (ideally > 100, but it will take some time) and run.
}
\value{
a numeric vector of size \code{max_threshold}. Each value represents the false positive fraction for a cutoff of adjusted shared visits that corresponds to the position (index) of the value in the vector.
}
\examples{
# load data
data("simulated_data")
db_dates <- prepare_db(your_database = simulated_data,
                      ids_column = "subject",
                      dates_column = "sim_dates")

# prepare 3 months periods for later shuffling
db_dates <- cbind(db_dates,month.day.year(db_dates$dates))
periodMo <- 3
db_dates <- cbind(db_dates,
                  fupdatePeriod=with(db_dates, year+round(floor((month-1)/periodMo)
                                                           *periodMo/12,
                                                          digits=2)))

# first get unadjusted pairs
unadjusted_observed_pairs <- get_observed_pairs(db_dates)
ids <- data.table(ids=as.character(names(table(db_dates$subject))), 
                  N_visits=as.character(as.numeric(table(db_dates$subject)))) 
setkey(ids, "ids") 

# now adjust the pairs
adjusted_observed <- adjust_visits(unadjusted_pairs = unadjusted_observed_pairs,
                                  ids = ids,
                                  prob = 1/75)
# number of simulations
n_shuffl <- 2 # two 
Shuffling_simulation_output <- replicate(n_shuffl,
                                         Shuffling_simulation(db_dates,
                                                              ids=ids,
                                                              adjusted_observed),
                                         simplify=FALSE)
head(Shuffling_simulation_output)
}