\name{Shuffling_simulation}
\alias{Shuffling_simulation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Shuffling simulation to determine shared visits false-positive threshold}
\description{
Shuffling simulation to determine shared visits false-positive threshold, one run, use "replicate" for several, see example
}
\usage{
Shuffling_simulation(db_dates, adjusted_observed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{db_dates}{
visits to be shuffled (long format) an output from prepare_db, see data("simulated_data)
}
  \item{adjusted_observed}{
The dataset with adjusted observed number of visits
}
}
\details{
Due to the multiple comparisons problem and a strongly simplifying assumptions about the uniform distribution of visits, the corrected number of shared visits is not an optimal test statistic. Instead, the visits are first shuffled within each quarter (such that the original distribution of the number of visits per individual was preserved) and the number of randomly collided shared visits per pair is counted and penalized using "adjust_visits". In other words, the observed visit dates from a given quarter were randomly re-assigned between the patients that attended during this quarter. 

Define the number of desired shuffling simulations (ideally >100, but it will take some time) and run

}
\value{
returns a numeric vector of size 10 with each value representing 
the false positive fraction for a cutoff of adjusted shared visits that corresponds 
to the position of the value in the order of the vector. 
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

# prepare 3 months periods for later shuffling
db_dates<-cbind(db_dates,month.day.year(db_dates$dates));periodMonths<-3
db_dates<-cbind(db_dates,fupdatePeriod=db_dates$year
                +round(floor((db_dates$month-1)/periodMonths) *periodMonths/12,2))

# first get unadjusted pairs
unadjusted_observed_pairs<-get_observed_pairs(db_dates)
ids<- data.table(ids=as.character(names(table(db_dates$subject))), 
                 N_visits=as.character(as.numeric(table(db_dates$subject)))) 
setkey(ids, "ids") 

# now adjust the pairs
adjusted_observed<-adjust_visits(unadjusted_pairs=unadjusted_observed_pairs,prob = 1/75)
# number of simulations
n_shuffl<-2 # two 
Shuffling_simulation_output<- replicate(n_shuffl,
                                        Shuffling_simulation(db_dates,adjusted_observed),
                                        simplify=FALSE)

head(Shuffling_simulation_output)
}