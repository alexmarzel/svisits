\name{find_pairs}
\alias{find_pairs}

\title{
Find significant pairs
}
\description{
Find significant pairs using the Bonferroni method and (optionally) with the shuffling algorithm to avoid False-Positives. This is the wrap-up of the following functions:
\itemize{
  \item{\code{\link{prepare_db}}}
  \item{\code{\link{get_observed_pairs}}}
  \item{\code{\link{Bonferroni_m}}}
  \item{\code{\link{adjust_visits}}} (only if the shuffling algorithm is performed)
  \item{\code{\link{Shuffling_simulation}}} (optional)
}
}
\usage{
find_pairs(your_database,
           ids_column,
           dates_column,
           prob = 1/75,
           alpha = 0.01,
           only_significant = TRUE,
           shuffling = FALSE,
           n_shuffl = 100,
           period_months = 3,
           FP_threshold = max,
           max_threshold = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{your_database}{a data frame with the clinic visit dates in a long format.}
  \item{ids_column}{character, column name with the person identifiers}
  \item{dates_column}{character, column name with the visit dates}
  \item{prob}{probability of sharing a single visit, default \code{1/75}}
  \item{alpha}{type I error}
  \item{only_significant}{returns only the significant pairs above the threshold (default) or all the pairs with the corresponding p-values}
  \item{shuffling}{logical indicating if the shuffling should be performed. By default (\code{FALSE}) only the Bonferroni method is used.}
  \item{n_shuffl}{number of repetitions of shuffling algorithm (default 100)}
  \item{period_months}{the length of the follow-up period in months (by default \code{3})}
  \item{FP_threshold}{the function used to get the False-Positive fraction from the shuffling repetitions at each position. By default the \code{max} is called at each threshold. Another possibilities are \code{mean} or \code{median}.}
  \item{max_threshold}{maximum false-positive threshold to be inspected. If it is \code{NULL}, the 99.995\%-quantile of the distribution of the unadjusted number of shared visits is used.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list containing
\describe{
  \item{\code{"unadjusted_pairs"}}{pairs with unadjusted number of shared visits}
  \item{\code{"Bonferroni_pairs"}}{pairs obtained with the Bonferroni correction. If \code{significant_only==TRUE} only those above the threshold are returned.}
  \item{\code{"adjusted_pairs"}}{pairs with adjusted number of shared visits returned only if \code{shuffling==TRUE}.}
  \item{\code{"Shuffled_pairs"}}{pairs obtained by the shuffling algorithm (returned only if \code{shuffling==TRUE}).}
  \item{\code{"shuffling_threshold"}}{threshold determined by the shuffling algorithm}
  \item{\code{"shuffling_simulation_output"}}{output from the shuffling function if \code{shuffling==TRUE}}
}
}
\seealso{
\code{\link{prepare_db}}, \code{\link{get_observed_pairs}}, \code{\link{Bonferroni_m}}, \code{\link{adjust_visits}}, \code{\link{Shuffling_simulation}}
}
\examples{
data("simulated_data")
pairs <- find_pairs(simulated_data, 
                    ids_column = "subject",
                    dates_column = "sim_dates",
                    shuffling = TRUE, n_shuffl = 2)
}