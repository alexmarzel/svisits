\name{simulated_data}
\alias{simulated_data}
\docType{data}
\title{
Simulated dataset of visit dates based on the SHCS
}
\description{
It contains 20 pseudo transmission pairs that should be detected with the method.
}
\usage{data("simulated_data")}
\format{
  A data frame with 90266 observations on the following 3 variables:
  \describe{
    \item{\code{X}}{row's indexing}
    \item{\code{subject}}{a numeric vector with patient identifiers}
    \item{\code{sim_dates}}{a factor with the visit dates}
  }
}

\examples{
data("simulated_data")
head(simulated_data)
}
