\name{getPairs}
\alias{getPairs}

\title{
Pairs with shared visit on a certain date}
\description{
From a given list of patient with the clinic visit on a certain date the function returns all the pairs.}
\usage{
getPairs(ids)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ids}{
a vector of patient identifiers
}
}
\details{
This function is called by \code{\link{get_observed_pairs}} which date-wise finds all the pairs with at least one shared visit (for all the dates that appear in your database).
}
\value{
a character vector containing the \code{\link{barcode}}s of the found pairs on a certain date
}
\examples{
dates_long <- data.frame(subject=c(12,15,2,202,31,3), 
                         dates=as.Date(c(rep("2001-03-16", 5), "2011-03-17")))
# subjects with a visit on "2001-03-16"
id_s <- dates_long$subject[which(dates_long$dates==as.Date("2001-03-16"))]
getPairs(id_s)
}
