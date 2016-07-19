\name{barcode}
\alias{barcode}
\title{
Unique pair identifier}
\description{
It makes a unique identifier for each pair. Each member`s id should be numeric! For example for \code{id_1=6758} and \code{id_2=3456} the barcode will be \code{"3456_6578"}.}
\usage{
barcode(raw_for_eval)
}
\arguments{
  \item{raw_for_eval}{a vector with two ids
}
}
\details{
}
\value{
string of a form id1_id2
}

\examples{
barcode(c(id_2=6578,id_1=3456))
barcode(c(3456, 6578))
}
