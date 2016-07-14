\name{barcode}
\alias{barcode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Unique pair identifier}
\description{
Makes a unique identifier for each pair. Each member`s id should be numeric! For example  id_1=6758 , id_2=3456 the barcode will be 3456_6578}
\usage{
barcode(raw_for_eval)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raw_for_eval}{
%%     ~~Describe \code{raw_for_eval} here~~
}
}
\details{
}
\value{
string of a form id1_id2
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

barcode(c(id_2=6578,id_1=3456))

## The function is currently defined as
function (raw_for_eval) 
{
    return(as.character(paste0(min(as.numeric(raw_for_eval)), 
        "_", max(as.numeric(raw_for_eval)))))
  }
}
