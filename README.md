## svisits
Using mathematical approximations and a shuffling approach to derive the probabilities of sharing a given number of 
visits by chance the package helps to identify couples that may represent either transmission pairs or serosorting couples 
in a stable relationship, for biological and epidemiologcial HIV research.

## Installation
Make sure you the most updated version of R is installed 
The packages "compiler" and "parallel" should be included by default.

## Install dependencies
install.packages(c("Rcpp", "chron","data.table")) 
install.packages("svisits_1.2.0.tar.gz",
                 repos = NULL, type = "source")  # Install the "svisits" package

## Fastest use
library(svisits) 
data("simulated_data") 
pairs <- find_pairs(simulated_data, 
                    ids_column = "subject",
                    dates_column = "sim_dates")
head(pairs$Bonferroni_pairs)
