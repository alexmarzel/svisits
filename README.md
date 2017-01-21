# R package: svisits
Using mathematical approximations and a shuffling approach to derive the probabilities of sharing a given number of 
visits by chance the package helps to identify couples that may represent either transmission pairs or serosorting couples 
in a stable relationship, for biological and epidemiologcial HIV research.

# Installation
Make sure the most updated version of R is installed.

The packages "compiler" and "parallel" should be included by default. If not, update R. 

## Install dependencies
```r
install.packages(c("Rcpp", "chron","data.table")) 

```
## Install the "svisits" package
```r
install.packages("svisits_1.2.0.tar.gz",
                 repos = NULL, type = "source")  
```              
Note: Place "svisits_1.2.0.tar.gz" into the working directory or provide the full path.

# Fastest use, an example

```r
library(svisits) 

data("simulated_data") 

pairs <- find_pairs(simulated_data, 
                    ids_column = "subject",
                    dates_column = "sim_dates")
                    
head(pairs$Bonferroni_pairs)
```   
