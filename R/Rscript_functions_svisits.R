
# Load packages
library("compiler");library("parallel"); 
library("chron"); library("data.table"); library("Rcpp");
enableJIT(3);


# prepare the database
prepare_db<-function(your_database, ids_column,dates_column)
{
  colnames(your_database)[which(colnames(your_database) %in% ids_column)]<-"subject"
  colnames(your_database)[which(colnames(your_database) %in% dates_column)]<-"dates"
  # order
  your_database[,c("subject","dates")]<-your_database[ order(as.numeric(your_database$subject),
                                                             as.numeric(your_database$dates)), c("subject","dates") ]
  return(your_database[,c("subject","dates")])
}


# Function that makes a unique identifier for each pair (ids have to be numeric)
barcode<-function(raw_for_eval){
  return(as.character(paste0(min(as.numeric(raw_for_eval)),"_",max(as.numeric(raw_for_eval)))))
}; barcode<-cmpfun(barcode)

# Corrected without duplicates
getPairs<-function(ids){
  if(length(ids)<2) return(c())
  comb<-expand.grid(ids,ids); comb<-comb[ comb$Var1!= comb$Var2,]
  comb<-cbind(comb,barcode=apply(X=comb[,1:2],1,barcode))
  return(unique(as.character(comb$barcode)))
}: getPairs<-cmpfun(getPairs)

# Real dataset for each center - unadjusted shared visits
get_observed_pairs<-function(tested_db){
  #Unique dates and periods
  datesUnique=unique(tested_db$dates); 
  # apply the above function to all unique dates ---> this gives all combination of ids that visited on the same date
  allPairs<-unlist(mclapply(datesUnique,mc.cores=detectCores(),FUN = function(x) getPairs(tested_db$subject[which(tested_db$dates==x)] )))
  return(table(allPairs))
};get_observed_pairs<-cmpfun(get_observed_pairs)

# make the Adjustemnt for the overall number of visits per pair

adjust_visits<-function(unadjusted_pairs, ids, prob=1/75){
  real_db = data.table::data.table(allPairs= unlist(names(unadjusted_pairs)), Freq=unlist(as.numeric(unadjusted_pairs)))
  real_db_1<-real_db[real_db$Freq >=1, ]
  real_db_1<-cbind(real_db_1,do.call("rbind",strsplit(as.character(real_db_1$allPairs),"_")))
  setnames(real_db_1,3:4,c("id_1","id_2")); setkey(real_db_1, "id_1");
  real_db_1<-na.omit(real_db_1[ids,allow.cartesian=TRUE]); setkey(real_db_1, "id_2");
  real_db_1<-na.omit(real_db_1[ids,allow.cartesian=TRUE]); setnames(real_db_1,5:6,c("N_visits.x","N_visits.y"))
  real_db_1$Freq<-as.numeric(real_db_1$Freq); real_db_1$N_visits.x<-as.numeric(real_db_1$N_visits.x);real_db_1$N_visits.y<-as.numeric(real_db_1$N_visits.y)  
  real_db_1<-cbind(real_db_1, Corrected=svisits::adjust_Rcpp_min(mat = as.matrix( real_db_1[,c(2,5,6),with = FALSE] ), prob=prob ))
  real_db_1<-as.data.frame(real_db_1)
  return(real_db_1)}
adjust_visits<-compiler::cmpfun(adjust_visits)

# Shuffling simulatoins
Shuffling_simulation<-function(db_dates,ids,adjusted_observed, max_threshold=NULL){
  tested_db_random<-db_dates
  #Unique dates and periods
  datesUnique=unique(db_dates$dates); fupdatePeriodUnique<-unique(db_dates$fupdatePeriod)
  # Shuffle within each 3 months period
  for(fs in fupdatePeriodUnique){
    tested_db_random[tested_db_random$fupdatePeriod==fs, "dates"]<-sample(tested_db_random[tested_db_random$fupdatePeriod==fs, "dates"])
  }
  # Get the random collisions
  RandomPairs<-unlist(mclapply(datesUnique,function(x) getPairs( tested_db_random$subject[which(tested_db_random$dates==x)] ),mc.cores = detectCores())) 
  # Adjust the randomized
  Adjusted_randomized<-adjust_visits(unadjusted_pairs =table(RandomPairs),ids=ids)
  # return vector of false-positive thresholds
  if (is.null(max_threshold)) {
    max_threshold <- ceiling(quantile(adjusted_observed$Freq, probs = 0.99995))
  }
  return(sapply(1:max_threshold,
                function(x) sum(Adjusted_randomized$Corrected>x)/sum(adjusted_observed$Corrected>x) )) 
}; Shuffling_simulation <-cmpfun(Shuffling_simulation)

# Bonfrerroni method
Bonferroni_m<-function(unadjusted_pairs,ids=ids,prob=1/75,alpha=0.01,only_significant=TRUE){
  real_db = data.table(allPairs= unlist(names(unadjusted_pairs)), Freq=unlist(as.numeric(unadjusted_pairs)))
  real_db_1<-real_db[real_db$Freq >=1, ]
  real_db_1<-cbind(real_db_1,do.call("rbind",strsplit(as.character(real_db_1$allPairs),"_")))
  setnames(real_db_1,3:4,c("id_1","id_2")); setkey(real_db_1, "id_1");
  real_db_1<-na.omit(real_db_1[ids,allow.cartesian=TRUE]); setkey(real_db_1, "id_2");
  real_db_1<-na.omit(real_db_1[ids,allow.cartesian=TRUE]); setnames(real_db_1,5:6,c("N_visits.x","N_visits.y"))
  real_db_1$Freq<-as.numeric(real_db_1$Freq); real_db_1$N_visits.x<-as.numeric(real_db_1$N_visits.x);real_db_1$N_visits.y<-as.numeric(real_db_1$N_visits.y)  
  real_db_1<-cbind(real_db_1, Prob_for_Bonferr=svisits::adjust_R_Rcpp_short_binomial(mat = as.matrix( real_db_1[,c(2,5,6),with = FALSE] ), prob=prob ))
  real_db_1<-as.data.frame(real_db_1)
  
  real_db_1<-cbind(real_db_1,BP=1:length(real_db_1[,1]))
  real_db_1<-cbind(real_db_1, ltp=-log10(real_db_1$Prob_for_Bonferr))
  
  selected<-real_db_1[ real_db_1$ltp>=-log10(alpha/length(real_db_1[,1])), ]
  
  if(only_significant==TRUE) {return(selected)} else {return(real_db_1)}
}
Shuffling_simulation <- compiler::cmpfun(Shuffling_simulation)

# wrap-up function
find_pairs <- function(your_database,
                       ids_column,
                       dates_column,
                       prob = 1/75,
                       alpha = 0.01,
                       only_significant = TRUE,
                       shuffling = FALSE,
                       n_shuffl = 100,
                       period_months = 3,
                       FP_threshold = max,
                       max_threshold=NULL) {
  
  # preprare the database
  db_dates <- prepare_db(your_database = your_database,
                         ids_column = ids_column,
                         dates_column = dates_column)
  
  # if shuffling is performed, add the periods
  if (shuffling) {
    db_dates <- cbind(db_dates, month.day.year(db_dates$dates))
    db_dates <- cbind(db_dates, fupdatePeriod = db_dates$year + round(floor((db_dates$month-1)/period_months)*period_months/12, 
                                                                      digits=2))
  }
  
  # unadjusted observed pairs:
  unadjusted_observed_pairs <- get_observed_pairs(db_dates)
  
  # construct ids:
  table_ids <- table(db_dates$subject)
  ids <- data.table(ids=as.character(names(table_ids)),
                    N_visits=as.character(as.numeric(table_ids)))
  setkey(ids, "ids")
  
  # perform the Bonferroni
  Bonferroni_m_output <- Bonferroni_m(unadjusted_observed_pairs,
                                      ids = ids,
                                      prob = prob,
                                      alpha = alpha,
                                      only_significant = only_significant)
  
  # output:
  pairs <- list("unadjusted_pairs"=unadjusted_observed_pairs,
                "Bonferroni_pairs"=Bonferroni_m_output) 
  
  # if numbers of repetitions of shuffling is more than 0 the shuffling should be performed
  if (shuffling) {
    adjusted_observed <- adjust_visits(unadjusted_pairs = unadjusted_observed_pairs,
                                       ids = ids,
                                       prob = prob)
    Shuffling_simulation_output <- replicate(n_shuffl,
                                             Shuffling_simulation(db_dates,
                                                                  ids=ids,
                                                                  adjusted_observed = adjusted_observed),
                                             simplify=FALSE)
    # determining the threshold
    for_threshold <- do.call("rbind", Shuffling_simulation_output)
    
    # take FP value:
    for_threshold_max <- apply(X = for_threshold,
                               MARGIN=2,
                               FUN=FP_threshold,
                               na.rm=TRUE)
    
    # find lowest threshold below alpha
    Threshold_lowest <- match(for_threshold_max[for_threshold_max<alpha][1], for_threshold_max)
    
    # extract the pairs above the threshold
    select_shuffled <- adjusted_observed[which(adjusted_observed$Corrected >Threshold_lowest),,drop=FALSE]
    
    pairs[["adjusted_pairs"]] <- adjusted_observed
    pairs[["Shuffled_pairs"]] <- select_shuffled
    pairs[["shuffling_threshold"]] <- Threshold_lowest
    pairs[["shuffling_simulation_output"]] <- Shuffling_simulation_output
  }
return(pairs)
}
find_pairs <- compiler::cmpfun(find_pairs)


