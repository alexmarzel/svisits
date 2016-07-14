
# Load packages
library("compiler");library("parallel"); 
library("chron"); library("data.table"); library("Rcpp");
enableJIT(3);


# prepare the database
prepare_db<-function(your_database=simulated_data, ids_column="subject",dates_column="sim_dates")
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
adjust_visits<-function(unadjusted_pairs, prob=1/75){
  real_db = data.table::data.table(allPairs= unlist(names(unadjusted_pairs)), Freq=unlist(as.numeric(unadjusted_pairs)))
  real_db_1<-real_db[real_db$Freq >=1, ]
  real_db_1<-cbind(real_db_1,do.call("rbind",strsplit(as.character(real_db_1$allPairs),"_")))
  setnames(real_db_1,3:4,c("id_1","id_2")); setkey(real_db_1, "id_1");
  real_db_1<-na.omit(real_db_1[ids,allow.cartesian=TRUE]); setkey(real_db_1, "id_2");
  real_db_1<-na.omit(real_db_1[ids,allow.cartesian=TRUE]); setnames(real_db_1,5:6,c("N_visits.x","N_visits.y"))
  real_db_1$Freq<-as.numeric(real_db_1$Freq); real_db_1$N_visits.x<-as.numeric(real_db_1$N_visits.x);real_db_1$N_visits.y<-as.numeric(real_db_1$N_visits.y)  
  real_db_1<-cbind(real_db_1, Corrected=svisits::adjust_Rcpp_min(mat = as.matrix( real_db_1[, .(Freq, N_visits.x,N_visits.y)] ), prob=prob ))
  real_db_1<-as.data.frame(real_db_1)
  return(real_db_1)}
adjust_visits<-compiler::cmpfun(adjust_visits)

# Shuffling simulatoins
Shuffling_simulation<-function(db_dates,adjusted_observed){
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
  Adjusted_randomized<-adjust_visits(unadjusted_pairs =table(RandomPairs))
  # return vector of false-positive thresholds
  return(sapply(1:10,
                function(x) sum(Adjusted_randomized$Corrected>x)/sum(adjusted_observed$Corrected>x) )) 
}; Shuffling_simulation <-cmpfun(Shuffling_simulation)

# Bonfrerroni method
Bonferroni_m<-function(unadjusted_pairs,prob=1/75,alpha=0.01,only_significant=TRUE){
  real_db = data.table(allPairs= unlist(names(unadjusted_pairs)), Freq=unlist(as.numeric(unadjusted_pairs)))
  real_db_1<-real_db[real_db$Freq >=1, ]
  real_db_1<-cbind(real_db_1,do.call("rbind",strsplit(as.character(real_db_1$allPairs),"_")))
  setnames(real_db_1,3:4,c("id_1","id_2")); setkey(real_db_1, "id_1");
  real_db_1<-na.omit(real_db_1[ids,allow.cartesian=TRUE]); setkey(real_db_1, "id_2");
  real_db_1<-na.omit(real_db_1[ids,allow.cartesian=TRUE]); setnames(real_db_1,5:6,c("N_visits.x","N_visits.y"))
  real_db_1$Freq<-as.numeric(real_db_1$Freq); real_db_1$N_visits.x<-as.numeric(real_db_1$N_visits.x);real_db_1$N_visits.y<-as.numeric(real_db_1$N_visits.y)  
  real_db_1<-cbind(real_db_1, Prob_for_Bonferr=svisits::adjust_R_Rcpp_short_binomial(mat = as.matrix( real_db_1[, .(Freq, N_visits.x,N_visits.y)]), prob=prob ))
  real_db_1<-as.data.frame(real_db_1)
  
  real_db_1<-cbind(real_db_1,BP=1:length(real_db_1[,1]))
  real_db_1<-cbind(real_db_1, ltp=-log10(real_db_1$Prob_for_Bonferr))
  
  selected<-real_db_1[ real_db_1$ltp>=-log10(alpha/length(real_db_1[,1])), ]
  
  if(only_significant==TRUE) {return(selected)} else {return(real_db_1)}
}



