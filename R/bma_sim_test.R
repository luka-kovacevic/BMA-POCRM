source("bma_pocrm_post.R")

# Number of treatment combinations
d<-6

# Number of candidate toxicity orderings
s<-3  

# Specify the possible toxicity orderings of the drug combinations
orders<-matrix(nrow=s,ncol=d)
orders[1,]<-c(1,2,3,4,5,6)  
orders[2,]<-c(2,1,3,4,5,6)  
orders[3,]<-c(1,2,3,4,6,5)  

# Skeleton probabilties
skeleton<-c(0.07,0.18,0.33,0.49,0.63,0.75)

# reordering of skeleton probabilities to correspond to candidate orderings
p.skel<-matrix(0,nrow=s,ncol=d)
for(j in 1:s){
  p.skel[j,]<-skeleton[order(orders[j,])]
}

# True toxicity probability scenarios
p1<-c(0.08,0.09,0.23,0.45,0.63,0.64)
p2<-c(0.08,0.09,0.23,0.45,0.66,0.82)
p3<-c(0.04,0.18,0.31,0.47,0.66,0.82)

ttr=0.4     # Target toxicity rate 
cohortsize=1 # Cohort size for each inclusion
ncohort=22   # Number of cohorts
start_comb=1  # Starting combination
n_stop=100    # Number of patients needed on one combination to stop the trial
ntrial=1000   # Number of simulated trials 

set.seed(580)
bma_pocrm.sim(p0 = p1,
              p.skel = p.skel,
              ttl = ttl,
              start_comb = start_comb,
              cohortsize = cohortsize,
              ncohort = ncohort,
              n_stop = n_stop,
              ntrial = ntrial,
              method = "Wages")

bma_pocrm.sim(p0 = p1,
              p.skel = p.skel,
              ttl = ttl,
              start_comb = start_comb,
              cohortsize = cohortsize,
              ncohort = ncohort,
              n_stop = n_stop,
              ntrial = ntrial,
              method = "BMA Point")

bma_pocrm.sim(p0 = p1,
              p.skel = p.skel,
              ttl = ttl,
              start_comb = start_comb,
              cohortsize = cohortsize,
              ncohort = ncohort,
              n_stop = n_stop,
              ntrial = ntrial,
              method = "BMA Mixture")

bma_pocrm.sim(p0 = p1,
              p.skel = p.skel,
              ttl = ttl,
              start_comb = start_comb,
              cohortsize = cohortsize,
              ncohort = ncohort,
              n_stop = n_stop,
              ntrial = ntrial,
              method = "BMA Mixture2")
