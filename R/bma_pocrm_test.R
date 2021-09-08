source("bma_pocrm_post.R")

# Numver of combinations 
d<-6

# Number of potential orderings
s<-6   

# Specify the possible toxicity orderings of the drug combinations
orders<-matrix(nrow=s,ncol=d)
orders[1,]<-c(1,2,3,4,5,6)  
orders[2,]<-c(1,2,3,5,4,6)  
orders[3,]<-c(1,2,3,5,6,4)  
orders[4,]<-c(1,2,5,3,4,6)  
orders[5,]<-c(1,2,5,3,6,4)  
orders[6,]<-c(1,2,5,6,3,4)

# Skeleton Toxicity Probabilities
skeleton<-c(0.07,0.18,0.33,0.49,0.63,0.75)

# The following will adjust the location of each skeleton value to correspond 
# to the 's' possible orderings specified above. 
# 'p.skel' is the matrix of toxicity skeleton values
p.skel<-matrix(0,nrow=s,ncol=d)
for(j in 1:s){
  p.skel[j,]<-skeleton[order(orders[j,])]
}

# target toxicity rate
ttr=0.33

# confidence level for the confidence interval
cs=0.90	

# number of DLT's in each combo
y=c(0,0,0,1,1,0)   

# number of patients treated on each combo
n=c(1,1,1,3,2,0)   

results <- pocrm.imp(p.skel = p.skel, ttl = ttr, y = y, n = n)
results

results <- bma_pocrm_point.imp(p.skel = p.skel, ttl = ttr, y = y, n = n)
results

results <- bma_pocrm_post.imp(p.skel = p.skel, ttl = ttr, y = y, n = n)
results

results <- bma_pocrm_post2.imp(p.skel = p.skel, ttl = ttr, y = y, n = n)
results


