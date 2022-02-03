source("bma_pocrm_post.R")

# Number of combinations 
d <- 6

# Number of potential orderings
s <- 3

# Specify the possible toxicity orderings of the drug combinations
orders<-matrix(nrow=s,ncol=d)
orders[1,]<-c(1,2,3,4,5,6)
orders[2,]<-c(1,3,2,4,6,5)
orders[3,]<-c(1,4,3,2,6,5)

# skeleton probability specification
skeleton <- getprior(halfwidth = 0.08, target = 0.4, nu = 3, nlevel = d)

# reordering of skeleton probabilities for each candidate
p_skel <- matrix(nrow = s, ncol = d)
for(j in 1:s){
  p_skel[j,] <- skeleton[order(orders[j,])]
}

# guiding true toxicity probabilities
p1 <- c(0.3,0.5,0.6,0.4,0.7,0.8)

# trial data
n1 <- c(0, 2, 1, 0, 0, 0)
y1 <- c(0, 1, 1, 0, 0, 0)

n2 <- c(1, 2, 1, 0, 0, 0)
y2 <- c(0, 1, 1, 0, 0, 0)

# first run

pocrm_run1 <- pocrm.imp(p.skel = p_skel, 
          ttl = 0.4, 
          y = y1, 
          n = n1)

bma_run1 <- bma_pocrm_post2.imp(p.skel = p_skel,
                    ttl = 0.4,
                    y = y1,
                    n = n1)

pocrm_run2 <- pocrm.imp(p.skel = p_skel, 
          ttl = 0.4, 
          y = y2, 
          n = n2)

bma_run2 <- bma_pocrm_post2.imp(p.skel = p_skel,
                    ttl = 0.4,
                    y = y2,
                    n = n2)

# Test:
# escalation: 3 -> 2 -> 2 -> 1
# ntest <- c(0, 2, 1, 0, 0, 0)
# ytest <- c(0, 1, 1, 0, 0, 0)
# pocrm.imp(p.skel = p_skel,
#           ttl = 0.4, 
#           y = ytest,
#           n = ntest)
