source("bma_pocrm_post.R")

## Simulation 4.2 Setup --------------------------------------------------------------
  
  # number of doses 
  d <- 6
  # number of partial orderings
  s <- 6
  
  # partial orderings 
  orders<-matrix(nrow=s,ncol=d)
  orders[1,] <- c(1,2,3,4,5,6)  
  orders[2,] <- c(1,2,3,5,4,6)  
  orders[3,] <- c(1,2,3,5,6,4)  
  orders[4,] <- c(1,2,5,3,4,6)  
  orders[5,] <- c(1,2,5,3,6,4)  
  orders[6,] <- c(1,2,5,6,3,4)
  
  # obtaining skeleton probabilities
  skeleton <- round(getprior(halfwidth = 0.08, target = 0.33, nu = 3, nlevel = 6), 2)
  
  # reordering skeleton probabilities according to partial orderings
  p.skel<-matrix(0,nrow=s,ncol=d)
  for(j in 1:s){
    p.skel[j,] <- skeleton[order(orders[j,])]
  }
  
  # true toxicity probabilities for various scenarios 
  p1 <- c(0.26,0.33,0.51,0.62,0.78,0.86)
  p2 <- c(0.12,0.21,0.34,0.50,0.66,0.79)
  p3 <- c(0.04,0.07,0.20,0.33,0.55,0.70)
  p4 <- c(0.01,0.04,0.05,0.17,0.33,0.67)
  p5 <- c(0.01,0.02,0.05,0.15,0.20,0.33)
  p6 <- c(0.50,0.60,0.70,0.80,0.90,0.95)
  
  ttl <- 0.33   ##target toxicity rate 
  cohortsize <- 1 ##cohort size for each inclusion
  ncohort <- 26  ##number of cohorts
  start_comb <- 3  ##starting combination
  n_stop <- 100   ##Number of patients needed on one combination to stop the trial
  ntrial <- 1000   ##number of simulated trials 
  cs <- 0.90        ##confidence level for the confidence interval 

## 4.2 Simulation --------------------------------------------------------------

# Scenario 1
set.seed(580)
bma_pocrm.sim(p0 = p1, 
           p.skel = p.skel, 
           ttl = ttl, 
           cohortsize = cohortsize, 
           ncohort = ncohort, 
           n_stop = n_stop, 
           ntrial = ntrial, 
           start_comb = start_comb,
           method = "Wages",
           cs = cs,
           conserv = T)

# Scenario 2
set.seed(580)
bma_pocrm.sim(p0 = p2, 
           p.skel = p.skel, 
           ttl = ttl, 
           cohortsize = cohortsize, 
           ncohort = ncohort, 
           n_stop = n_stop, 
           ntrial = ntrial, 
           start_comb = start_comb,
           method = "Wages",
           cs = cs,
           conserv = T)

# Scenario 3
set.seed(580)
bma_pocrm.sim(p0 = p3, 
           p.skel = p.skel, 
           ttl = ttl, 
           cohortsize = cohortsize, 
           ncohort = ncohort, 
           n_stop = n_stop, 
           ntrial = ntrial,
           start_comb = start_comb,
           method = "Wages",
           cs = cs,
           conserv = T)

# Scenario 4
set.seed(580)
bma_pocrm.sim(p0 = p4, 
           p.skel = p.skel, 
           ttl = ttl, 
           cohortsize = cohortsize,
           ncohort = ncohort,
           n_stop = n_stop, 
           ntrial = ntrial, 
           start_comb = start_comb,
           method = "Wages",
           cs = cs,
           conserv = T)

# Scenario 5
set.seed(580)
bma_pocrm.sim(p0 = p5, 
           p.skel = p.skel,
           ttl = ttl, 
           cohortsize = cohortsize, 
           ncohort = ncohort,
           n_stop = n_stop, 
           ntrial = ntrial,
           start_comb = start_comb,
           method = "Wages",
           cs = cs,
           conserv = T)

set.seed(580)
bma_pocrm.sim(p0 = p6, 
           p.skel = p.skel,
           ttl = ttl, 
           cohortsize = cohortsize, 
           ncohort = ncohort,
           n_stop = n_stop, 
           ntrial = ntrial,
           start_comb = start_comb,
           method = "Wages",
           cs = cs,
           conserv = T)

## Simulation 4.3 Setup --------------------------------------------------------------

ttl1 <- 0.3 ## target toxicity rate for scenarios 1-4
ttl5 <- 0.4 ## target toxicity rate for scenarios 5-6
cohortsize <- 3
ncohort <- 20
start_comb1 <- 8 ## starting combination for scenarios 1-4
start_comb5 <- 10 ## starting combination for scnearios 5-6
n_stop <- 100
ntrial <- 2000
cs <- 0.9
delta <- 0.08

d1 <- 16 ## number of doses in scenarios 1-2
d3 <- 15 ## number of doses in scenarios 3-4
d5 <- 20 ## number of doses in scenarios 5-6
s <- 3 ## number of partial orderings

## partial orders for scenarios 1-2
orders1 <- matrix(nrow = s, ncol = d1)
orders1[1,] <- c(1, 2, 5, 3, 6, 9, 4, 7, 10, 13, 8, 11, 14, 12, 15, 16)  
orders1[2,] <- c(1, 5, 2, 3, 6, 9, 13, 10, 7, 4, 8, 11, 14, 15, 12, 16)
orders1[3,] <- c(1, 5, 2, 9, 6, 3, 13, 10, 7, 4, 14, 11, 8, 15, 12, 16) 

## partial orders for scenarios 3-4
orders3 <- matrix(nrow = s, ncol = d3)
orders3[1,] <- c(1, 2, 4, 3, 5, 7, 6, 8, 10, 9, 11, 13, 12, 14, 15)  
orders3[2,] <- c(1, 2, 4, 7, 5, 3, 6, 8, 10, 13, 11, 9, 12, 14, 15)  
orders3[3,] <- c(1, 4, 2, 7, 5, 3, 10, 8, 6, 13, 11, 9, 14, 12, 15)  

orders5<-matrix(nrow=s, ncol=d5)
orders5[1,] <- c(1, 2, 5, 3, 6, 9, 4, 7, 10, 13, 8, 11, 14, 17, 12, 15, 18, 16, 19, 20)  
orders5[2,] <- c(1, 2, 5, 9, 6, 3, 4, 7, 10, 13, 17, 14, 11, 8, 12, 15, 18, 19, 16, 20) 
orders5[3,] <- c(1, 5, 2, 9, 6, 3, 13, 10, 7, 4, 17, 14, 11, 8, 18, 15, 12, 19, 16, 20)  

# Skeleton probabilities for scenarios 1-2.
skeleton1 <- getprior(halfwidth = 0.03, target = ttl1, nu = start_comb1, nlevel = d1)

# Skeleton probabilities for scenarios 3-4.
skeleton3 <- getprior(halfwidth = 0.03, target = ttl1, nu = start_comb1, nlevel = d3)

# Skeleton probabilities for scenarios 5-6.
skeleton5 <- getprior(halfwidth = 0.03, target = ttl5, nu = start_comb5, nlevel = d5)

# reordering skeleton probabilities according to partial orderings

p.skel1 <- matrix(0, nrow=s, ncol=d1)
for(j in 1:s){
  p.skel1[j,] <- skeleton1[order(orders1[j,])]
}

p.skel3 <- matrix(0, nrow = s, ncol = d3)
for(j in 1:s){
  p.skel3[j,] <- skeleton3[order(orders3[j,])]
}

p.skel5 <- matrix(0, nrow = s, ncol = d5)
for(j in 1:s){
  p.skel5[j,] <- skeleton5[order(orders5[j,])]
}

# true toxicity probabilities for scenarios 1-6
p1 <- c(0.06, 0.1, 0.15, 0.5, 
        0.08, 0.12, 0.3, 0.55,
        0.1,  0.3, 0.5,  0.6, 
        0.15, 0.45, 0.6, 0.7)
p2 <- c(0.08, 0.1, 0.12, 0.3,  0.12, 0.15, 0.3,  0.5,
        0.16, 0.3, 0.5,  0.55, 0.18, 0.45, 0.55, 0.6)
p3 <- c(0.06, 0.12, 0.15, 0.08, 0.16, 0.3,
        0.1,  0.3,  0.5,  0.3,  0.5,  0.55,
        0.5,  0.55, 0.6)
p4 <- c(0.06, 0.1,  0.5, 0.1, 0.3, 0.6,
        0.15, 0.5,  0.7, 0.3, 0.7, 0.8,
        0.5,  0.8,  0.9)
p5 <- c(0.24, 0.4,  0.48, 0.54,
        0.4,  0.45, 0.59, 0.67,
        0.47, 0.59, 0.68, 0.75,
        0.56, 0.67, 0.75, 0.81,
        0.64, 0.74, 0.81, 0.86)
p6 <- c(0.18, 0.27, 0.40, 0.49,
        0.29, 0.40, 0.49, 0.58,
        0.40, 0.45, 0.59, 0.68,
        0.47, 0.59, 0.68, 0.75,
        0.56, 0.67, 0.75, 0.81)
## Simulation 4.3 --------------------------------------------------------------

# Scenario 1
set.seed(580)
bma_pocrm.sim(p0 = p1, 
           p.skel = p.skel1, 
           ttl = ttl1, 
           cohortsize = cohortsize, 
           ncohort = ncohort, 
           n_stop = n_stop, 
           ntrial = ntrial, 
           start_comb = start_comb1,
           method = "Wages",
           cs = cs,
           conserv = T)

# Scenario 2
set.seed(580)
bma_pocrm.sim(p0 = p2, 
           p.skel = p.skel1, 
           ttl = ttl1, 
           cohortsize = cohortsize, 
           ncohort = ncohort, 
           n_stop = n_stop, 
           ntrial = ntrial, 
           start_comb = start_comb1,
           method = "Wages",
           cs = cs,
           conserv = T)

# Scenario 3 
set.seed(580)
bma_pocrm.sim(p0 = p3, 
           p.skel = p.skel3, 
           ttl = ttl1, 
           cohortsize = cohortsize, 
           ncohort = ncohort, 
           n_stop = n_stop, 
           ntrial = ntrial, 
           start_comb = start_comb1,
           method = "Wages",
           cs = cs,
           conserv = T)

# Scenario 4
set.seed(580)
bma_pocrm.sim(p0 = p4, 
           p.skel = p.skel3, 
           ttl = ttl1, 
           cohortsize = cohortsize,
           ncohort = ncohort, 
           n_stop = n_stop, 
           ntrial = ntrial, 
           start_comb = start_comb1,
           method = "Wages",
           cs = cs,
           conserv = T)

# Scenario 5
set.seed(580)
bma_pocrm.sim(p0 = p5, 
           p.skel = p.skel5, 
           ttl = 0.4, 
           cohortsize = 3, 
           ncohort = 20, 
           n_stop = 100, 
           ntrial = 2000, 
           start_comb = 10,
           method = "Wages",
           cs = cs,
           conserv = T)

# Scenario 6
set.seed(580)
bma_pocrm.sim(p0 = p6, 
           p.skel = p.skel5,
           ttl = ttl5, 
           cohortsize = cohortsize,
           ncohort = ncohort,
           n_stop = n_stop, 
           ntrial = ntrial, 
           start_comb = 10,
           method = "Wages",
           cs = cs,
           conserv = T)

## Question 7  -----------------------------------------------------------------

# Effect of change in dose allocation on results

# Scenario 1
set.seed(580)
bma_pocrm.sim_(p0 = p1, 
               p.skel = p.skel, 
               ttl = 0.33, 
               cohortsize = 1, 
               ncohort = 22, 
               n_stop = 100, 
               ntrial = 1000, 
               start_comb = 3,
               method = "Wages",
               cs = 0.9, 
               conserv = T)

# Scenario 2
set.seed(580)
bma_pocrm.sim(p0 = p2, 
               p.skel = p.skel, 
               ttl = 0.33, 
               cohortsize = 1, 
               ncohort = 26, 
               n_stop = 100, 
               ntrial = 1000, 
               start_comb = 3,
               method = "Wages",
               cs = 0.9, 
               conserv = T)
# Scenario 3
set.seed(580)
bma_pocrm.sim(p0 = p3, 
               p.skel = p.skel, 
               ttl = 0.33, 
               cohortsize = 1, 
               ncohort = 29, 
               n_stop = 100, 
               ntrial = 1000, 
               start_comb = 3,
               method = "Wages",
               cs = 0.9, 
               conserv = T)

# Scenario 4
set.seed(580)
bma_pocrm.sim(p0 = p4, 
               p.skel = p.skel, 
               ttl = 0.33, 
               cohortsize = 1, 
               ncohort = 29, 
               n_stop = 100, 
               ntrial = 1000, 
               start_comb = 3,
               method = "Wages",
               cs = 0.9, 
               conserv = T)

# Scenario 5
set.seed(580)
bma_pocrm.sim(p0 = p5, 
               p.skel = p.skel, 
               ttl = 0.33, 
               cohortsize = 1, 
               ncohort = 27, 
               n_stop = 100, 
               ntrial = 1000, 
               start_comb = 3,
               method = "Wages",
               cs = 0.9, 
               conserv = T)

## Question 8 ------------------------------------------------------------------

# Check different values of c e.g. 0.10, 0.25, 0.40

# Scenario 1
set.seed(580)
bma_pocrm.sim(p0 = p1, 
               p.skel = p.skel, 
               ttl = 0.33, 
               cohortsize = 1, 
               ncohort = 22, 
               n_stop = 100, 
               ntrial = 1000, 
               start_comb = 3,
               method = "Wages",
               cs = 0.9, 
               c_od = .5)

# Scenario 2
set.seed(580)
bma_pocrm.sim(p0 = p2, 
               p.skel = p.skel, 
               ttl = 0.33, 
               cohortsize = 1, 
               ncohort = 26, 
               n_stop = 100, 
               ntrial = 1000, 
               start_comb = 3,
               method = "Wages",
               cs = 0.9, 
               c_od = .5)
# Scenario 3
set.seed(580)
bma_pocrm.sim(p0 = p3, 
               p.skel = p.skel, 
               ttl = 0.33, 
               cohortsize = 1, 
               ncohort = 29, 
               n_stop = 100, 
               ntrial = 1000, 
               start_comb = 3,
               method = "Wages",
               cs = 0.9, 
               c_od = .5)

# Scenario 4
set.seed(580)
bma_pocrm.sim(p0 = p4, 
               p.skel = p.skel, 
               ttl = 0.33, 
               cohortsize = 1, 
               ncohort = 29, 
               n_stop = 100, 
               ntrial = 1000, 
               start_comb = 3,
               method = "Wages",
               cs = 0.9, 
               c_od = .5)

# Scenario 5
set.seed(580)
bma_pocrm.sim(p0 = p5, 
               p.skel = p.skel, 
               ttl = 0.33, 
               cohortsize = 1, 
               ncohort = 27, 
               n_stop = 100, 
               ntrial = 1000, 
               start_comb = 3,
               method = "Wages",
               cs = 0.9, 
               c_od = .5)

# Scenario 6
set.seed(580)
bma_pocrm.sim(p0 = p6, 
               p.skel = p.skel, 
               ttl = 0.33, 
               cohortsize = 1, 
               ncohort = 27, 
               n_stop = 100, 
               ntrial = 1000, 
               start_comb = 3,
               method = "Wages",
               cs = 0.9, 
               c_od = .5)
