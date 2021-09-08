source("bma_pocrm_post.R")
library(rjson)

task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id <- as.numeric(task_id_string)

# possible combinations of scenario, delta (as indicies of scenario matrix, delta vector)
comb <- matrix(0, nrow = 72, ncol = 2)
for(i in 1:72){
  if(i < 19){
    comb[i,] <- c(1,i %% 18)
  } else if (i < 37){
    comb[i,] <- c(2,i %% 18)
  } else if (i < 55){
    comb[i,] <- c(3,i %% 18)
  } else {
    comb[i,] <- c(4,i %% 18)
  }
}
comb[which(comb[,2] == 0),2] <- 18

sim_comb <- comb[task_id,]

## Simulation Parameters -------------------------------------------------------

ttl <- 0.3    # Target toxicity rate 
cohortsize <- 1 # Cohort size for each inclusion
ncohort <- 20   # Number of cohorts
start_comb <- 5 # Starting combination
n_stop <- 100    # Number of patients needed on one combination to stop the trial
ntrial <- 1000   # Number of simulated trials 

# number of potential orderings
s <- 6

# number of doses
d <- 9

# half-width used to generate skeleton probabilities
delta <- seq(0.01,0.6*ttl,0.01)[sim_comb[2]]

# true toxicity probabilities for scenarios 1, 8, 10 and 13
p1 <- c(0.05, 0.10, 0.15, 0.10, 0.15, 0.20, 0.15, 0.20, 0.30)
p8 <- c(0.05, 0.10, 0.20, 0.10, 0.20, 0.30, 0.30, 0.45, 0.55)
p10 <- c(0.15, 0.30, 0.45, 0.30, 0.45, 0.55, 0.45, 0.55, 0.65)
p13 <- c(0.30, 0.45, 0.50, 0.45, 0.50, 0.55, 0.50, 0.55, 0.60)
scenario <- as.matrix(rbind(p1,p8,p10,p13))

# generating skeleton
skeleton <- getprior(halfwidth = delta, target = ttl, nu = 5, nlevel = 9)

# potential orderings
orders <- matrix(0, nrow = s, ncol = d)
orders[1,] <- c(1,2,3,4,5,6,7,8,9)
orders[2,] <- c(1,4,7,2,5,8,3,6,9)
orders[3,] <- c(1,4,2,7,5,3,8,6,9)
orders[4,] <- c(1,2,4,3,5,7,6,8,9)
orders[5,] <- c(1,2,4,7,5,3,6,8,9)
orders[6,] <- c(1,4,2,3,5,7,8,6,9)

p.skel <- matrix(0, nrow = s, ncol = d)
for(j in 1:s){
  p.skel[j,] <- skeleton[order(orders[j,])]
}

## Simulations -----------------------------------------------------------------

result <- c()

# Wages et al. Simulation
set.seed(580)
result[[1]] <- bma_pocrm.sim(p0 = scenario[sim_comb[1],],
                             p.skel = p.skel,
                             ttl = ttl,
                             start_comb = start_comb,
                             cohortsize = cohortsize,
                             ncohort = ncohort,
                             n_stop = n_stop,
                             ntrial = ntrial,
                             method = "Wages",
                             sim_method = "Mozgunov",
                             output = "list")

# BMA method with alloction via densities
set.seed(580)
result[[2]] <- bma_pocrm.sim(p0 = scenario[sim_comb[1],],
                             p.skel = p.skel,
                             ttl = ttl,
                             start_comb = start_comb,
                             cohortsize = cohortsize,
                             ncohort = ncohort,
                             n_stop = n_stop,
                             ntrial = ntrial,
                             method = "BMA Mixture",
                             sim_method = "Mozgunov",
                             output = "list")

# BMA method with allocation via point estimates
set.seed(580)
result[[3]] <- bma_pocrm.sim(p0 = scenario[sim_comb[1],],
                             p.skel = p.skel,
                             ttl = ttl,
                             start_comb = start_comb,
                             cohortsize = cohortsize,
                             ncohort = ncohort,
                             n_stop = n_stop,
                             ntrial = ntrial,
                             method = "BMA Mixture2",
                             sim_method = "Mozgunov",
                             output = "list")

# saving to file
filename <- paste0("output_calib_", task_id, ".Rout")
cat(toJSON(result), file = filename)
