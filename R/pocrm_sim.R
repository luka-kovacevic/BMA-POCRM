source("bma_pocrm_post.R")
library(rjson)

task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id <- as.numeric(task_id_string)

# number of potential orderings
s <- 6

# number of doses
d <- 9

# half-width used to generate skeleton probabilities
delta <- 0.08

# true toxicity probabilities for scenarios 1-15
p1 <- c(0.05, 0.10, 0.15, 0.10, 0.15, 0.20, 0.15, 0.20, 0.30)
p2 <- c(0.05, 0.10, 0.15, 0.10, 0.20, 0.30, 0.20, 0.30, 0.45)
p3 <- c(0.02, 0.05, 0.10, 0.10, 0.15, 0.20, 0.20, 0.30, 0.45)
p4 <- c(0.05, 0.10, 0.15, 0.10, 0.20, 0.30, 0.20, 0.45, 0.60)
p5 <- c(0.02, 0.05, 0.15, 0.20, 0.30, 0.45, 0.45, 0.55, 0.65)
p6 <- c(0.10, 0.15, 0.30, 0.15, 0.30, 0.45, 0.30, 0.45, 0.60)
p7 <- c(0.10, 0.20, 0.45, 0.15, 0.30, 0.50, 0.30, 0.50, 0.60)
p8 <- c(0.05, 0.10, 0.20, 0.10, 0.20, 0.30, 0.30, 0.45, 0.55)
p9 <- c(0.10, 0.15, 0.30, 0.30, 0.40, 0.50, 0.40, 0.50, 0.60)
p10 <- c(0.15, 0.30, 0.45, 0.30, 0.45, 0.55, 0.45, 0.55, 0.65)
p11 <- c(0.02, 0.05, 0.10, 0.30, 0.45, 0.60, 0.45, 0.60, 0.75)
p12 <- c(0.20, 0.30, 0.45, 0.45, 0.50, 0.55, 0.65, 0.70, 0.75)
p13 <- c(0.30, 0.45, 0.50, 0.45, 0.50, 0.55, 0.50, 0.55, 0.60)
p14 <- c(0.45, 0.50, 0.55, 0.50, 0.55, 0.60, 0.55, 0.60, 0.65)
p15 <- rep(0.10, 9)

scenario <- as.matrix(rbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15))

skeleton <- getprior(halfwidth = delta, target = 0.3, nu = 5, nlevel = 9)

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

## Simulation Parameters -------------------------------------------------------
ttl <- 0.3    # Target toxicity rate 
cohortsize <- 1 # Cohort size for each inclusion
ncohort <- 20   # Number of cohorts
start_comb <- 5 # Starting combination
n_stop <- 100    # Number of patients needed on one combination to stop the trial
ntrial <- 1000   # Number of simulated trials 

result <- c()

set.seed(580)
result[[1]] <- bma_pocrm.sim(p0 = scenario[task_id,],
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

set.seed(580)
result[[2]] <- bma_pocrm.sim(p0 = scenario[task_id,],
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

set.seed(580)
result[[3]] <- bma_pocrm.sim(p0 = scenario[task_id,],
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

filename <- paste0("output_scenario_", task_id, ".Rout")
cat(toJSON(result), file = filename)




