## Information required for analysis -------------------------------------------
geom_mean <- function(x) { 
    return(exp(mean(log(x))))
  }

library(rjson)
library(ggplot2)
library(tidyverse)

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

p1 <- c(0.05, 0.10, 0.15, 0.10, 0.15, 0.20, 0.15, 0.20, 0.30)
p8 <- c(0.05, 0.10, 0.20, 0.10, 0.20, 0.30, 0.30, 0.45, 0.55)
p10 <- c(0.15, 0.30, 0.45, 0.30, 0.45, 0.55, 0.45, 0.55, 0.65)
p13 <- c(0.30, 0.45, 0.50, 0.45, 0.50, 0.55, 0.50, 0.55, 0.60)
scenario <- as.matrix(rbind(p1,p8,p10,p13))

scenario_int <- c(1,8,10,13)

ttl <- 0.3
delta <- seq(0.01,0.6*ttl,0.01)

# values of delta to be analysed
# other values excluded due to errors in analysis
keep <- which(comb[,2] %in% c(1:8))

selection_data <- c()

# obtaining PCS results for all simulations 
for(id in c(1:67,69:72)){
  results <-fromJSON(file = paste0("output_calib_",id,".Rout"))

  selection_data <- rbind(sum(results[[1]]$select_percent[which(scenario[comb[id,1],] == 0.3)]), 
               sum(results[[2]]$select_percent[which(scenario[comb[id,1],] == 0.3)]), 
               sum(results[[3]]$select_percent[which(scenario[comb[id,1],] == 0.3)])) %>%
    cbind(c("POCRM", "BMA Mix. Dens.","BMA Mix. Mean"), rep(scenario_int[comb[id,1]],3), rep(delta[comb[id,2]],3)) %>%
    rbind(selection_data)
  
}

# formatting data
colnames(selection_data) <- c("PCS", "Method", "Scenario", "Delta")
selection_data <- transform(selection_data, PCS = as.numeric(levels(PCS))[PCS], Scenario = as.numeric(levels(Scenario))[Scenario], Delta = as.numeric(levels(Delta))[Delta])

mean_data <- selection_data %>% 
  group_by(Method, Delta) %>% 
  summarise(geom_mean = geom_mean(PCS))

ggplot(mean_data, mapping = aes(x = Delta, y = geom_mean, color = Method)) + 
  geom_line() +
  labs(y = "Geometric Mean PCS", title = "Geometric mean of proportion of correct selections vs. Delta")

