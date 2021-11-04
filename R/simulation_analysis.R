library(ggplot2)

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

scenario <- as.matrix(rbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13))

data <- c()

for(i in 1:13){
   result <- fromJSON(file = paste0("output_scenario_",i,".Rout"))
   
   # percent correct selections 
   pcs_scen <- rbind(sum(result[[1]]$select_percent[which(scenario[i,] == 0.3)]), 
                     sum(result[[2]]$select_percent[which(scenario[i,] == 0.3)]), 
                     sum(result[[3]]$select_percent[which(scenario[i,] == 0.3)]))
   
   # percent accetable selections 
   pas_scen <- rbind(sum(result[[1]]$select_percent[which(scenario[i,] == 0.3 | scenario[i,] == 0.2)]), 
                       sum(result[[2]]$select_percent[which(scenario[i,] == 0.3 | scenario[i,] == 0.2)]), 
                       sum(result[[3]]$select_percent[which(scenario[i,] == 0.3 | scenario[i,] == 0.2)]))
   
   # percent overly toxic selections 
   pots_scen <- rbind(sum(result[[1]]$select_percent[which(scenario[i,] > 0.33)]), 
                         sum(result[[2]]$select_percent[which(scenario[i,] > 0.33)]), 
                         sum(result[[3]]$select_percent[which(scenario[i,] > 0.33)]))

   # number of patients treated at doses considered to be overly toxic
   unsafe_scen <- rbind(sum(result[[1]]$num_pat[which(scenario[i,] > 0.33)])*0.2, 
                          sum(result[[2]]$num_pat[which(scenario[i,] > 0.33)])*0.2, 
                          sum(result[[3]]$num_pat[which(scenario[i,] > 0.33)])*0.2)
   
   data <- cbind(pcs_scen,pas_scen,pots_scen, unsafe_scen, rep(i,3), c("POCRM","BMA Mix. Dens.","BMA Mix. Mean")) %>%
      rbind(data)
}

colnames(data) <- c("PCS", "PAS", "POTS", "over_tox","Scenario","Method")
data <- transform(data, PCS = as.numeric(levels(PCS))[PCS], PAS = as.numeric(levels(PAS))[PAS], POTS = as.numeric(levels(POTS))[POTS], over_tox = as.numeric(levels(over_tox))[over_tox])

mean_data <- data %>% 
   group_by(Method) %>% 
   summarise(PCS = mean(PCS), PAS = mean(PAS), POTS = mean(POTS), over_tox = mean(over_tox)) %>%
   cbind(Scenario = rep("Mean",3))

data <- rbind(data, mean_data)

ggplot(data, mapping = aes(x = Scenario, y = PCS, fill = Method)) +
   geom_col(position = "dodge") +
   scale_x_discrete(limits = factor(c(seq(1,13,1), "Mean"))) +
   ylim(0,100) +
   labs(y = "%", title = "Proportion of correct selections vs. scenario")

ggplot(data, mapping = aes(x = Scenario, y = PAS, fill = Method)) +
   geom_col(position = "dodge") +
   scale_x_discrete(limits = factor(c(seq(1,13,1), "Mean"))) +
   ylim(0,100) + 
   labs(y = "%", title = "Proportion of acceptable selections vs. scenario")

ggplot(data, mapping = aes(x = Scenario, y = POTS, fill = Method)) +
   geom_col(position = "dodge") +
   scale_x_discrete(limits = factor(c(seq(1,13,1), "Mean"))) +
   ylim(0,100) +
   labs(y = "%", title = "Proportion of overly toxic selections vs. scenario")

ggplot(data, mapping = aes(x = Scenario, y = over_tox, fill = Method)) +
   geom_col(position = "dodge") +
   scale_x_discrete(limits = factor(c(seq(1,13,1), "Mean"))) +
   labs(y = "Count", title = "Number of patients treated at overly toxic doses vs. scenario")
