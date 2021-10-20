library(nnet)
library(dfcrm)
library(numDeriv)

#' Expression for posterior density of a 
#' 
#' @param a value of the parameter a
#' @param p skeleton probabilities
#' @param y number of observed toxicities at each dose
#' @param n number of observed patients at each dose
#' 
bcrmh<-function(a,p,y,n){
  s2=1.34
  lik=exp(-0.5*a*a/s2)
  for(j in 1:length(p)){ # for each dose level
    pj=p[j]**exp(a)
    lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]); # uses binomial distribution
    # i.e. bernoulli distribution for each individual patient in the trial
  }
  return(lik);
}

#' Used to calculate the posterior mean of a 
#' 
bcrmht<-function(a,p,y,n){
  s2=1.34
  lik=a*exp(-0.5*a*a/s2)
  for(j in 1:length(p)){
    pj=p[j]**exp(a)
    lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
  }
  return(lik);
}

#' Used to calculate the posterior mean of a^2
#' 
bcrmht2<-function(a,p,y,n){
  s2=1.34
  lik=a^2*exp(-0.5*a*a/s2)
  for(j in 1:length(p)){
    pj=p[j]**exp(a)
    lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
  }
  return(lik);
}

#' Uses POCRM proposed by Wages et al. (2011) to recommend next dose in trial
#' based on available data
#' 
#' @param p.skel matrix of skeleton probabilities for all possible orderings
#' @param ttl target toxicity level
#' @param y number of DLTs at each dose
#' @param n number of patients treated at each dose
#' @param cs in (0,1): confidence level for confidence intervals
#' @param m_prior prior model probabilities, if unspecified uniform probabilities are assigned
#' @param c_od constant for overdose controlling, if unspecified overdose controlling not applied
#' @param conserv boolean for whether dose allocation is conservative:
#'  - if TRUE, assign only doses with probability of toxicity less than the target toxicity rate
#' 
#' @return list containing the recommended next dose, posterior model probabilities,
#' indicator for stopping due to safety concerns, estimated probability of toxicities for each dose
#' 
pocrm.imp<-function(p.skel, ttl, y, n, cs = 0.9, m_prior = NA, c_od = NA, conserv = F){
  
  # if single ordering is given, convert to matrix
  
  if (is.vector(p.skel)) p.skel <- t(as.matrix(p.skel))
  
  # number of potential orderings
  
  nord.tox <- nrow(p.skel)
  
  # assigns default, uniform prior for model probabilities if necessary
  
  if (is.na(m_prior)) {
    mprior.tox <- rep(1/nord.tox, nord.tox);
  } else {
    mprior.tox <- mprior
  }
  
  # number of drug combinations
  
  ncomb = ncol(p.skel)
  
  # indicator for stopping due to safety
  # stop = 1 indicates safety concern
  
  stop <- 0
  
  # initialising vector for integral results for each ordering
  
  marginal.tox <- est.tox <- e2.tox <- rep(0, nord.tox);
  
  for (k in 1:nord.tox) {
    
    # normalising constant for each ordering
    
    marginal.tox[k] <- integrate(bcrmh, lower = -Inf, upper = Inf, p = p.skel[k,], y = y, n = n, abs.tol = 0)$value;
    
    # posterior mean of a for each ordering
    
    est.tox[k] <- integrate(bcrmht, lower=-10, upper=10, p.skel[k,], y, n, abs.tol = 0)$value/marginal.tox[k]
    
    # posterior mean of a^2 for each ordering
    
    e2.tox[k] <- integrate(bcrmht2, lower=-10, upper=10, p.skel[k,], y, n, abs.tol = 0)$value/marginal.tox[k]
  }		
  
  # posterior model probabilities
  
  postprob.tox <- (marginal.tox*mprior.tox)/sum(marginal.tox*mprior.tox)
  
  # model selection, identify the model with the highest posterior probability
  
  if (nord.tox > 1) { 
    mtox.sel <- which.is.max(postprob.tox)
  } else {
    mtox.sel <- 1
  }
  
  # estimated probability of toxicity for each dose
  
  ptox_hat <- p.skel[mtox.sel,]**exp(est.tox[mtox.sel])
  
  # checks if overdose controlling is applied
  
  if(is.na(c_od)){
    
    # checks if safety constraint is applied
    
    if (!is.na(cs)) {
      
      # variance of posterior distribution of a
      
      post.var.tox <- e2.tox[mtox.sel] - (est.tox[mtox.sel])^2
      
      # critical value for given confidence level 
      
      crit.tox <- qnorm(0.5+cs/2)
      
      # lower and upper bound for a based on critical value
      
      lb.tox <- est.tox[mtox.sel]-crit.tox*sqrt(post.var.tox)
      ub.tox <- est.tox[mtox.sel]+crit.tox*sqrt(post.var.tox)
      
      # lower and upper bound for probability of toxicity
      
      ptoxL <- p.skel[mtox.sel,]**exp(ub.tox)
      ptoxU <- p.skel[mtox.sel,]**exp(lb.tox)
      
      # checks safety contraint
      
      if (ptoxL[1]>ttl) {
        stop <- 1
      }	
    }
    
    # checks if allocations should be conservative
    
    if(!conserv){
      
      # assigns dose with estimated probability of toxicity closest to ttl
      
      distance <- abs(ptox_hat-ttl)
      comb_curr <- which.is.max(-distance)

    } else {
      
      # assigns dose with estimated probabilitiy of toxicity closest to ttl while being lower than ttl
      
      distance <- ptox_hat - ttl
      distance[which(distance >= 0)] <- -10000
      comb_curr <- which.is.max(distance)

    }
    
  } else {
    
    # posterior density of a for a particular ordering
    
    post <- function(a, model){
      return(bcrmh(a = a, p = p.skel[model,], y = y, n = n)/marginal.tox[model])
    }
    
    # overdose controlling 
    
    prob_od <- numeric(ncomb)
    for(k in 1:ncomb){
      a_star <- log(log(ttl)/log(p.skel[mtox.sel, k]))
      prob_od[k] <- integrate(post, model = mtox.sel, lower = -Inf, upper = a_star)$value
    }
    
    # checks if dose allocation should be conservative
    
    if(!conserv){
      
    # assigns dose with estimated probability of toxicity closest to ttl which satisfies overdose control
      
    distance <- abs(ptox_hat - ttl)
    distance[which(prob_od >= c_od)] <- 10000
    comb_curr <- which.is.max(-distance)
    
    } else {
      
      # assigns dose with estimated probability of toxicity closest to ttl which satisfies overdose control and is lower than the ttl
      
      distance <- ptox_hat - ttl
      distance[which(prob_od >= c_od | distance >= 0)] <- -10000
      comb_curr <- which.is.max(distance)
      
    }
    
  }
  
  return(list(comb_curr = comb_curr,
              m_post = postprob.tox,
              stop = stop,
              ptox_hat = ptox_hat))
}

#' Uses Bayesian model averaging in the POCRM to recommend next dose in trial based
#' on available data
#' 
#' @param p.skel matrix of skeleton probabilities for all possible orderings
#' @param ttl target toxicity level
#' @param y number of DLTs at each dose
#' @param n number of patients treated at each dose
#' @param cs in (0,1): confidence level for confidence intervals
#' @param m_prior prior model probabilities, if unspecified uniform probabilities are assigned
#' @param c_od constant for overdose controlling, if unspecified overdose controlling not applied
#' @param conserv boolean for whether dose allocation is conservative:
#'  - if TRUE, assign only doses with probability of toxicity less than the target toxicity rate
#' 
#' @return list containing the recommended next dose, posterior model probabilities,
#' indicator for stopping due to safety concerns, estimated probability of toxicities for each dose
#' 
bma_pocrm_point.imp <- function(p.skel, ttl, y, n, m_prior = NA, c_od = NA, conserv = F){
  
  # if single ordering is given, convert vector to matrix
  
  if(is.vector(p.skel)){ 
    p.skel <- t(as.matrix(p.skel)) 
  }
  
  # number of potential orderings
  
  n_ord <- nrow(p.skel)
  
  # number of drug combinations
  
  n_dose <- ncol(p.skel)
  
  # assigns default, uniform prior for model probabilities if necessary
  
  if(is.na(m_prior)){
    m_prior <- rep(1/n_ord, n_ord)
  }
  
  # indicator for stopping due to safety
  # stop = 1 indicates safety concern
  
  stop <- 0
  
  # initialising vector for integral results for each ordering
  
  marginal_tox <- rep(0, n_ord)
  
  for(k in 1:n_ord){
    
    # obtaining normalising constant for posterior distributions
    
    marginal_tox[k] <- integrate(f = bcrmh, lower = -Inf, upper = Inf, p = p.skel[k,], y = y, n = n, abs.tol = 0)$value
    
  }
  
  # posterior model probabilities
  
  m_post <- (marginal_tox*m_prior)/sum(marginal_tox*m_prior)

  # posterior density of a under a specific ordering
  
  post <- function(a, model){
    return(bcrmh(a = a, p = p.skel[model, ], y = y, n = n)/marginal_tox[model])
  }
  
  # working model for a particular dose, a, and model
  
  wm <- function(dose, a, model){
    return(p.skel[model, dose]**exp(a))
  }
  
  # product of working model and posterior density of a 
  
  prod_post <- function(dose, a, model){
    return(wm(dose, a, model)*post(a, model))
  }
  
  # posterior mean of the working model \Psi_m(d_i,a) across the density of a
  
  marginalised_ptox <- function(dose, model){
    return(integrate(f = prod_post, lower = -Inf, upper = Inf, dose = dose, model = model, abs.tol = 0)$value)
  }
  
  # Applies Bayesian model averaging to combine posterior means for working model
  # under all orderings
  
  bma_ptox <- function(dose){
    bma <- 0
    for(k in 1:n_ord){
      
      # multiplies each posterior mean by corresponding posterior model probability
      
      bma <- bma + m_post[k] * marginalised_ptox(dose = dose, model = k)
    }
    return(bma)
  }
  
  # initialises vector for estimates of probability of toxicity 
  
  ptox_hat <- numeric(ncol(p.skel))
  
  for(i in 1:ncol(p.skel)){
    
    # obtains the estimated probability of toxicity for each dose
    
    ptox_hat[i] <- bma_ptox(i)
    
  }
  
  # checks if overdose controlling is applied
  
  if(is.na(c_od)){
    
    # checks if dose allocation is conservative
    
    if(conserv){
      
      # allocates dose with estimated probability of toxicity closest to TTL which 
      # is strictly below the TTL
      
      distance <- ptox_hat - ttl
      distance[which(distance >= 0)] <- -10000
      comb_curr <- which.is.max(distance)

    } else {
      
      # allocates dose with estimated probability of toxicity closest to TTL
      
      distance <- abs(ptox_hat - ttl)
      comb_curr <- which.is.max(-distance)

    }
    
  } else { 
    
    # initialises vector for probabilities of overdose
    
    prob_overdose <- numeric(n_ord)
    
    # iterates through each dose
    
    for(k in 1:ncol(p.skel)){
      
      # iterates through each potential ordering
      
      for(m in 1:length(m_prior)){
        a_star <- log(log(ttl)/log(p.skel[m,k]))
        
        # obtains BMA estimates for probability of overdose for each dose
        
        prob_overdose[k] <- prob_overdose[k] + m_post[m]*integrate(f = post, model = m, lower = -Inf, upper = a_star, abs.tol = 0)$value
      }
    }
    
    # checks if dose allocation is conservative
    
    if(conserv){
      
      # selects dose closest to TTL which satisfies the overdose control and is 
      # lower than the TTL
      
      distance <- ptox_hat - ttl
      distance[which(prob_overdose >= c_od | distance >= 0)] <- -10000
      comb_curr <- which.is.max(distance)

    } else {
      
      # selects dose closest to TTL which satisfies the overdose control
      
      distance <- abs(ptox_hat - ttl)
      distance[which(prob_overdose >= c_od)] <- 10000
      comb_curr <- which.is.max(-distance)

    }
    
    # checks if selected dose satisfies the overdose control
    # if it does not, this implies no dose level satisfies the overdose control
    # the trial is stopped for safety concerns
    
    if(distance[comb_curr] == 10000){
      stop <- 1
    }
    
  }
  
  return(list(comb_curr = comb_curr,
              m_post = m_post,
              ptox_hat = ptox_hat,
              stop = stop))
  
}

#' Uses Bayesian model averaging in the POCRM to recommend next dose using densities and
#' based on available data
#' 
#' @param p.skel matrix of skeleton probabilities for all possible orderings
#' @param ttl target toxicity level
#' @param y number of DLTs at each dose
#' @param n number of patients treated at each dose
#' @param m_prior prior model probabilities, if unspecified uniform probabilities are assigned
#' @param c_od constant for overdose controlling, if unspecified overdose controlling not applied
#' 
#' @return list containing the recommended next dose, posterior model probabilities,
#' indicator for stopping due to safety concerns, estimated probability of toxicities for each dose
#' 
bma_pocrm_post.imp <- function(p.skel, ttl, y, n, m_prior = NA,  c_od = NA){
  
  # if single ordering is given, convert vector to matrix
  
  if(is.vector(p.skel)){ 
    p.skel <- t(as.matrix(p.skel)) 
  }
  
  # number of potential orderings
  
  n_ord <- nrow(p.skel)
  
  # number of drug combinations
  
  n_dose <- ncol(p.skel)
  
  # assigns default, uniform prior for model probabilities if necessary
  
  if(is.na(m_prior)){
    m_prior <- rep(1/n_ord, n_ord)
  }
  
  # indicator for stopping due to safety
  # stop = 1 indicates safety concern
  
  stop <- 0
  
  # initialising vector for integral results for each ordering
  
  marginal_tox <- rep(0, n_ord)
  
  for(k in 1:n_ord){
    
    # obtaining normalising constant for posterior distributions
    
    marginal_tox[k] <- integrate(f = bcrmh, lower = -Inf, upper = Inf, p = p.skel[k,], y = y, n = n, abs.tol = 0)$value
    
  }
  
  # posterior model probabilities
  
  m_post <- (marginal_tox*m_prior)/sum(marginal_tox*m_prior)
  
  # posterior density for a
  
  post <- function(a, model){
    return(bcrmh(a = a, p = p.skel[model,], y = y, n = n)/marginal_tox[model])
  }
  
  # posterior probability Pr[\Psi(d_i, a) <  \theta], probability for overdose controlling
  
  post_prob_od <- function(tox, dose, model){
    
    # a_star is value of a corresponding to particular toxicity under specific dose, model
    
    a_star <- log(log(tox)/log(p.skel[model,dose]))
    
    # initialising vector for integral values
    
    int <- numeric(length(a_star))
    
    # for loop allows for multiple toxicity levels to be evaluated simultaneously 
    
    for(i in 1:length(a_star)){
      int[i] <- integrate(post, lower = a_star[i], upper = Inf, model = model, abs.tol = 0)$value
    }
    
    return(int)
    
  }
  
  # applies Bayesian model averaging to obtain the weighted average for the 
  # probability of overdose, allows for overdose controlling to be implemented
  
  comb_prob_od <- function(tox, dose){
    total <- 0
    for(model in 1:n_ord){
      
      # obtains probability for each ordering and combines uses posterior model probabilities
      
      total <- total + m_post[model]*post_prob_od(dose = dose, model = model, tox = tox)
      
    }
    return(total)
  }
  
  # combined posterior density of a at a particular toxicity
  
  comb_post <- function(tox, dose){
    total <- 0 
    for(model in 1:n_ord){
      
      # obtains value of a which corresponds to toxicity level of interest
      
      a_star <- log(log(tox)/log(p.skel[model,dose]))
      
      # combines the density for a_star using the posterior model probabilities
      
      total <- total + m_post[model]*post(a = a_star, m = model)
      
    }
    return(total)
  }
  
  # combined posterior density of the probability of toxicity
  
  comb_post_tox <- function(tox, dose){
    
    return(1/abs(tox*log(tox)) * comb_post(tox = tox, dose = dose))
  }
  
  # initialising vectors for density at TTL and probability of overdose for each dose level
  
  pdf <- od_control <- numeric(n_dose)

  for(k in 1:n_dose){
    
    # combined posterior density of TTL for specific dose
    
    pdf[k] <- comb_post_tox(dose = k, tox = ttl)
    
    # probability of overdose at TTL for specific dose
    
    od_control[k] <- 1 - comb_prob_od(dose = k, tox = ttl)
  }
  
  # checks if overdose control is applied
  
  if(is.na(c_od)){
    
    # next dose is recommended by maximising density at TTL for each dose
    
    comb_curr <- which.is.max(pdf)

  } else {
    
    # checks which doses do not satisfy overdose control

    pdf[which(od_control >= c_od)] <- -10000
    
    # selects next dose by maximising density at TTL 
    
    comb_curr <- which.is.max(pdf)
    
    # checks if recommended dose satisfies overdose control
    # if overdose control is not satisfied for this dose, trial is stopped
    
    if(pdf[comb_curr] == -10000){
      stop <- 1
    }
    
  }
  
  return(list(comb_curr = comb_curr,
              m_post = m_post,
              stop = stop,
              ptox_hat = pdf))
}

#' Uses Bayesian model averaging in the POCRM to recommend next dose using estimates
#' of probability of toxicity and based on available data
#' 
#' @param p.skel matrix of skeleton probabilities for all possible orderings
#' @param ttl target toxicity level
#' @param y number of DLTs at each dose
#' @param n number of patients treated at each dose
#' @param m_prior prior model probabilities, if unspecified uniform probabilities are assigned
#' @param c_od constant for overdose controlling, if unspecified overdose controlling not applied
#' 
#' @return list containing the recommended next dose, posterior model probabilities,
#' indicator for stopping due to safety concerns, estimated probability of toxicities for each dose
#' 
bma_pocrm_post2.imp <- function(p.skel, ttl, y, n, m_prior = NA,  c_od = NA){
  
  # if single ordering is given, convert vector to matrix
  
  if(is.vector(p.skel)){ 
    p.skel <- t(as.matrix(p.skel)) 
  }
  
  # number of potential orderings
  
  n_ord <- nrow(p.skel)
  
  # number of drug combinations
  
  n_dose <- ncol(p.skel)
  
  # assigns default, uniform prior for model probabilities if necessary
  
  if(is.na(m_prior)){
    m_prior <- rep(1/n_ord, n_ord)
  }
  
  # indicator for stopping due to safety
  # stop = 1 indicates safety concern
  
  stop <- 0
  
  # initialising vector for integral results for each ordering
  
  marginal_tox <- rep(0, n_ord)
  
  for(k in 1:n_ord){
    
    # obtaining normalising constant for posterior distributions
    
    marginal_tox[k] <- integrate(f = bcrmh, lower = -Inf, upper = Inf, p = p.skel[k,], y = y, n = n, abs.tol = 0)$value
    
  }
  
  # posterior model probabilities
  
  m_post <- (marginal_tox*m_prior)/sum(marginal_tox*m_prior)
  
  # posterior density of a
  
  post_a <- function(a, model){
    return(bcrmh(a = a, p = p.skel[model,], y = y, n = n)/marginal_tox[model])
  }
  
  # posterior density of the probability of toxicity under a specific dose and ordering
  
  post_tox <- function(tox, dose, model){
    return(1/abs(tox*log(tox)) * post_a(a = log(log(tox)/log(p.skel[model,dose])), model = model))
  }
  
  # combined posterior density of the probability of toxicity under a specific dose and ordering
  
  comb_post_tox <- function(tox, dose){
    total <- 0
    for(m in 1:n_ord){
      total <- total + m_post[m] * post_tox(tox = tox, dose = dose, model = m)
    }
    return(total)
  }
  
  # function to be integrated to obtain the expected value of the probability of toxicity
  # for a specific dose
  
  integrand <- function(tox, dose){
    return(tox*comb_post_tox(tox = tox, dose = dose))
  }
  
  # posterior probability Pr[\Psi(d_i, a) <  \theta], probability for overdose controlling
  
  post_prob_od <- function(tox, dose, model){
    
    # a_star is value of a corresponding to particular toxicity under specific dose, model
    
    a_star <- log(log(tox)/log(p.skel[model,dose]))
    
    # initialising vector for integral values
    
    int <- numeric(length(a_star))
    
    # for loop allows for multiple toxicity levels to be evaluated simultaneously 
    
    for(i in 1:length(a_star)){
      int[i] <- integrate(post_a, lower = a_star[i], upper = Inf, model = model, abs.tol = 0)$value
    }
    
    return(int)
    
  }
  
  # applies Bayesian model averaging to obtain the weighted average for the 
  # probability of overdose, allows for overdose controlling to be implemented
  
  comb_prob_od <- function(tox, dose){
    total <- 0
    for(model in 1:n_ord){
      
      # obtains probability for each ordering and combines uses posterior model probabilities
      
      total <- total + m_post[model]*post_prob_od(dose = dose, model = model, tox = tox)
      
    }
    return(total)
  }
  
  exp <- numeric(n_dose)
  
  for(k in 1:n_dose){
    
    # integrating across possible values of probability of toxicity to determine expected value
    
    exp[k] <- integrate(integrand, lower = 0, upper = 1, dose = k, abs.tol = 0)$value
  }

  # checks if overdose controlling is applied
    
  if(is.na(c_od)){
    
    # next dose is assigned by minimising difference between expected value of 
    # probability of toxicity and TTL
    
    distance <- abs(exp - ttl)
    comb_curr <- which.is.max(-distance)
  } else {
    
    prob_od <- numeric(n_dose)
    
    for(k in 1:n_dose){
      
      prob_od[k] <- 1 - comb_prob_od(tox = ttl, dose = k)
      
    }
    
    distance <- abs(exp - ttl)
    distance[which(prob_od > c_od)] <- 10000
    comb_curr <- which.is.max(-distance)
    
    if(distance[comb_curr] == 10000){
      stop <- 1
    }
    
  }
  
  return(list(comb_curr = comb_curr,
              m_post = m_post,
              ptox_hat = exp,
              stop = stop))
}

#' Uses Bayesian model averaging in the POCRM to recommend next dose using estimates
#' of probability of toxicity and based on available data
#' 
#' @param p.skel matrix of skeleton probabilities for all possible orderings
#' @param ttl target toxicity level
#' @param y number of DLTs at each dose
#' @param n number of patients treated at each dose
#' @param m_prior prior model probabilities, if unspecified uniform probabilities are assigned
#' @param c_od constant for overdose controlling, if unspecified overdose controlling not applied
#' @param num_rand number of doses for which randomisation should be applied, rounds down from given number
#' 
#' @return list containing the recommended next dose, posterior model probabilities,
#' indicator for stopping due to safety concerns, estimated probability of toxicities for each dose
#' 
bma_pocrm_rand.imp <- function(p.skel, ttl, y, n, m_prior = NA,  c_od = NA, num_rand = NA){
  
  # if single ordering is given, convert vector to matrix
  
  if(is.vector(p.skel)){ 
    p.skel <- t(as.matrix(p.skel)) 
  }
  
  # ensuring that num_rand is an integer
  
  num_rand <- floor(num_rand)
  
  # number of potential orderings
  
  n_ord <- nrow(p.skel)
  
  # number of drug combinations
  
  n_dose <- ncol(p.skel)
  
  if (is.na(num_rand)) {
  
    num_rand <- n_dose
    
  } else if (num_rand < 2){
    
    warning("Randomisation must occur between at least 2 doses. num_rand has defaulted to 2.")
    num_rand <- 2
    
  }
  
  # assigns default, uniform prior for model probabilities if necessary
  
  if(is.na(m_prior)){
    m_prior <- rep(1/n_ord, n_ord)
  }
  
  # indicator for stopping due to safety
  # stop = 1 indicates safety concern
  
  stop <- 0
  
  # initialising vector for integral results for each ordering
  
  marginal_tox <- rep(0, n_ord)
  
  for(k in 1:n_ord){
    
    # obtaining normalising constant for posterior distributions
    
    marginal_tox[k] <- integrate(f = bcrmh, lower = -Inf, upper = Inf, p = p.skel[k,], y = y, n = n, abs.tol = 0)$value
    
  }
  
  # posterior model probabilities
  
  m_post <- (marginal_tox*m_prior)/sum(marginal_tox*m_prior)
  
  m_sel <- sample.int(n = num_rand, size = 1, prob = m_post[1:num_rand])
  
  # posterior density of a
  
  post_a <- function(a, model){
    return(bcrmh(a = a, p = p.skel[model,], y = y, n = n)/marginal_tox[model])
  }
  
  # posterior density of the probability of toxicity under a specific dose and ordering
  
  post_tox <- function(tox, dose, model){
    return(1/abs(tox*log(tox)) * post_a(a = log(log(tox)/log(p.skel[model,dose])), model = model))
  }
  
  # combined posterior density of the probability of toxicity under a specific dose and ordering
  
  # function to be integrated to obtain the expected value of the probability of toxicity
  # for a specific dose
  
  integrand <- function(tox, dose){
    return(tox*post_tox(tox = tox, dose = dose, model = m_sel))
  }
  
  # posterior probability Pr[\Psi(d_i, a) <  \theta], probability for overdose controlling
  
  post_prob_od <- function(tox, dose, model){
    
    # a_star is value of a corresponding to particular toxicity under specific dose, model
    
    a_star <- log(log(tox)/log(p.skel[model,dose]))
    
    # initialising vector for integral values
    
    int <- numeric(length(a_star))
    
    # for loop allows for multiple toxicity levels to be evaluated simultaneously 
    
    for(i in 1:length(a_star)){
      int[i] <- integrate(post_a, lower = a_star[i], upper = Inf, model = model, abs.tol = 0)$value
    }
    
    return(int)
    
  }
  
  # applies Bayesian model averaging to obtain the weighted average for the 
  # probability of overdose, allows for overdose controlling to be implemented
  
  exp <- numeric(n_dose)
  
  for(k in 1:n_dose){
    
    # integrating across possible values of probability of toxicity to determine expected value
    
    exp[k] <- integrate(integrand, lower = 0, upper = 1, dose = k, abs.tol = 0)$value
  }
  
  # checks if overdose controlling is applied
  
  if(is.na(c_od)){
    
    # next dose is assigned by minimising difference between expected value of 
    # probability of toxicity and TTL
    
    distance <- abs(exp - ttl)
    comb_curr <- which.is.max(-distance)
    
  } else {
    
    prob_od <- numeric(n_dose)
    
    for(k in 1:n_dose){
      
      prob_od[k] <- 1 - post_prob_od(tox = ttl, dose = k, model = m_sel)
      
    }
    
    distance <- abs(exp - ttl)
    distance[which(prob_od > c_od)] <- 10000
    comb_curr <- which.is.max(-distance)
    
    if(distance[comb_curr] == 10000){
      stop <- 1
    }
    
  }
  
  return(list(comb_curr = comb_curr,
              m_post = m_post,
              ptox_hat = exp,
              stop = stop))
}

#' Runs a single trial simulation
#' 
#' @param p0 true toxicity probabilities
#' @param p.skel skeleton probabilities
#' @param ttl target toxicity level
#' @param cohortsize size of each cohort (either number or vector)
#' @param ncohort number of cohorts
#' @param n_stop number of patients assigned to single combination to stop trial
#' @param no_skip boolean which determines if dose skipping is allowed at allocation
#' @param start_comb starting combination
#' @param c_od constant used for overdose controlling
#' @param cs confidence level used by POCRM method
#' @param method specifies which method is used for recommendation from "Wages", "BMA Mixture", "BMA Point", "BMA Mixture2"
#' @param m_prior  prior model probabilities, if unspecified uniform probabilities are assigned
#' @param conserv boolean for whether dose allocation is conservative:
#' @param sim_method specifies which method is used for data generation from "Default", "Mozgunov"
#' @param num_rand number of doses for which randomisation should be carried out, minimum of 2
#' 
#' @return list containing the recommended dose, toxicity data, patient treatment data, 
#' the safety stopping indicator, posterior model probabilities and historical data describing the full trial
#' 
bma_pocrm <- function(p0, p.skel, ttl, cohortsize, ncohort, n_stop, no_skip = F, start_comb, c_od = NA, cs = NA, method = "BMA Mixture", m_prior = NA, conserv = F, sim_method = "Default", num_rand = NA){
  
  # if single ordering is passed, convert to matrix
  
  if(is.vector(p.skel)){
    p.skel <- as.matrix(p.skel)
  }
  
  # check that format of cohort size is correct
  
  if(length(cohortsize) != 1 && length(cohortsize) != ncohort){
    stop("Cohort sizes must be constant or a cohort size must be specified for each cohort.")
  }
  
  # number of drug combinations 
  
  ncomb <- ncol(p.skel)
  
  # initialising vectors for patient data, estimated probability of toxicities, and combination selections
  
  y <- n <- ptox_hat <- comb_select <- comb_select_hist <- hist <- trial <- numeric(ncomb)
  
  # initialising starting combination
  
  comb_curr <- start_comb

  # safety stopping indicator
  # stop = 1 indicates safety concern, trial is stopped
  
  stop <- 0

  # generate matrices for storing historical posterior model probability and estimates
  # for the probability of toxicity at each step of trial
  
  m_post_hist <- matrix(0, nrow = ncohort, ncol = nrow(p.skel))
  ptox_hat_hist <- matrix(0, nrow = ncohort, ncol = ncol(p.skel))
  
  # formatting of cohortsize if single size is specified 
  
  if(length(cohortsize) == 1){
    cohortsize <- rep(cohortsize, ncohort)
  }
  
  # generates tolerance level for each patient depending on simulation method
  
  if(sim_method == "Mozgunov"){
    patient_tol <- runif(ncohort*cohortsize)
  }
  
  i <- 1
  
  while(i <= ncohort){
    
    # generates DLT data for current trial run, depending on method specified
    
    # default method: uses draws from binomial distribution with p = true probability of 
    # toxicity for current dose and n = cohortsize for current cohort
    
    # mozgunov method: compares tolerance for each patient in cohort to true probability
    # of toxicity for current dose, counts number of DLTs observed
    
    trial[i] <- switch(sim_method,
                       "Default" = rbinom(1, cohortsize[i], p0[comb_curr]),
                       "Mozgunov" = sum(p0[comb_curr]  > patient_tol[(sum(n)+1):(sum(n)+cohortsize[i])]))
    
    # iterate DLT and patient treatment data under each dose
    
    y[comb_curr] <- y[comb_curr] + trial[i]
    n[comb_curr] <- n[comb_curr] + cohortsize[i]

    # checks if maximum number of patients under one dose is reached
    # stops trial for non-safety concerns if satisfied
    
    if(any(n > n_stop)){
      stop <- 0
      break
    }
    
    # obtains dose recommendation based on selected method
    
    results <- switch(method,
                      "BMA Mixture" = bma_pocrm_post.imp(p.skel = p.skel, ttl = ttl, y = y, n = n, m_prior = m_prior, c_od = c_od),
                      "BMA Point" = bma_pocrm_point.imp(p.skel = p.skel, ttl = ttl, y = y, n = n, m_prior = m_prior, c_od = c_od, conserv = conserv),
                      "Wages" = pocrm.imp(p.skel = p.skel, ttl = ttl, y = y, n = n, cs = cs, m_prior = m_prior, c_od = c_od, conserv = conserv),
                      "BMA Mixture2" = bma_pocrm_post2.imp(p.skel = p.skel, ttl = ttl, y = y, n = n, m_prior = m_prior,  c_od = c_od),
                      "Randomised" = bma_pocrm_rand.imp(p.skel = p.skel, ttl = ttl, y = y, n = n, m_prior = m_prior, c_od = c_od, num_rand = num_rand))
    
    comb_best <- results$comb_curr
    m_post <- results$m_post
    stop <- results$stop
    ptox_hat_hist[i,] <- results$ptox_hat
    
    m_post_hist[i,] <- m_post
    
    # checks if dose skipping is allowed
    
    if(no_skip){
      # chooses lower dose between the next highest dose, and the selected dose
      # prevents dose skipping
      
      comb_curr <- min(comb_best, comb_curr + 1)
    } else {
      comb_curr <- comb_best
    }
    
    # checks if trial is stopped for safety
    
    if(stop == 1){
      break
    } else {
      
      # historical combination selection data updated
      
      comb_select_hist[i] <- comb_curr
    }
    
    i = i + 1
  }

  comb_select[comb_curr] <- comb_select[comb_curr] + 1
  
  hist <- cbind(cohortsize, trial, m_post_hist, ptox_hat_hist, comb_select_hist)
  colnames(hist) <- c("n", "y", paste0("m",1:nrow(p.skel)), paste0("d",1:ncol(p.skel)), "comb_select")
  
  output <- list(comb_select = comb_select,
                 tox_data = y,
                 n_patients = n,
                 stop = stop,
                 m_post = m_post,
                 hist = hist)
  
  return(output)
}

#' Runs repeated trial simulation
#' 
#' @param p0 true toxicity probabilities
#' @param p.skel skeleton probabilities
#' @param ttl target toxicity level
#' @param cohortsize size of each cohort (either number or vector)
#' @param ncohort number of cohorts
#' @param n_stop number of patients assigned to single combination to stop trial
#' @param no_skip boolean which determines if dose skipping is allowed at allocation
#' @param start_comb starting combination
#' @param c_od constant used for overdose controlling
#' @param cs confidence level used by POCRM method
#' @param method specifies which method is used for recommendation from "Wages", "BMA Mixture", "BMA Point", "BMA Mixture2"
#' @param m_prior  prior model probabilities, if unspecified uniform probabilities are assigned
#' @param conserv boolean for whether dose allocation is conservative:
#' @param output specifies whether output is printed or returned in a list
#' @param sim_method specifies which method is used for data generation from "Default", "Mozgunov"
#' @param num_rand number of doses for which randomisation should be applied, minimum of 2
#' 
#' @return list containing the recommended dose, toxicity data, patient treatment data, 
#' the safety stopping indicator, posterior model probabilities and historical data describing the full trial
#' 
bma_pocrm.sim <- function(p0, p.skel, ttl, start_comb, cohortsize, ncohort, n_stop, ntrial, c_od = NA, no_skip = F, cs = 0.9, method = "BMA Mixture", m_prior = NA, conserv = F, output = "print", sim_method = "Default", num_rand = NA){
  
  # check that valid allocation method is selected
  
  if(method != "BMA Mixture" && method != "BMA Point" && method != "Wages" && method != "BMA Mixture2"){
    stop("Invalid method specified. The options are BMA Mixture, BMA Point, and Wages.")
  }
  
  ncomb <- length(p0)
  
  comb_select <- y <- n <- matrix(nrow = ntrial, ncol = ncomb)
  m_post <- matrix(nrow = ntrial, ncol = nrow(p.skel))
  nstop = 0
  
  # contains appendices of trials that have been stopped 
  trial_stop <- c()
  
  hist <- matrix(0, nrow = ncohort, ncol = (nrow(p.skel)+ncol(p.skel)+3))

  for (i in 1:ntrial){
    result <- bma_pocrm(p0 = p0, 
                        p.skel = p.skel, 
                        cohortsize = cohortsize, 
                        ttl = ttl, 
                        ncohort = ncohort,
                        n_stop = n_stop,
                        start_comb = start_comb,
                        cs = cs,
                        c_od = c_od,
                        no_skip = no_skip,
                        method = method, 
                        m_prior = m_prior,
                        conserv = conserv,
                        sim_method = sim_method,
                        num_rand = num_rand)
    
    comb_select[i,] <- result$comb_select
    hist <- hist + result$hist
    if(result$stop == 1){
      trial_stop <- append(trial_stop, i)
    }
    nstop <- nstop + result$stop
    y[i,] <- result$tox_data
    n[i,] <- result$n_patients
    m_post[i,] <- result$m_post
  }
  
  # removes recommendations based on trials that have been stopped
  if(length(trial_stop) != 0){
    comb_select <- comb_select[-trial_stop,]
  }
  
  if(output == "list"){
    output <- list(true_tox = round(p0,4),
                   select_percent = colMeans(comb_select)*100,
                   num_dlt = colMeans(y)/(ncohort*cohortsize)*100,
                   num_pat = colMeans(n)/(ncohort*cohortsize)*100,
                   percent_stop = nstop/ntrial*100,
                   m_post = colMeans(m_post))
                  #hist = round(hist/ntrial,4)
        return(output)
  } else {
    cat("True tox probability: ", round(p0,3), sep="    ",  "\n");
    cat("selection percentage: ", formatC(colMeans(comb_select)*100, digits=1, format="f"), sep="    ",  "\n");
    cat("number of toxicities:    ", formatC(colMeans(y)/(cohortsize*ncohort)*100, digits=1, format="f"), sep="    ",   "\n");
    cat("number of patients treated:     ", formatC(colMeans(n)/(cohortsize*ncohort)*100, digits=1, format="f"), sep="    ",   "\n");
    cat("percentage of stop:    ", nstop/ntrial*100, "\n");
    cat("posterior model probabilities: ", formatC(colMeans(m_post), digits = 4, format = "f"), sep = "    ",   "\n")
  }
}

#' checks if the value in the index of a particular vector is unique in the vector
#'
check_rand <- function(index, vector, location){
  if(vector[index] %in% vector[duplicated(vector)]){
    print(location)
    return(T)
  }
  return(F)
}
