#JAGS code to run the analysis from the barred owl application.

model {
  #Priors
  lambda ~ dunif(0,10) # initial abundance
  p.occ ~ dunif(0, 1) # detection
  b0 ~ dnorm(0,0.1) # intercept on survival
  b1 ~ dnorm(0,0.1) # slope on survival
  a0 ~ dnorm(0,0.1) # intercept on effort (p.count)
  a1 ~ dnorm(0,0.1) # slope on effort (p.count)
  g0 ~ dnorm(0,0.1) # intercept on gamma
  g1 ~ dnorm(0,0.1) # slope on mean(N)
  g2 ~ dnorm(0,0.1) # squared term on mean(N)
  
  #Likelihood - Biological process model
  for(i in 1:nSites) {
    
    #First year of sampling - process and observation components
    N[i,1] ~ dpois(lambda)
    
    #All other years of sampling - process and observation components   
    for(t in 2:nYears) {
      logit(omega[i,t-1]) <- b0 + b1*PHAB.std[i,t-1]
      S[i,t-1] ~ dbin(omega[i,t-1], N[i,t-1])
      G[i,t-1] ~ dpois(gamma[t-1])
      N[i,t] <- S[i,t-1] + G[i,t-1] 
    }}
  
  #Detection model for occupancy data
  for (k in 1:occ.end){
    p.site[k] <- 1-pow((1-p.occ),N[site[k],year[k]]) 
    y[k] ~ dbern(p.site[k])                           #y is a vector of each observation
    p.count[k]<-0
  }
  #Detection for the count model 
  for (k in count.start:count.end){
    logit(p.count[k]) <- a0 + a1*OV[k]
    y[k] ~ dbin(p.count[k],N[site[k],year[k]])
  }
  
  # covariate for gamma
  for (t in 2:nYears){
    log(gamma[t-1]) <- g0 + g1*N.mean[t-1] + g2*N.mean[t-1]*N.mean[t-1]  # covariates on gains 
    N.mean[t-1] <- mean(N[,t-1]) - 1   # per site mean from previous year
  }
#Close the model file.  
}
