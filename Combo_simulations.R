#######################################################################################
#######################################################################################
#Scenario 1: Basic version of the model combining count and detection/nondetection data
#######################################################################################
#######################################################################################


#########################################################################
# Data Generation
########################################################################

# Simulate data on N under the so-called "constant" model
# Basic birth-death process but birth rate isn't affected by
# abundance in previous year. 

#"True" values
lambda <- runif(1,0.5,3)
omega <- runif(1,0,1)
gamma <- runif(1,0,2.5)
p <- runif(1,0,1)

#Specify the number of sites, years, and reps
nYears <- 10
nReps <- 3
nCount<-sample(c(30, 15, 5),1)           # number of sites with count data
nOcc<-sample(c(25, 75, 150),1)        # number of sites with 
#detection/nondetection data
nSites<-nCount + nOcc                    # total sites

#Simulate true abundances, N, for each location
N <- matrix(NA, nSites, nYears)
S <- G <- matrix(NA, nSites, nYears-1)

#First year of sampling follows a Pois distribution
N[,1] <- rpois(nSites, lambda)

#Subsequent years follow the birth-death-immigration process
for(t in 2:nYears) {
  S[,t-1] <- rbinom(nSites, N[,t-1], omega)
  G[,t-1] <- rpois(nSites, gamma)
  N[,t] <- S[,t-1] + G[,t-1] 
}

#Generate data vector y for the counts
y <- array(NA, c(nSites, nYears, nReps))
for(t in 1:nYears) {
  for(j in 1:nReps) {
    y[,t,j] <- rbinom(nSites, N[,t], p)
  }
}

#Assume that there are a vector of sites, x, have only have #detection/nondetection data
#And change their data to detection/nondetection data
x=1:nOcc
for (i in x) {
  for (j in 1:nReps) {
    a = which(y[i,,j]>0)
    y[i,a,j] = 1
  }
}

#Divide the data into two datasets
#y1 = detection/nondetection data
y1 = y[1:nOcc,,]

#y2 = count data
y2 = y[(nOcc+1):nSites,,]

#########################################################################
# Create the JAGS file
########################################################################

sink("combo_model.R")
cat("
    model {
    #Priors
    lambda ~ dunif(0,10) # initial abundance
    gamma ~  dunif(0,10) # gains
    omega ~ dunif(0, 1)  # survival
    p ~ dunif(0, 1)      # detection
    
    #Likelihood - Biological process model
    for(i in 1:nSites) {
    #First year of sampling 
    N[i,1] ~ dpois(lambda)
    
    #All other years of sampling 
    for(t in 2:nYears) {
    S[i,t-1] ~ dbin(omega, N[i,t-1])
    G[i,t-1] ~ dpois(gamma)
    N[i,t] <- S[i,t-1] + G[i,t-1]
    }
    }
    
    #Detection process model for detection/nondetection data
    for (i in 1:nOcc) {
    for (t in 1:nYears) {
    p.site[i,t] <- 1-pow( (1-p),N[i,t] )
    for (j in 1:nReps) {  
    y1[i,t,j] ~ dbern(p.site[i,t])
    }}}
    
    #Detection process model for count data
    for (i in 1:nCount) {
    for (j in 1:nReps) {  
    for (t in 1:nYears) {
    y2[i,t,j] ~ dbin(p, N[i+nOcc,t])
    }}}
    }
    ",fill = TRUE)
sink()

#########################################################################
# Run the JAGS code
########################################################################

#Format data 
jags.data <- list(nSites=nSites, 
                  nOcc=nOcc, 
                  nCount=nCount, 
                  nYears=nYears, 
                  y1=y1,y2=y2,
                  nReps=nReps)

#Initial values
#Note, JAGS will throw an error if the initial values aren't in agreement
#with the data. It helps to start N at large values. 

#Parameters monitored
params <- c("lambda", "gamma", "omega", "p", "N")

#Generate inits
Ni <- y[,,1]+20
Si <- S
Si[] <- 2
Gi<-matrix(10,nrow=nSites,ncol=(nYears-1))
Ni[,-1] <- NA

#Path to model file
model=normalizePath("combo_model.R")

#Load the correct library
#library("jagsUI")
library(jagsUI)
#Compile the model, and ensure correct inits are found (may take multiple tries)
inits <- function() list(N=Ni,
                         S=Si,
                         G=Gi)
#Run the model
jags.out<-jags(jags.data, inits, params, model.file=model, store.data=TRUE,
               n.chains=3, n.iter=10000, n.burnin=2000, n.thin=10, 
               n.adapt=200)




#############################################################
#############################################################
#Scenario 2: Model containing a covariate effect on survival
#############################################################
#############################################################


#########################################################################
# Data Generation
########################################################################

# Simulate data for N under the covariate model

#"True" values
lambda <- 4
gamma <- 2
p <- 0.5
b0 <- 0.5
b1 <- 0.7

#Specify the number of sites, years, and reps
nYears <- 10
nReps <- 3
nCount<-40
nOcc<-100

#Generate covariate values across the four sampling scenarios
scenario<-sample(1:4,1)
if (scenario==1){
  occ.covar <- runif(nOcc,-3,3)
  count.covar <- runif(nCount,-3,3)
}
if (scenario==2){
  occ.covar <- runif(nOcc,-3,1)
  count.covar <- runif(nCount,1,3)
}
if (scenario==3){
  occ.covar <- sample(c(runif(10000, -3,-2),runif(10000, 2,3)),nOcc)
  count.covar <- runif(nCount,-2,2)
}
if (scenario==4){
  occ.covar <- runif(nOcc,1,3)
  count.covar <- runif(nCount,1,3)
}

covariate<-c(occ.covar,count.covar)
nSites <- nOcc + nCount

#Calculate survival probabilities for each site
omega<-plogis(b0 + b1*covariate)

#Simulate true abundances, N, for each location
N <- matrix(NA, nSites, nYears)
S <- G <- matrix(NA, nSites, nYears-1)

#First year of sampling follows a Pois distribution
N[,1] <- rpois(nSites, lambda)

#Subsequent years follow the birth-death-immigration process
for(t in 2:nYears) {
  for (i in 1:nSites) {
    S[i,t-1] <- rbinom(1, N[i,t-1], omega[i])
    G[,t-1] <- rpois(nSites, gamma)
    N[,t] <- S[,t-1] + G[,t-1]
  }}

#Generate data vector y for the counts
y <- array(NA, c(nSites, nYears, nReps))
for(t in 1:nYears) {
  for(j in 1:nReps) {
    y[,t,j] <- rbinom(nSites, N[,t], p)
  }
}

#The data, y, will be converted to a vector from 1:nSamples. 
#Year and site for each observation will be provided in an accompanying 
#vector "site" or "year"

year<-array(NA, dim=dim(y))
site<-array(NA, dim=dim(y))
for (i in 1:nSites){
  for (t in 1:nYears){
    year[i,t,]<-t
    site[i,t,]<-i
  }
}

#Convert arrays to vectors
y<-c(y)
site<-c(site)
year<-c(year)
nSamples<-length(y)

#Separate detection/nondetection data from count
occ.samples<-c()
count.samples<-c()
for (i in 1:nSamples){
  if(site[i]<=nOcc){                 # determine if site falls within 
    #occupancy site
    occ.samples<-c(occ.samples, i)   # instances from 1:nSamples that will be analyzed as occ data
    if(y[i] > 0){y[i]=1}             # convert counts to detection/non-detection
  }else {
    count.samples<-c(count.samples,i)  # instances from 1:nSamples that will be analyzed as count data
  }}

#########################################################################
# Create the JAGS file
########################################################################

sink("combo_covar_model.R")
cat("
    model {
    #Priors
    lambda ~ dunif(0,10)
    gamma ~  dunif(0,10)
    p ~ dunif(0, 1)
    b0 ~ dnorm(0,0.1)
    b1 ~ dnorm(0,0.1)
    
    #Likelihood - Biological process model
    for(i in 1:nSites) {
    #First year of sampling 
    N[i,1] ~ dpois(lambda)
    logit(omega[i]) <- b0 + b1*covariate[i]
    
    #All other years of sampling 
    for(t in 2:nYears) {
    S[i,t-1] ~ dbin(omega[i], N[i,t-1])
    G[i,t-1] ~ dpois(gamma)
    N[i,t] <- S[i,t-1] + G[i,t-1] 
    }}
    
    #Detection model for detection/nondetection data
    for (k in 1:length.occ.samples){
    occ.p[k] <- 1-pow( (1-p),N[site[occ.samples[k]],year[occ.samples[k]]] )
    y[occ.samples[k]] ~ dbern(occ.p[k])
    }
    
    #Detection model for count data
    for (k in 1:length.count.samples){
    y[count.samples[k]] ~ 
    dbin(p,N[site[count.samples[k]],year[count.samples[k]]])
    }
    }
    ",fill = TRUE)
sink()

#########################################################################
# Run the JAGS code
########################################################################

#Format data 
jags.data <- list(nSites=nSites, 
                  nYears=nYears,
                  length.count.samples=length(count.samples),
                  length.occ.samples=length(occ.samples),
                  y=y,
                  site=site,
                  year=year,
                  covariate=covariate,
                  occ.samples=occ.samples,
                  count.samples=count.samples)

#Parameters to monitor
params<-c("b1","b0","gamma","p","lambda","N")

#Path to model file
model<- normalizePath("combo_covar_model.R")

#Generate inits
Ni <- N[,]
Si <- S
Si[] <- 1
Gi<-matrix(1,nrow=nSites,ncol=(nYears-1))
Ni[,-1] <- NA

#Run the model
library(jagsUI)

#Finding the right inits can be difficult, this while loop was designed to 
#automate the process

#If after 15 attempts, no suitable inits were generated change the values on 
#Ni, Si and Gi

jags.out<-NA
class(jags.out)<-"try-error"
counter=0
while (class(jags.out)=="try-error"){
  inits <- function() list(N=Ni+rpois(1,20),
                           S=Si+rpois(1,8),
                           G=Gi+rpois(1,8))
  counter=counter+1
  if(counter>15){break}
  
  jags.out<-try(jags(jags.data, inits, params, model.file=model, 
                     store.data=TRUE, n.chains=3, n.iter=10, n.burnin=2, 
                     n.thin=1, n.adapt=0))
}




############################################################
############################################################
#Scenario 3: Model containing a covariate effect on survival
############################################################
############################################################


#########################################################################
# Data Generation
########################################################################

#Simulate data on N under assuming that detection varies by survey type

#"True" values
lam <- 1
omega <- 0.7
gamma <- 1.5
scenario<-sample(1:4,1)

#Generate detection probabilities according to the four specified scenarios
if (scenario %in% c(1,3)){
  p.count <- 0.5
  p.occ <- 0.3
}else{
  p.count <- 0.3
  p.occ <- 0.5
}

#Specify the number of sites, years, and reps
nYears <- 10
nReps <- 3
nOcc <- sample(c(25,75,150),1)    # sites with detection/nondetection data
nCount <- sample(c(5,15,75),1)    # sites with count data
nSites<-nCount + nOcc

#Simulate true abundances, N, for each location
N <- matrix(NA, nSites, nYears)
S <- G <- matrix(NA, nSites, nYears-1)

#First year of sampling follows a Pois distribution
N[,1] <- rpois(nSites, lam)

#Subsequent years follow the birth-death-immigration process
for(t in 2:nYears) {
  S[,t-1] <- rbinom(nSites, N[,t-1], omega)
  G[,t-1] <- rpois(nSites, gamma)
  N[,t] <- S[,t-1] + G[,t-1] 
}

#Generate data vector y for the counts
y <- array(NA, c(nSites, nYears, nReps))
for(t in 1:nYears) {
  for(j in 1:nReps) {
    for (i in 1:nOcc){
      y[i,t,j] <- rbinom(1, N[i,t], p.occ)
    }
    for (i in 1:nCount){
      y[i+nOcc,t,j] <- rbinom(1, N[i+nCount,t], p.count)
    }}}

#Assume that there are a vector of sites, x, have only have #detection/nondetection data
#Change their data to detection/nondetection data
x=1:nOcc
for (i in x) {
  for (j in 1:nReps) {
    a = which(y[i,,j]>0)
    y[i,a,j] = 1
  }
}
#y1 = detection/nondetection data
y1 = y[1:nOcc,,]
y2 = y[(nOcc+1):nSites,,]

#########################################################################
# Create the JAGS file
########################################################################

sink("combo_model_2p.R")
cat("
    model {
    #Priors
    lambda ~ dunif(0,10)
    gamma ~  dunif(0,10)
    omega ~ dunif(0, 1)
    p.occ ~ dunif(0, 1)
    p.count ~ dunif(0, 1)
    
    #Likelihood - Biological process model
    for(i in 1:nSites) {
    #First year of sampling 
    N[i,1] ~ dpois(lambda)
    
    #All other years of sampling 
    for(t in 2:nYears) {
    S[i,t-1] ~ dbin(omega, N[i,t-1])
    G[i,t-1] ~ dpois(gamma)
    N[i,t] <- S[i,t-1] + G[i,t-1]
    }}
    
    #Detection models
    #Detection model for detection/nondetection data
    for (i in 1:nOcc) {
    for (t in 1:nYears) {
    p.site[i,t] <- 1-pow( (1-p.occ),N[i,t] )
    for (j in 1:nReps) {  
    y1[i,t,j] ~ dbern(p.site[i,t])
    }}}
    
    #Detection model for count data
    for (i in 1:nCount) {
    for (j in 1:nReps) {  
    for (t in 1:nYears) {
    y2[i,t,j] ~ dbin(p.count, N[i+nOcc,t])
    }}}
    }
    ",fill = TRUE)
sink()

#########################################################################
# Run the JAGS code
########################################################################

#Format data 
jags.data <- list(nSites=nSites, 
                  nOcc=nOcc, 
                  nCount=nCount, 
                  nYears=nYears, 
                  y1=y1,y2=y2,
                  nReps=nReps)

#Parameters monitored
if (scenario%in%c(1,2)){
  params <- c("lambda", "gamma", "omega", "p.count","p.occ", "N")
}else{
  params <- c("lambda", "gamma", "omega", "p", "N")
}

#Initial values
#Note, JAGS will throw an error if the initial values aren't in agreement
#with the data. It helps to start N at large values
Ni <- y[,,1]+20
Si <- S
Si[] <- 2
Gi<-matrix(10,nrow=nSites,ncol=(nYears-1))
Ni[,-1] <- NA

if(scenario %in% 1:2){
  model=normalizePath("combo_model_2p.R")
}else {
  model=normalizePath("combo_model.R")  #This is the original model from 
  #Scenario 1
}

#Load the correct library
#library("jagsUI")
library(jagsUI)

#Compile the model, and insure correct inits are found 
#(may take multiple tries)

inits <- function() list(N=Ni,
                         S=Si,
                         G=Gi)

jags.out<-jags(jags.data, inits, params, model.file=model, store.data=TRUE,
               n.chains=3, n.iter=10, n.burnin=2, n.thin=1, n.adapt=0)
