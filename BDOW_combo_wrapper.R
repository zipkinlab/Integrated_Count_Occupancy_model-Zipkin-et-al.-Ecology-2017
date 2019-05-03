rm(list=ls())

#Note: the results for the analysis in Zipkin et al. (2017) are mistakenly 
#reported for a model in which the detection model for the count data was:
#plogis(p.count[k]) <- a0 + a1*OV[k]
#when it should have been:
# cloglog(p.count[k]) <- a0 + log(OV[k])
#as specified in the paper.

#This code contains the corrected model.  The only substantial difference
#is in left handside of Figure 4e, which shows the relationship between
#percent of area surveyed and the detection probability.

#We included an updated pdf of Figure 4e in the Github repo for reference.


##########################################

#Below is the code to read in the data, reshape it, and run the JAGS model.

########################################################################
# Reshape data
#######################################################################
# Read in data 
model.file<-file.choose() # select BDOW_combo_model_cloglog.R

#set working directory (all files must be here)
setwd(dirname(model.file))

# import data
raw.data<-read.csv("BDOW_combo_data.csv",header=TRUE) 

# Replace unsampled entries with NA
raw.data[raw.data=="."]<-NA

# Barred Owl data
BO<-as.matrix(raw.data[,grep("BO",colnames(raw.data))])

# Years
year<-as.numeric(gsub("[^0-9]","",gsub("_.*$","",colnames(BO))))
# Adjust so first year = 1
year[year<90]<-year[year<90]+6
year[year>90]<-year[year>90]-94

# Pull out habitat 
AHAB<-as.matrix(raw.data[,grep("AH",colnames(raw.data))]) # PHAB = area older riparian growth forest (hectare)


# Pull out effort covariate for count data
OV<-data.matrix(cbind(matrix(NA,nrow=nrow(BO),ncol=(ncol(BO)-ncol(raw.data[,grep("OV",colnames(raw.data))]))),
          raw.data[,grep("OV",colnames(raw.data))]))


#Standardize habitat covariates
AHAB.std<- (AHAB - mean(AHAB)) / sd(AHAB)

# Data reshape (matrix to vector)
sample.matrix<-matrix(1:length(c(as.matrix(BO))),ncol=ncol(BO),nrow=nrow(BO))

# Make a dummy variable for type of data (occ vs. count)
type<-matrix(0,nrow=nrow(BO),ncol=ncol(BO))
type[,c(grep("BO15",colnames(BO)),grep("BO16",colnames(BO)))]<-1

# loop to reshape data and create dummy variables
obs<-years<-sites<-effort.std<-data.type<-OV.jags<-NA
for (i in 1:ncol(BO)){
  for (j in 1:nrow(BO)){
    obs[sample.matrix[j,i]]<-as.numeric(BO[j,i])
    years[sample.matrix[j,i]]<-year[i]
    OV.jags[sample.matrix[j,i]]<-OV[j,i]
    sites[sample.matrix[j,i]]<-j
    if(year[i]<21){             # count data starts in year 21
    type[sample.matrix[j,i]]<-0
    }else{
      type[sample.matrix[j,i]]<-1
      }
  }}

########################################################################
# Prepare data for JAGS
#######################################################################
# Observations are fed into the model as a vector object from 1:nSamples
nSites<-max(sites)
nYears<-max(years)
jags.data<-list(y=obs,
                nYears=nYears,
                nSites=nSites,
                count.start=grep(1,type)[1],
                count.end=length(obs),
                occ.end=grep(1,type)[1]-1 ,
                site=sites,
                year=years,
                OV=OV.jags,
                AHAB.std=AHAB.std)

# Inits only need for N in first year, S and G
inits<-function() list(N=matrix(c(rep(20,nSites),rep(NA,nSites*(nYears-1))),nrow=nSites,ncol=nYears),
                       S=matrix(10,nrow=nSites,ncol=nYears-1),
                       G=matrix(10,nrow=nSites,ncol=nYears-1))

# omitted omega, and gamma to save memory can be derived after, also omitted N
params<-c("a0","b0","b1","g0","g1","g2","p.occ","p.count","lambda","N", "gamma", "omega")

nc=3
nb=25000
ni=75000
nt=50
na=1000

library(jagsUI)
jags.out<-jags(jags.data, inits, params, model.file=model.file,store.data = TRUE, parallel=TRUE,
               n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, n.adapt=na)
traceplot(jags.out)

#Examine average annual N values
apply(jags.out$mean$N,2, mean)


