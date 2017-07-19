rm(list=ls())

########################################################################
# Reshape data
#######################################################################

# Read in data 
model.file<-file.choose() # select BDOW_combo_model.R or the excel file

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
PHAB<-as.matrix(raw.data[,grep("PH",colnames(raw.data))]) # PHAB = percent older riparian growth forest)

# Pull out effort covariate for count data
OV<-data.matrix(cbind(matrix(NA,nrow=nrow(BO),ncol=(ncol(BO)-ncol(raw.data[,grep("OV",colnames(raw.data))]))),
          raw.data[,grep("OV",colnames(raw.data))]))

# standardize effort
OV.std<-OV
OV.std[OV>0&!is.na(OV)]<-(OV[OV>0&!is.na(OV)] - mean(OV[OV>0&!is.na(OV)],na.rm=TRUE))/sd(OV[OV>0&!is.na(OV)], na.rm=TRUE)

#Standardize habitat covariates
PHAB.std<- (PHAB - mean(PHAB)) / sd(PHAB)

# Data reshape (matrix to vector)
sample.matrix<-matrix(1:length(c(as.matrix(BO))),ncol=ncol(BO),nrow=nrow(BO))
# Make a dummy variable for type of data (occ vs. count)
type<-matrix(0,nrow=nrow(BO),ncol=ncol(BO))
type[,c(grep("BO15",colnames(BO)),grep("BO16",colnames(BO)))]<-1

# loop to reshape data and create dummy variables
obs<-years<-sites<-effort.std<-data.type<-NA
for (i in 1:ncol(BO)){
  for (j in 1:nrow(BO)){
    obs[sample.matrix[j,i]]<-as.numeric(BO[j,i])
    years[sample.matrix[j,i]]<-year[i]
    effort.std[sample.matrix[j,i]]<-as.numeric(OV.std[j,i])
    sites[sample.matrix[j,i]]<-j
    if(year[i]<21){             # count data starts in year 21
    type[sample.matrix[j,i]]<-0
    }else{
      type[sample.matrix[j,i]]<-1
      PHAB[sample.matrix[j,i]]<-NA
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
                OV=effort.std,
                PHAB.std=PHAB.std)

# Inits only need for N in first year, S and G
inits<-function() list(N=matrix(c(rep(20,nSites),rep(NA,nSites*(nYears-1))),nrow=nSites,ncol=nYears),
                       S=matrix(10,nrow=nSites,ncol=nYears-1),
                       G=matrix(10,nrow=nSites,ncol=nYears-1))

# omitted omega, and gamma to save memory can be derived after, also omitted N
params<-c("a0","a1","b0","b1","g0","g1","g2","p.occ","p.count","lambda")

nc=3
nb=25000
ni=75000
nt=50
na=1000

###################################################################
#Run the model in JAGS
##################################################################

library(jagsUI)
jags.out<-jags(jags.data, inits, params, model.file=model.file,store.data = TRUE,# parallel=TRUE,
               n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, n.adapt=na)
