##R code to process results and create panels for figure 4

##############################
#gamma
#############################

#Figure 4b - gains related to abundance

# calculate gamma across range of covariate
nIter<-length(jags.out$sims.list$a0)
Ns<-seq(1,6,by=0.1)
gam<-matrix(NA, nrow=nIter, ncol=length(Ns))
for (i in 1:nIter){
  gam[i,]<-exp(jags.out$sims.list$g0[i] + jags.out$sims.list$g1[i]*Ns + jags.out$sims.list$g2[i]*Ns*Ns)
}
gam.q<-matrix(NA, ncol=5, nrow=ncol(gam))
for (i in 1:ncol(gam)){
  gam.q[i,]<-quantile(gam[,i],probs=c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE)
}
#plot
par(las=1)
xlims=c(1,6)
ylims=c(0,10)
plot(NA,NA,ylim=ylims,xlim=xlims,xlab="",ylab="", bty="n", main="",
     xaxt="n",yaxt="n",xaxs="i",yaxs="i", cex.lab=2, cex.main=2)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey100")#, border=NA)
# 95% CI
polygon(c(Ns,rev(Ns)),c(gam.q[,1],rev(gam.q[,5])),border=NA,col='grey80') 
#50% CI
polygon(c(Ns,rev(Ns)),c(gam.q[,2],rev(gam.q[,4])),border=NA,col='grey50') 
# Median
lines(Ns,gam.q[,3],lwd=3)
#label
axis(1,1:6,cex.axis=1.5)
axis(2,seq(0,100,by=2), cex.axis=1.5)
mtext(expression(bar(N)[t-1]), side=1, cex=1.5, padj =2.5)
par(las=3)
mtext(side=2,expression("  Gains ("*gamma*")"),cex=1.5, line=2.4)
par(las=1)



#########################################
#omega
########################################

#Figure 4c - survival related to older growth forest (ha)

# Calcualte survival across range of covariate (Area older riparian growth forrest)
range(AHAB, na.rm=TRUE)
AHAB.mean<-mean(AHAB[AHAB>0&!is.na(AHAB)])
AHAB.sd<-sd(AHAB[AHAB>0&!is.na(AHAB)])
o.steps<-seq(-2,3,0.01)
omega<-matrix(NA, ncol=length(o.steps), nrow=nIter)
for (i in 1:nIter){
  omega[i,]<-plogis(jags.out$sims.list$b0[i] + jags.out$sims.list$b1[i]*o.steps)
}
omega.q<-matrix(NA, nrow=length(o.steps), ncol=5)
for (i in 1:length(o.steps)){
  omega.q[i,]<-quantile(omega[,i],probs=c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE)
}
real.steps<-o.steps*AHAB.sd + AHAB.mean
par(las=1)
xlims=c(0,700)
ylims=c(0.7,1)
plot(NA,NA,ylim=ylims,xlim=xlims,xlab="",ylab="", bty="n", main="",
     xaxt="n",yaxt="n",xaxs="i",yaxs="i", cex.lab=2, cex.main=2)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey100")#, border=NA)
# 95% CI
polygon(c(real.steps,rev(real.steps)),c(omega.q[,1],rev(omega.q[,5])),border=NA,col='grey80') 
# 50% CI
polygon(c(real.steps,rev(real.steps)),c(omega.q[,2],rev(omega.q[,4])),border=NA,col='grey50') 
# Median 
lines(real.steps,omega.q[,3],lwd=3)
# Label
axis(1, at=seq(0,700,100),cex.axis=1.5)
axis(2,seq(0,1,by=0.1), cex.axis=1.5)
mtext("Riparian Old Growth (ha)", side=1, cex=1.5, padj =2.5)
par(las=3)
mtext(side=2,expression("  Survival ("*omega*")"),cex=1.5, line=2.6)
par(las=1)


##################################
#colonization extinction
##################################

#Figure 4d


# Calculate extinction and colonization for ever year/site 
col<-matrix(NA, ncol=nYears-1, nrow=nIter)
ext<-array(NA,dim=c(nIter, nSites, nYears-1))

# colonization (does not vary by site)
for (t in 2:nYears){
  col[,t-1]<-1-exp(-(jags.out$sims.list$gamma[,t-1]))
}

# Extinction (varies by site based on N in previous year)
for (t in 2:nYears){
  for (j in 1:nSites){
    for (k in 1:nIter){
      if (jags.out$sims.list$N[k,j,t-1]>0){
        ext[k,j,t-1]<-(1-jags.out$sims.list$omega[k,j,t-1])^jags.out$sims.list$N[k,j,t-1]*exp(-jags.out$sims.list$gamma[k,t-1])
      }}}}

# Calcualte quantiles
col.q<-matrix(NA,nrow=nYears-1,ncol=6)
ext.q<-matrix(NA,nrow=nYears-1,ncol=6)
for (t in 1:nYears-1){
  col.q[t,]<-c(c(1996:2016)[t],quantile(col[,t],probs=c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE))
  ext.q[t,]<-c(c(1996:2016)[t],quantile(ext[,,t],probs=c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE)) 
}

#plot
par(las=1)
xlims=c(1995.5,2016.5)
ylims=c(-0.01,1)
plot(NA,NA,ylim=ylims,xlim=xlims,xlab="",ylab="", bty="n", main="",
     xaxt="n",yaxt="n",xaxs="i",yaxs="i", cex.lab=2, cex.main=2)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey100")#, border=NA)
# 95% CI
segments(col.q[,1]-0.1,col.q[,2],col.q[,1]-0.1,col.q[,6],col="grey40",lwd=2)
segments(ext.q[,1]+0.1,ext.q[,2],ext.q[,1]+0.1,ext.q[,6],col="black",lwd=2)
# median
points(col.q[,1]-0.1,col.q[,4],pch=19,lwd=1,cex=2,col="grey40")
points(ext.q[,1]+0.1,ext.q[,4],pch=23,lwd=1,cex=2,col="black",bg="black")
#label
axis(1,c(1996,seq(1995,2015,by=5)),cex.axis=1.5)
axis(2,seq(0,100,by=0.2), cex.axis=1.5)
mtext("Year", side=1, cex=2, padj =2.5)
par(las=3)
mtext("Col/Ext Prob", side=2, cex=2, line=2.75)
par(las=1)
mtext("Colonization", col="grey40", line=-2, side=3)
mtext("Extinction", col="black", line=-1, side=3)


#######################################
# Detection
#######################################

#Figure 4e


# Vector to loop results across
nIter<-length(jags.out$sims.list$a0)

# Calculate detection across range of covariate for effort
eff.range<-seq(0.1,1,0.01)
eff<-matrix(NA, nrow=nIter, ncol=length(eff.range))
for (i in 1:length(jags.out$sims.list$a0)){
  eff[i,]<-1-exp(-exp(  jags.out$sims.list$a0[i] + log(eff.range) ) )   
}

# Calculate quantiles
eff.q<-matrix(NA, nrow=ncol(eff), ncol=6)
for (i in 1:nrow(eff.q)){
  eff.q[i,2:6]<-quantile(eff[,i],probs=c(0.025,0.25,0.5,0.75,0.975))
  eff.q[i,4]<-mean(eff[,i],probs=c(0.025,0.25,0.5,0.75,0.975))
}
eff.q[,1]<-eff.range

#plot
par(las=1)
xlims=c(0.05,1.5)
ylims=c(0,0.6)
plot(NA,NA,ylim=ylims,xlim=xlims,xlab="",ylab="", bty="n", main="",
     xaxt="n",yaxt="n",xaxs="i",yaxs="i", cex.lab=2, cex.main=2)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey100")#, border=NA)

# 95% CI
polygon(c(eff.q[,1],rev(eff.q[,1])), c(eff.q[,2], rev(eff.q[,6])), border=NA, col='grey80')
# 50% CI
polygon(c(eff.q[,1],rev(eff.q[,1])), c(eff.q[,3], rev(eff.q[,5])), border=NA, col='grey50')
# Median
lines(eff.q[,1], eff.q[,4], lwd=2, col="black")
# boxplot of occupancy detection
boxplot(jags.out$sims.list$p.occ,add=T,at=1.25,outline=F, yaxt="n",col="blue")

# Label
axis(1,seq(0,1,0.2), labels=seq(0,100,20),cex.axis=1.5)
axis(1, at=1.25, labels=c(expression("    P"[occ])),cex.axis=1.5)
axis(2,seq(0,1,by=0.1), cex.axis=1.5)
mtext("% Site Surveyed", side=1, cex=1.5, at=0.5, line=2.25)
mtext(expression("P"[count]), side=1, cex=1.5, at=0.5, line=3.65)
abline(v=1.025,lty=2)
par(las=3)
mtext("Detection", side=2, cex=1.5, line=2.75)
par(las=1)




