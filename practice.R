cat ("hello world","\n")

# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!

source("C:/Users/Milo/Dropbox/My Documents/Programming/Kruschke/DBDA2E-utilities.R") # Kruschke's utilities
require(rjags)
require(ggplot2)
require(tidyr)
require(dplyr)
fileNameRoot="HierLab" # For output file names.

npops <- 6
nyears <- 20
spawners <- array(dim=c(npops,nyears))
recruits <- array(dim=c(npops,nyears))
pre_recruits <- array(dim=c(npops,nyears))
sp <- seq(from=0,to=6000,length.out=nyears)

# simulate data
set.seed(72146)
#alphas <- c(log(5),log(10),log(5),log(10),log(5),log(10))
alphas <- rnorm(6,mean=log(7.5),sd=1)
areas <- c(1000.,1000.,1000.,3000.,3000.,3000.)
bs <- rnorm(npops,mean=1,sd=0.2)
betas <- bs*areas
RCV <- 0.5
hrate <- 0.2
for (ipop in 1:npops) {
  spawners[ipop,1] <- runif(1)*betas[ipop]
  for (iyear in 1:nyears) {
    recruits[ipop,iyear] <- spawners[ipop,iyear]*exp(alphas[ipop]*(1-spawners[ipop,iyear]/betas[ipop])+RCV*rnorm(1))
    # pre_recruits is for plotting only
    pre_recruits[ipop,iyear] <- sp[iyear]*exp(alphas[ipop]*(1-sp[iyear]/betas[ipop]))
    if(iyear<nyears) {
      spawners[ipop,iyear+1] <- recruits[ipop,iyear]* (1-hrate)
    } #if iyear
    
  } #iyear
} #ipop


# plot stock-recruit data sets
smelt <- gather(as.data.frame(spawners)) 
rmelt <- gather(as.data.frame(recruits)) 
pmelt <- gather(as.data.frame(pre_recruits)) 
dmelt <- data.frame(pop=rep(seq(1:npops),nyears),year=rep(1:nyears,each=npops),spawners=smelt[,2],
                    recruits=rmelt[,2],pre_recruits=pmelt[,2],sp=rep(sp,each=npops))


pp <- ggplot(dmelt)  + geom_point(aes(spawners,recruits))+ geom_line(aes(sp,pre_recruits),color="red")+ 
              facet_wrap(~pop)
pp




########## JAGS code #############################################
jagsData <- list(npops=npops,nyears=nyears,spawners=spawners,recruits=recruits,areas=areas) # data for JAGS


# Define the model that JAGS will use: this is then written to a text file that JAGS reads in
# Note - this looks like R code but it isn't. The syntax for JAGS is similar but it is not an R package 
# and some R expressions won't run.
modelString = "
# data likelihood 
model {
for ( i in 1:npops ) {
  for (j in 1:nyears) {
   recruits[i,j] ~ dlnorm(pre_lnR[i,j],invsig)
   pre_lnR[i,j] <- log(spawners[i,j]) + est_alpha[i]*(1-spawners[i,j]/est_beta[i])
  } # year loop j
  est_alpha[i] ~ dnorm(theta,invtau) #priors on alpha
  est_beta[i] ~ dnorm(2000,0.0001) I(0,) #priors on beta
} # population loop i

# priors on other parameters
invsig ~ dgamma( 0.001 , 0.001) # gamma often used for prior for the inverse of the variance
sig <- pow(invsig,-0.5) # JAGS uses precision = 1/variance instead of std dev as parameter of normal
invtau ~ dgamma( 0.001 , 0.001) # gamma for hyperprior variance
tau <- pow(invtau,-0.5) # JAGS uses precision = 1/variance instead of std dev as parameter of normal
theta ~ dnorm(1, 0.1) I(0,) # vague normal for hyperprior mean (precision = 0.1 implies variance = 10, or sd ~ 3)
}
" # close quote for modelString
writeLines( modelString , con="TEMPmodel.txt" )

# Generate random starting values for a MCMC chains 
# should be plausible but contrasting sets, sort of like different starting values for a nonlinear MLE estimation
# these are passed to JAGS in list form - here as a function that generates a list
initsList = function() {
  theta <- runif(n=1,min=1.0,max=10.0)
  sig <- runif(n=1,min=1.,max=3.0)
  invsig <- 1./sig^2
  tau <- runif(n=1,min=1.,max=3.0)
  invtau <- 1./tau^2
  est_alpha <- runif(6,min=0.1,max=2)
  est_beta <- runif(6,min=500,max=4000)
  cat(theta,sig,tau,"\n")
  return( list(theta=theta,invsig=invsig,invtau=invtau,est_beta=est_beta,est_alpha=est_alpha) )
}

######### Have JAGS run 3 MCMC chains###########################
# Step 1: initializing and adaptation 
# JAGS is optimizing the MCMC sampling scheme for the current problem
jagsModel = jags.model( file="TEMPmodel.txt" , data=jagsData , inits=initsList , 
                        n.chains=3 , n.adapt=5000 ) 
# Step 2: burn-in (convergence) period, these chain values not used
# An MCMC approach converges to the desired distribution, but it takes a while.
update( jagsModel , n.iter=10000 ) 
# Step 3: generate chains of parameter values that are saved as random (but autocorrelated)
# random draws from the posterior distribution
codaSamples = coda.samples( jagsModel , variable.names=c("theta","sig","tau","est_alpha","est_beta") ,
                            n.iter=10000 ) #saved samples from the chains
save( codaSamples , file=paste0(fileNameRoot,"Mcmc.Rdata") )
#####################################################################

######## diagnose the posterior ##############
diagMCMC(codaSamples,parName="theta")
diagMCMC(codaSamples,parName="sig")
diagMCMC(codaSamples,parName="tau")

s <- as.matrix(codaSamples)
m1 <- colMeans(s)
m1 #hierarchical means of posteriors
diffalpha <- m1[1:npops]-alphas
diffalpha #difference showing shrinkage
plot(alphas,diffalpha)
diffbeta <- m1[(npops+1):(2*npops)]-betas
diffbeta #difference showing shrinkage
plot(betas,diffbeta)


######### Kruschke's R function to plot the posterior #############
# shows the hstogram and 95% high density ("credible") interval
plotPost( codaSamples[,"est_alpha[1]"] , main="alpha[1]" , xlab=bquote(alpha) )
plotPost( codaSamples[,"est_alpha[2]"] , main="alpha[2]" , xlab=bquote(alpha) )
plotPost( codaSamples[,"est_alpha[3]"] , main="alpha[3]" , xlab=bquote(alpha) )
plotPost( codaSamples[,"est_alpha[4]"] , main="alpha[4]" , xlab=bquote(alpha) )

plotPost( codaSamples[,"est_beta[1]"] , main="beta[1]" , xlab=bquote(beta) )
plotPost( codaSamples[,"est_beta[2]"] , main="beta[2]" , xlab=bquote(beta) )
plotPost( codaSamples[,"est_beta[3]"] , main="beta[3]" , xlab=bquote(beta) )
plotPost( codaSamples[,"est_beta[4]"] , main="beta[4]" , xlab=bquote(beta) )



#@@@@@@@@@@@@@@@@@@@@@ PART 2 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

########## JAGS code #############################################

# Define the model that JAGS will use: this is then written to a text file that JAGS reads in
# Note - this looks like R code but it isn't. The syntax for JAGS is similar but it is not an R package 
# and some R expressions won't run.
modelString = "
# data likelihood 
model {
for ( i in 1:npops ) {
  for (j in 1:nyears) {
   recruits[i,j] ~ dlnorm(pre_lnR[i,j],invsig)
   pre_lnR[i,j] <- log(spawners[i,j]) + est_alpha[i]*(1-spawners[i,j]/est_beta[i])
  } # year loop j
  est_alpha[i] ~ dnorm(theta,invtau) #priors on alpha
  est_b[i] ~ dnorm(thetab,invtaub)  #priors on beta
  est_beta[i] <- est_b[i]*areas[i]
} # population loop i

# priors on other parameters
invsig ~ dgamma( 0.001 , 0.001) # gamma often used for prior for the inverse of the variance
sig <- pow(invsig,-0.5) # JAGS uses precision = 1/variance instead of std dev as parameter of normal
invtau ~ dgamma( 0.001 , 0.001) # gamma for hyperprior variance
tau <- pow(invtau,-0.5) # JAGS uses precision = 1/variance instead of std dev as parameter of normal
theta ~ dnorm(1, 0.1) I(0,) # vague normal for hyperprior mean (precision = 0.1 implies variance = 10, or sd ~ 3)
invtaub ~ dgamma( 0.001 , 0.001) # gamma for hyperprior variance
taub <- pow(invtaub,-0.5) # JAGS uses precision = 1/variance instead of std dev as parameter of normal
thetab ~ dnorm(1, 1) I(0,) # vague normal for hyperprior mean (precision = 1 implies variance = 1, or sd = 1)
}
" #  quote for modelString
writeLines( modelString , con="TEMPmodel.txt" )

# Generate random starting values for a MCMC chains 
# should be plausible but contrasting sets, sort of like different starting values for a nonlinear MLE estimation
# these are passed to JAGS in list form - here as a function that generates a list
initsList = function() {
  theta <- runif(n=1,min=1.0,max=10.0)
  thetab <- runif(n=1,min=0.5,max=1.5)
  sig <- runif(n=1,min=1.,max=3.0)
  invsig <- 1./sig^2
  tau <- runif(n=1,min=1.,max=3.0)
  invtau <- 1./tau^2
  taub <- runif(n=1,min=1.,max=3.0)
  invtaub <- 1./taub^2
  est_alpha <- runif(6,min=0.1,max=2)
  est_b <- runif(6,min=0.5,max=1.5)
  cat(theta,sig,tau,"\n")
  return( list(theta=theta,invsig=invsig,invtau=invtau,thetab=thetab,invtaub=invtaub,est_b=est_b,est_alpha=est_alpha) )
}

######### Have JAGS run 3 MCMC chains###########################
# Step 1: initializing and adaptation 
# JAGS is optimizing the MCMC sampling scheme for the current problem
jagsModel = jags.model( file="TEMPmodel.txt" , data=jagsData , inits=initsList , 
                        n.chains=3 , n.adapt=5000 ) 
# Step 2: burn-in (convergence) period, these chain values not used
# An MCMC approach converges to the desired distribution, but it takes a while.
update( jagsModel , n.iter=10000 ) 
# Step 3: generate chains of parameter values that are saved as random (but autocorrelated)
# random draws from the posterior distribution
codaSamples = coda.samples( jagsModel , variable.names=c("theta","sig","tau","thetab","taub","est_alpha","est_b") ,
                            n.iter=10000 ) #saved samples from the chains
save( codaSamples , file=paste0(fileNameRoot,"Mcmc.Rdata") )
#####################################################################

s <- as.matrix(codaSamples)
m1 <- colMeans(s)
m1 #hierarchical means of posteriors
plotPost( codaSamples[,"est_b[1]"] , main="b[1]" , xlab=bquote(beta) )
plotPost( codaSamples[,"est_b[2]"] , main="b[2]" , xlab=bquote(beta) )
plotPost( codaSamples[,"est_b[3]"] , main="b[3]" , xlab=bquote(beta) )
plotPost( codaSamples[,"est_b[4]"] , main="b[4]" , xlab=bquote(beta) )
