###################################################################################################
##                  Simulating Single Season Occupancy Data
###################################################################################################
library(R2jags)
library(ggplot2)

# Selecting number of sites visited and number of surveys done
R = 60
T = 7

# Setting occupancy probability and detection probability
p = 0.4
psi = 0.7

# Setting up Observation Matrix

y <-matrix(NA,nrow=R,ncol=T)

# Ecologicalprocess:Sample true occurrence(z,yes/no) from a Bernoulli (occurrenceprobability=psi)

z <-rbinom(n=R,size=1,prob=psi)

# Observation Process

for (j in 1:T){
  y[,j] <-rbinom(n = R ,size = 1, prob = z * p)
}

# Summarizing 'TRUE' data

sum(z) # Realized Occupancy

sum(apply(y,1,max)) # Observed occupancy

# Build Bayesian Model in jags

sink("model.txt")
cat("
    model {
    
    # Priors
      psi ~ dunif(0,1)
      p ~ dunif(0,1)
    
    # Likelihood
     for (i in 1:R) {
        # True state model for the partially observed true state
        z[i] ~ dbern(psi)		# True occurrence z at site i
    
        for (j in 1:T) {
          # Observation model for the actual observations
            y[i,j] ~ dbern(p.eff[i,j])	# Detection-nondetection at i and j
            p.eff[i,j] <- z[i] * p
        }
      }
    
    # Derived quantities
    occ.fs <- sum(z[])	# Number of occupied sites among 60 studied
    }
    ",fill=TRUE)
sink()

# data 
data <-list(y=yy,R=nrow(yy),T=ncol(yy))

# Initial values
zst <-apply(y,1,max)
zst[is.na(zst)] <- 1 #Observedoccurrenceasstartingvaluesforz
inits <-function()list(z=zst)

# Parametersmonitored
params <-c("psi","p","occ.fs")

# MCMC settings
ni <- 10000
nt <- 4
nb <- 500
nc <- 3

# Call JAGS from R (BRT < 1 min)
out <- jags(data, inits, params, "model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, digits = 2)

###################################################################################################
####                            Adding Covariate effects
###################################################################################################
data.fn <-function(R=200,T=3,xmin= -1, xmax=1,alpha.psi= -1, beta.psi =3,alpha.p=1,beta.p= -3) {

y <- array(dim=c(R,T)) # An array for the count data

X <- sort(runif(n = R, min = xmin, max = xmax))

psi <- plogis(alpha.psi + beta.psi*X) # apply inverse logit

z <- rbinom(n = R, size = 1, prob = psi)

occ.fs <- sum(z) # constrains sample to finite occupancy

p <- plogis(alpha.p + beta.p * X)

p.eff <- z*p

for (i in 1:T){
  y[,i] <- rbinom(n = R, size = 1, prob = p.eff)
}

naive.pred <- plogis(predict(glm(apply(y,1,max) ~ X + I(X^2), family = binomial)))

par(mfrow =c(2,2))
plot(X, psi,main="Expectedoccurrence",xlab="Covariate",ylab ="Occupancyprobability",las=1,type="l",col="red",lwd =3,frame.plot=FALSE)
plot(X, z,main="Realised(true)occurrence",xlab="Covariate",ylab ="Occurrence",las=1,frame.plot=FALSE,col="red")
plot(X, p,ylim=c(0,1),main="Detectionprobability",xlab ="Covariate",ylab="p",type="l",lwd=3,col="red",las =1,frame.plot=FALSE)
plot(X, naive.pred,main="Detection/nondetectionobservations\nand conventionalSDM",xlab="Covariate",ylab="Apparentoccupancy",ylim=c(min(y),max(y)),type="l",lwd=3,lty=2,col ="blue",las=1,frame.plot=FALSE)
points(rep(X,T),y)

return(list(R=R,T=T,X=X,alpha.psi=alpha.psi,beta.psi=
              beta.psi, alpha.p=alpha.p,beta.p=beta.p,psi=psi,z=z,
            occ.fs =occ.fs,p=p,y=y))
}

sink("model.txt")
cat("
    model{
      #PRIORS
        alpha.occ <- dunif(-10,10)
        beta.occ  <- dunif(-10,10)
        alpha.p   <- dunif(-10,10)
        beta.p    <- dunif(-10,10)
    }
    ")




##################################################################################################
##      Analyzing real data with Covariates
##################################################################################################

data <- read.table("~/Desktop/bluebug.txt", header=T)

y <-as.matrix(data[,4:9])#as.matrixessentialforWinBUGS
y[y>1] <-1#Reducecountsto0/1
edge <-data$forest_edge
dates <-as.matrix(data[,10:15])
hours <-as.matrix(data[,16:21])

# Standardizecovariates
mean.date <-mean(dates,na.rm=TRUE)
sd.date <-sd(dates[!is.na(dates)])
DATES <-(dates-mean.date)/sd.date#Standardisedate
DATES[is.na(DATES)]<-0#Imputezeroes(means)
mean.hour <-mean(hours,na.rm=TRUE)
sd.hour <-sd(hours[!is.na(hours)])
HOURS <-(hours-mean.hour)/sd.hour#Standardisehour
HOURS[is.na(HOURS)]<-0#Imputezeroes(means)


sink("model.txt")
cat("
    model {
      #Priors

      alpha.psi ~ dnorm(0, 0.01)
      beta.psi ~ dnorm(0, 0.01)
      alpha.p ~ dnorm(0, 0.01)
      beta1.p ~ dnorm(0, 0.01)
      beta2.p ~ dnorm(0, 0.01)
      beta3.p ~ dnorm(0, 0.01)
      beta4.p ~ dnorm(0, 0.01)

      # Likelihood

      for (i in 1:R){
        z[i] ~ dbern(psi[i])
        psi[i] <- 1 / (1 + exp(-lpsi.lim[i]))
        lpsi.lim[i] <- min(999, max(lpsi[i], -999))
        lpsi[i] <- alpha.psi + beta.psi * edge[i]
          for(j in 1:T){
            y[i,j] ~ dbern(mu.p[i,j])
            mu.p[i,j] <- z[i] * p[i,j]
            p[i,j] <- 1 / (1+ exp(-lp.lim[i,j]))
            lp.lim[i,j] <- min(999, max(-999, lp[i,j]))
            lp[i,j] <- alpha.p + beta1.p * DATES[i,j] + beta2.p * pow(DATES[i,j],2) + beta3.p * HOURS[i,j] + beta4.p * pow(HOURS[i,j],2)
      }
      }
occ.fs <-sum(z[]) # Numberofoccupiedsites
mean.p <-exp(alpha.p)/(1+exp(alpha.p)) #Averagedetection
}",fill =TRUE)
sink()

jdata <- list(y=y, R=nrow(y), T=ncol(y), edge=edge, DATES=DATES, HOURS=HOURS)

zst <-apply(y,1,max,na.rm=TRUE)#Goodstartingvaluescrucial
inits <-function(){list(z=zst,alpha.psi=runif(1, -3, 3),alpha.p=runif(1, -3, 3))}


params <-c("alpha.psi","beta.psi","mean.p","occ.fs","alpha.p","beta1.p","beta2.p","beta3.p","beta4.p")

ni <-30000
nt <-10
nb <-20000
nc <-3

out <- jags(jdata, inits, params, "model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

hist(out$BUGSoutput$sims.list$occ.fs,nclass=30,col="gray",main="",xlab="Number ofoccupiedwoodpiles(occ.fs)",xlim=c(9,27))
abline(v =10,lwd=2)#Theobservednumber

Pstar <-array(NA,dim=c(out$BUGSoutput$n.sims,10))
x <-cbind(rep(1,3000),rep(2,3000),rep(3,3000),rep(4,3000),rep
          (5, 3000),rep(6,3000),rep(7,3000),rep(8,3000),rep(9,3000),
          rep(10, 3000))
for (i in 1:out$BUGSoutput$n.sims){
  for (j in 1:10){
    Pstar[i,j]<-1 - (1 - out$BUGSoutput$sims.list$mean.p[i])^j
  } #j
} #i

boxplot(Pstar~x,col="gray",las=1,ylab="Pstar",xlab="Numberofsurveys", outline=FALSE)
abline(h =0.95,lty=2,lwd=2)
