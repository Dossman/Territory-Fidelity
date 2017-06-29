###################################################################################################
##                  Territory Fidelity - Dynamic Occupancy Approach
###################################################################################################
library(R2jags)
library(ggplot2)

## read in real data & mess with it

data <- read.csv("./data/bw_resights.csv", sep="\t")

data$Habitat <- substr(data$Plot, 1,1)

data$Date <- as.POSIXct(data$Date, "%d-%b-%y", tz="EST")

data$doy <- as.numeric(as.character(format(data$Date, "%j")))

data$Seen <- as.character(data$Seen)

data$Seen <- ifelse(data$Seen == "0*", 0, data$Seen) #--> these represent detections outside detection window

data$Seen <- as.numeric(data$Seen)

data <- data[,-c(10,12)]

data <- data[order(data$Band, data$doy),]

## Creating Data Array for Analysis

y <- array(NA, dim = c(61, 5, 6))  # Creating array with 61 Territories, Search 5 Consec Days, 
                                   # over 6 bi weekly periods

for (i in 1:5){
  for (j in 1:6){
    y[,i,j] <- data[data$Period==j & data$Day == i, "Seen"]
  }
}


## Determining how many birds were resighted at least once per period

tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == "-Inf"] <- NA
apply(tmp, 2, sum, na.rm = TRUE)

## Bundling Data for JAGS

#################################################################################################
###     Simplest Dynamic Occupancy Model - No Covariates
#################################################################################################


# Specify model in BUGS language
sink("Dynocc.jags")
cat("
    model {
    
    # Specify priors
    psi1 ~ dunif(0, 1)

    for (k in 1:(nyear-1)){
      phi[k] ~ dunif(0, 1)
      gamma[k] ~ dunif(0, 1)
      p[k] ~ dunif(0, 1) 
    }

    p[nyear] ~ dunif(0, 1)
    
    # Ecological submodel: Define state conditional on parameters

    for (i in 1:nsite){
      z[i,1] ~ dbern(psi1)
      for (k in 2:nyear){
        muZ[i,k]<- z[i,k-1]*phi[k-1] + (1-z[i,k-1])*gamma[k-1]
        z[i,k] ~ dbern(muZ[i,k])
      } #k
    } #i
    
    # Observation model
    for (i in 1:nsite){
      for (j in 1:nrep){
        for (k in 1:nyear){
          muy[i,j,k] <- z[i,k]*p[k]
          y[i,j,k] ~ dbern(muy[i,j,k])
        } #k
      } #j
    } #i
    
    # Derived parameters: Sample and population occupancy and turnover
    psi[1] <- psi1
    for (k in 2:nyear){
      psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
      turnover[k-1] <- (1 - psi[k-1]) * gamma[k-1]/psi[k]
    }
}",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = y, nsite = dim(y)[1], nrep = dim(y)[2], nyear = dim(y)[3])

# Initial values
zest <- apply(y, c(1, 3), max)
zest[is.na(zest)] <- 1
inits <- function(){ list(z = zest)}

# Parameters monitored
params <- c("psi", "phi", "gamma", "p", "turnover")  

# MCMC settings
ni <- 50000
nt <- 3
nb <- 2000
nc <- 3

# Call JAGS from R (BRT 1 min)
out <- jags(win.data, inits, params, "Dynocc.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(out, digits=2)

Period <- cbind(rep(1, out$BUGSoutput$n.sims), rep(2, out$BUGSoutput$n.sims), rep(3, out$BUGSoutput$n.sims), rep(4, out$BUGSoutput$n.sims), rep(5, out$BUGSoutput$n.sims), rep(6, out$BUGSoutput$n.sims))
boxplot(out$BUGSoutput$sims.list$gamma ~ Period, col = "gray", ylab = "Detection probability", xlab = "Period", las = 1, frame.plot = FALSE)


#################################################################################################
###     Simplest Dynamic Occupancy Model - No Covariates
#################################################################################################
tmp <- group_by(data, Band) %>% summarise(Habitat = Habitat[1]) %>% data.frame()

HABITAT <- ifelse(tmp$Habitat=="M", 1,0)

# Specify model in BUGS language
sink("Dynocc2.txt")
cat("
    model {
    
    # Specify priors
    psi1 ~ dunif(0, 1)
    

    for (k in 1:(nyear-1)){
    phi[k] ~ dunif(0, 1)
    gamma[k] ~ dunif(0, 1)
    }
    
    alpha.p ~ dnorm(0,0.1)
    beta.p ~ dnorm(0,0.1)

    # Ecological submodel: Define state conditional on parameters
    
    for (i in 1:nsite){
      z[i,1] ~ dbern(psi1)
      for (k in 2:nyear){
        muZ[i,k]<- z[i,k-1]*phi[k-1] + (1-z[i,k-1])*gamma[k-1]
        z[i,k] ~ dbern(muZ[i,k])
      } #k
    } #i
    
    # Observation model
    for (i in 1:nsite){
       lp[i] <- alpha.p + beta.p * HABITAT[i]
       lp.lim[i] <- min(999, max(-999, lp[i]))
       p[i] <- 1 / (1+ exp(-lp.lim[i]))
      for (j in 1:nrep){
        for (k in 1:nyear){
          muy[i,j,k] <- z[i,k]*p[i]
          y[i,j,k] ~ dbern(muy[i,j,k])
        } #k
      } #j
    } #i
    
    # Derived parameters: Sample and population occupancy and turnover
    psi[1] <- psi1
    for (k in 2:nyear){
    psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
    turnover[k-1] <- (1 - psi[k-1]) * gamma[k-1]/psi[k]
    }
    }",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = y, nsite = dim(y)[1], nrep = dim(y)[2], nyear = dim(y)[3], HABITAT = HABITAT)

# Initial values
zest <- apply(y, c(1, 3), max)
zest[is.na(zest)] <- 1
inits <- function(){ list(z = zest)}

# Parameters monitored
params <- c("psi", "phi", "gamma", "alpha.p","beta.p", "turnover")  

# MCMC settings
ni <- 50000
nt <- 3
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 1 min)
out <- jags(win.data, inits, params, "Dynocc2.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(out, digits=2)

par(mfrow =c(1,1))
boxplot(out$BUGSoutput$sims.list$phi, col = "gray", ylab = "Territory Survival Probability", xlab = "Period", las = 1, frame.plot = FALSE)

par(mfrow =c(2,1))
hist(plogis(out$BUGSoutput$sims.list$alpha.p),nclass=40,col="gray",main=
       "Scrub",xlab="Detection Probability",xlim=c(0.1,0.5))
hist(plogis(out$BUGSoutput$sims.list$alpha.p+out$BUGSoutput$sims.list$beta.p),
     nclass =40,col="gray",main="Mangrove",xlab="Detection
     Probability",xlim=c(0.1,0.5))


