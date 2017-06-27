###################################################################################################
##                  Simulating Dynamic Occupancy Process
###################################################################################################
library(R2jags)
library(ggplot2)

## Function to simulate dynamic occupancy data

data.fn <- function(R = 250, J = 3, K = 10, psi1 = 0.4, range.p = c(0.2, 0.4), range.phi = c(0.6, 0.8), range.gamma = c(0, 0.1)) {
  # Function to simulate detection/nondetection data for dynamic site-occ model
  # Annual variation in probabilities of patch survival, colonization and 
  # detection is specified by the bounds of a uniform distribution.
  
  # Function arguments:
  # R - Number of sites
  # J - Number of replicate surveys
  # K - Number of years
  # psi1 - occupancy probability in first year
  # range.p - bounds of uniform distribution from which annual p drawn 
  # range.psi and range.gamma - same for survival and colonization probability
  
  # Set up some required arrays
  site <- 1:R					# Sites
  year <- 1:K					# Years
  psi <- rep(NA, K)				# Occupancy probability
  muZ <- z <- array(dim = c(R, K))	# Expected and realized occurrence
  y <- array(NA, dim = c(R, J, K))	# Detection histories
  
  # Determine initial occupancy and demographic parameters
  psi[1] <- psi1				# Initial occupancy probability
  p <- runif(n = K, min = range.p[1], max = range.p[2])
  phi <- runif(n = K-1, min = range.phi[1], max = range.phi[2])
  gamma <- runif(n = K-1, min = range.gamma[1], max = range.gamma[2])
  
  # Generate latent states of occurrence
  # First year
  z[,1] <- rbinom(R, 1, psi[1])		# Initial occupancy state
  # Later years
  for(i in 1:R){				# Loop over sites
    for(k in 2:K){				# Loop over years
      muZ[k] <- z[i, k-1]*phi[k-1] + (1-z[i, k-1])*gamma[k-1] # Prob for occ.
      z[i,k] <- rbinom(1, 1, muZ[k])
    }
  }
  
  # Plot realised occupancy
  plot(year, apply(z, 2, mean), type = "l", xlab = "Year", ylab = "Occupancy or Detection prob.", col = "red", xlim = c(0,K+1), ylim = c(0,1), lwd = 2, lty = 1, frame.plot = FALSE, las = 1)
  lines(year, p , type = "l", col = "red", lwd = 2, lty = 2)
  
  # Generate detection/nondetection data
  for(i in 1:R){
    for(k in 1:K){
      prob <- z[i,k] * p[k]
      for(j in 1:J){
        y[i,j,k] <- rbinom(1, 1, prob)
      }
    }
  }
  
  # Compute annual population occupancy
  for (k in 2:K){
    psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
  }
  
  # Plot apparent occupancy
  psi.app <- apply(apply(y, c(1,3), max), 2, mean)
  lines(year, psi.app, type = "l", col = "black", lwd = 2)
  text(0.85*K, 0.06, labels = "red solid - true occupancy\n red dashed - detection\n black - observed occupancy")
  
  # Return data
  return(list(R = R, J = J, K = K, psi = psi, psi.app = psi.app, z = z, phi = phi, gamma = gamma, p = p, y = y))
}


data <-data.fn(R=60,J=5,K=6,psi1=0.7,range.p=c(0.1,0.9),range.phi =c(0.6,0.9),range.gamma=c(0.1,0.5))

str(data)



## Bayesian Model

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
    
    # Derived parameters: Sample and population occupancy, growth rate and turnover
    psi[1] <- psi1
    n.occ[1]<-sum(z[1:nsite,1])
    for (k in 2:nyear){
      psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
      n.occ[k] <- sum(z[1:nsite,k])
      growthr[k-1] <- psi[k]/psi[k-1]                         # originally we had growthr[k]. JAGS seem to dislike vectoring going from 2..K.
      turnover[k-1] <- (1 - psi[k-1]) * gamma[k-1]/psi[k]
    }
}
    ",fill = TRUE)
sink()

y <- data$y

jdata <- list(y=y, nsite = dim(y)[1], nrep = dim(y)[2], nyear = dim(y)[3])

zst <-apply(y,c(1,3),max) #Observedoccurrenceasinitsforz

inits <-function(){list(z=zst)}

# Parametersmonitored
params <- c("psi", "phi", "gamma", "p", "n.occ", "growthr", "turnover") 



# MCMC settings
ni <- 2500
nt <- 4
nb <- 500
nc <- 3

# Call JAGS from R (BRT 3 min)
out <- jags(jdata, inits, params, "Dynocc.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())


# Summarize posteriors
print(out, dig = 2)
psiall <- paste("psi[", 1:K, "]", sep="")
print(cbind(data$psi, out$BUGSoutput$summary[psiall, c(1, 2, 3, 7)]), dig = 3)
phiall <- paste("phi[", 1:(K-1), "]", sep="")
print(cbind(data$phi, out$BUGSoutput$summary[phiall, c(1, 2, 3, 7)]), dig = 3)
gammaall <- paste("gamma[", 1:(K-1), "]", sep="")
print(cbind(data$gamma, out$BUGSoutput$summary[gammaall, c(1, 2, 3, 7)]), dig = 3)
pall <- paste("p[", 1:K, "]", sep="")
print(cbind(data$p, out$BUGSoutput$summary[pall, c(1, 2, 3, 7)]), dig = 3)

plot(1:K, data$psi, type = "l", xlab = "Year", ylab = "Occupancy probability", col = "red", xlim = c(0,K+1), ylim = c(0,1), lwd = 2, lty = 1, frame.plot = FALSE, las = 1)
lines(1:K, data$psi.app, type = "l", col = "black", lwd = 2)
points(1:K, out$BUGSoutput$mean$psi, type = "l", col = "blue", lwd = 2)
segments(1:K, out$BUGSoutput$summary[psiall,3], 1:K, out$BUGSoutput$summary[psiall,7], col = "blue", lwd = 1)


### Analyzing Real Data

bdat <- read.table(file = "~/Desktop/burnet.txt", header = T)
str(bdat)

y <- array(NA, dim = c(95, 2, 7))	# 95 sites, 2 reps, 7 days

for (i in 1:7){
  sel.rows <- bdat$day == i
  y[,,i] <- as.matrix(bdat)[sel.rows, 3:4]
}
str(y)

# Convert counts to detection/nondetection data
y[y>0] <- 1

# Look at the number of sites with detections for each day
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == "-Inf"] <- NA
apply(tmp, 2, sum, na.rm = TRUE)

# Bundle data
win.data <- list(y = y, nsite = dim(y)[1], nrep = dim(y)[2], nyear = dim(y)[3])

# Initial values
inits <- function(){ list(z = apply(y, c(1, 3), max))}

# Parameters monitored
params <- c("psi", "phi", "gamma", "p", "n.occ", "growthr", "turnover")  

# MCMC settings
ni <- 5000
nt <- 4
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 1 min)
out1 <- jags(win.data, inits, params, "Dynocc.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out1, dig = 3)

# Specify model in BUGS language
sink("Dynocc2.jags")
cat("
    model {
    
    # Specify priors
    psi1 ~ dunif(0, 1)
    for (k in 1:(nyear-1)){
    phi[k] ~ dunif(0, 1)
    gamma[k] ~ dunif(0, 1)
    }
    p ~ dunif(0, 1)
    
    # Both models at once
    for (i in 1:nsite){
    z[i,1] ~ dbern(psi1)     # State model 1: Initial state
    muy[i,1] <- z[i,1]*p
    y[i,1] ~ dbin(muy[i,1], 2)
    for (k in 2:nyear){      # State model 2: State dynamics
    muZ[i,k] <- z[i,k-1]*phi[k-1] + (1-z[i,k-1])*gamma[k-1]
    z[i,k] ~ dbern(muZ[i,k])
    
    # Observation model
    muy[i,k] <- z[i,k]*p
    y[i,k] ~ dbin(muy[i,k], 2)
    } #k
    } #i
    
    # Derived parameters: Sample and population occupancy, growth rate and turnover
    psi[1] <- psi1
    n.occ[1] <- sum(z[1:nsite,1])
    for (k in 2:nyear){
    psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
    n.occ[k] <- sum(z[1:nsite,k])
    growthr[k-1] <- psi[k]/psi[k-1]
    turnover[k-1] <- (1 - psi[k-1]) * gamma[k-1]/psi[k]
    }
    } # end of model
    ",fill = TRUE)
sink()

# Aggregate detections over reps within a day and bundle data
yy <- apply(y, c(1, 3), sum, na.rm = TRUE)
win.data <- list(y = yy, nsite = dim(yy)[1], nyear = dim(yy)[2])

# Initial values
inits <- function(){list(z = apply(y, c(1, 3), max))}

# Parameters monitored
params <- c("psi", "phi", "gamma", "p", "n.occ", "growthr", "turnover")  

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R (BRT 1 min)
out2 <- jags(win.data, inits, params, "Dynocc2.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posterior
print(out2, dig = 3)

DAY <- cbind(rep(1, out2$BUGSoutput$n.sims), rep(2, out2$BUGSoutput$n.sims), rep(3, out2$BUGSoutput$n.sims), rep(4, out2$BUGSoutput$n.sims), rep(5, out2$BUGSoutput$n.sims), rep(6, out2$BUGSoutput$n.sims), rep(7, out2$BUGSoutput$n.sims))
boxplot(out2$BUGSoutput$sims.list$psi ~ DAY, col = "gray", ylab = "Occupancy probability", xlab = "Day of survey", las = 1, frame.plot = FALSE)
apply(apply(y, c(1, 3), max), 2, function(x){sum(!is.na(x))})


