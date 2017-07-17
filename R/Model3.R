library(ggplot2)
library(dplyr)
library(R2jags)

y <- readRDS("./data/territory_data.rds")
data <- readRDS("./data/full_data.rds")

tmp <- group_by(data, Band) %>% summarise(Habitat = Habitat[1], Sex = Sex[1]) %>% data.frame()

HABITAT <- ifelse(tmp$Habitat=="M", 1,0)
SEX <- ifelse(tmp$Sex == "M", 1, 0)

rain <- read.csv("~/Desktop/2017_RAIN.csv", header=T)

rain <- group_by(rain, Period) %>% summarise(mean = sum(Mean, na.rm=T))%>% data.frame()

RAIN <- rain[2:6,"mean"]

RAIN <- (RAIN - mean(RAIN))/sd(RAIN)

## Specify model in BUGS language
sink("Dynocc3.txt")
cat("
    model {
    
    # Specify priors
    psi1 ~ dunif(0, 1)

    for (k in 2:nyear){
    phi[k-1] ~ dunif(0,1)
    }

    alpha.gamma ~ dnorm(0,0.1)
    beta.gamma ~ dnorm(0,0.1)

    alpha.p ~ dnorm(0,0.1)
    beta.p ~ dnorm(0,0.1)
    
    # Ecological submodel: Define state conditional on parameters
    for (i in 1:nsite){
      for (k in 2:nyear){
        logit(gamma[i,k-1]) <- alpha.gamma + beta.gamma*RAIN[i,k-1] + beta.gamma2*HABITAT[i,k]
      }    
    }

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
    }
}",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = y, nsite = dim(y)[1], nrep = dim(y)[2], nyear = dim(y)[3], HABITAT = HABITAT, RAIN = RAIN)

# Initial values
zest <- apply(y, c(1, 3), max)
zest[is.na(zest)] <- 1
inits <- function(){ list(z = zest)}

# Parameters monitored
params <- c("phi","psi", "gamma", "alpha.p","beta.p","alpha.gamma","beta.gamma")  

# MCMC settings
ni <- 10000
nt <- 3
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 1 min)
out <- jags(win.data, inits, params, "Dynocc3.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(out,  digits=2)

par(mfrow =c(1,1))
boxplot(out$BUGSoutput$sims.list$gamma, col = "gray", ylab = "Territory Survival Probability", xlab = "Period", las = 1, frame.plot = FALSE)

plot(out$BUGSoutput$mean$gamma~RAIN)

par(mfrow =c(2,1))
hist(plogis(out$BUGSoutput$sims.list$alpha.p),nclass=40,col="gray",main=
       "Scrub",xlab="Detection Probability",xlim=c(0.1,0.5))
hist(plogis(out$BUGSoutput$sims.list$alpha.p+out$BUGSoutput$sims.list$beta.p),
     nclass =40,col="gray",main="Mangrove",xlab="Detection Probability",xlim=c(0.1,0.5))


