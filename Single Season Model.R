## Simulating Single Season Occupancy Data

# Selecting number of sites visited and number of surveys done
R = 60
T = 10

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
data <-list(y=y,R=nrow(y),T=ncol(y))

# Initial values
zst <-apply(y,1,max) #Observedoccurrenceasstartingvaluesforz
inits <-function()list(z=zst)
# Parametersmonitored
params <-c("psi","p","occ.fs")

# MCMC settings
ni <- 2000
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R (BRT < 1 min)
out <- jags(data, inits, params, "model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, digits = 2)

plot(out)




