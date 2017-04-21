library(runjags)
library(R2jags)
library(rjags)

##################################################################################################
## Simulate Data for a single season occupancy model
##################################################################################################

## Select Sample Sizes (spatial replicates-R and temporal replicates-T)
R = 200
T = 6

## Generating the process parameters
psi = 0.5 #occupancy
p = 0.9 #detection probability

# Create data matrix
y <-matrix(NA,nrow=R,ncol=T)

# Ecological process: Sample true occurrence(z,yes/no) from a Bernoulli (occurrence probability=psi)
z = rbinom(n = R, size = 1, prob = psi)

# Observation process: Sample detection/nondetection observations from a Bernoulli(with p) if z=1

for (j in 1:T){
  y[,j] <- rbinom(n = R, size = 1, prob = z*p)
}

##################################################################################################
## Analyze Data with JAGS (runjags)
##################################################################################################

## Specify model in BUGS language
sink("model.jags")
cat("
    model{
    
    ## Specify Priors
    psi ~ dunif(0,1)
    p   ~ dunif(0,1)
    
    ##Likelihood
    ## Ecological model for true occurence
    for (i in 1:R){
    z[i] ~ dbern(psi)
    p.eff[i] <- z[i] * p
    
    ## Observation model for replicated detection/nondetection
    
    for (j in 1:T){
    y[i,j] ~ dbern(p.eff[i])
    }
    }
    
    ## Derived Quantities
    occ.fs <-sum(z[]) # Number of occupied sites among the R sites
    }
    ",fill =TRUE)
sink()

# Bundle data
data = list(y = y, R = nrow(y), T = ncol(y))

zst <-apply(y,1,max)  #Observed occurrence as starting values for z
inits <-function()list(z=zst)

params <-c("psi","p","occ.fs")

# MCMCsettings
ni <-1500
nt <-2
nb <-500
nc <-3

out <- jags(data, inits, params, "model.jags", 
            n.chains = nc, n.thin = nt, n.iter = ni, 
            n.burnin = nb, working.directory = getwd())

print(out, dig = 2)
