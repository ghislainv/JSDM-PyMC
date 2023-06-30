#==================
#== Data simulation
#==================
#= Number of species
nsp <- 30
#= Number of sites
nsite <- 100
#= Number of latent variables
n_latent <- 2
#= Set seed for repeatability
seed <- 123
set.seed(seed)

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- cbind(rep(1,nsite),x1,x2)
colnames(X) <- c("Intercept","x1","x2")
np <- ncol(X)
#= Latent variables W
W <- matrix(rnorm(nsite*n_latent,0,1), nrow=nsite, ncol=n_latent)
#= Fixed species effect beta 
beta.target <- t(matrix(runif(nsp*np,-1.5,1.5), byrow=TRUE, nrow=nsp))
#= Factor loading lambda  
lambda.target <- matrix(0, n_latent,nsp)
mat <- t(matrix(runif(nsp*n_latent,-1.5,1.5), byrow=TRUE, nrow=nsp))
lambda.target[upper.tri(mat, diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
diag(lambda.target) <- rep(-0.1, n_latent)
lambda.target[1,15] <- 3
lambda.target[2,15] <- 0
lambda.target[2,25] <- 3
#= Variance of random site effect 
V_alpha.target <- 0.5
#= Random site effect alpha
alpha.target <- rnorm(nsite,0,sqrt(V_alpha.target))

# Simulation of response data with probit link
probit_theta <- X %*% beta.target + W %*% lambda.target + alpha.target
theta <- pnorm(probit_theta)
e <- matrix(rnorm(nsp*nsite,0,1),nsite,nsp)
# Latent variable Z 
Z_true <- probit_theta + e
# Presence-absence matrix Y
Y <- matrix(NA, nsite,nsp)
for (i in 1:nsite){
  for (j in 1:nsp){
    if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
    else {Y[i,j] <- 0}
  }
}
hist(probit_theta)
