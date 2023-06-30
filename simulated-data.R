## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

## Libraries
library(jSDM)

## Output directory
out_dir <- "outputs/R_sim_data/"
dir.create(out_dir)

## ==================
## Data simulation
## ==================

## Set seed
seed <- 1234
set.seed(seed)

## Variables
n_species <- 30  ## Number of species
n_sites <- 100  ## Number of sites
n_p <- 3  ## Number of explanatory variables (including intercept)
n_q <- 3  ## Number of latent variables

## Explanatory variables X
x <- matrix(rnorm(n_sites * (n_p - 1), 0, 1), nrow=n_sites, ncol=(n_p - 1))
X <- cbind(rep(1, n_sites), x)
colnames(X) <- c("Intercept", "x1", "x2")

## Latent factors W
W <- matrix(rnorm(n_sites * n_q, 0, 1), nrow=n_sites, ncol=n_q)

## Fixed species effect beta 
beta_target <- matrix(runif(n_species * n_p, -1, 1), nrow=n_p, ncol=n_species)

## Factor loadings lambda  
mat0 <- matrix(runif(n_species * n_q, -1, 1), nrow=n_q, ncol=n_species)
axis_imp <- seq(n_q, 1, -1)  ## Decreasing importance of the latent axis
mat <- mat0 * axis_imp
lambda_target <- matrix(0, n_q, n_species)
lambda_target[upper.tri(mat, diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
# Small negative values on the diagonal
diag(lambda_target) <- rep(-0.1, n_q)
# Large positive values for some species
for (i in 1:n_q) {
  lambda_target[i, n_q + i] <- n_q - i + 1
  if (i + 1 <= n_q) {
    lambda_target[(i + 1):n_q, n_q + i] <- 0
  }
}

## Variance of random site effects
V_alpha_target <- 0.1
## Random site effect alpha
alpha_target <- rnorm(n_sites, mean=0, sd=sqrt(V_alpha_target))

## Simulation of response data with probit link
## Latent variable Z
e <- matrix(rnorm(n_species * n_sites, 0, 1), n_sites, n_species)
Z <- alpha_target + X %*% beta_target + W %*% lambda_target + e
# Presence-absence matrix Y
Y <- matrix(0, n_sites, n_species)
Y[Z > 0] <- 1

## Save data
Y <- data.frame(Y)
colnames(Y) <- paste("sp", 1:n_species, sep="_")
rownames(Y) <- paste("site", 1:n_sites, sep="_")
write.csv(Y, file=file.path(out_dir, "Y.csv"), row.names=FALSE)

## =======================================
## Model 1
## =======================================

## mod_1
mod_1 <- jSDM_binomial_probit(
  # Iteration
  burnin=1000,
  mcmc=1000,
  thin=1,
  # Response variable
  presence_data=Y,
  # Explanatory variables
  site_formula=~x1+x2,
  site_data=X,
  n_latent=3,
  site_effect="random",
  # Starting values
  alpha_start=0,
  beta_start=0,
  lambda_start=0,
  W_start=0,
  V_alpha=1,
  # Priors
  shape_Valpha=0.5,
  rate_Valpha=0.0005,
  mu_beta=0, V_beta=1,
  mu_lambda=0, V_lambda=1,
  seed=1234, verbose=1)

## beta_j
mean_beta <- matrix(0, n_species, ncol(X))
for (j in 1:n_species) {
  mean_beta[j,] <- apply(mod_1$mcmc.sp[[j]]
                         [,1:ncol(X)], 2, mean)  
}

## lambda_j
mean_lambda <- matrix(0, n_species, n_q)
for (j in 1:n_species) {
  mean_lambda[j,] <- apply(mod_1$mcmc.sp[[j]]
                           [,(ncol(X)+1):(ncol(X)+n_q)], 2, mean)
}

## Species effects beta
pdf(file=file.path(out_dir, "pred_obs_beta_mod_1.pdf"), width=5, height=5)
plot(t(beta_target), mean_beta,
     main="species effect beta",
     xlab ="obs", ylab ="fitted")
abline(a=0,b=1,col='red')
dev.off()

## Factor loadings lambda
pdf(file=file.path(out_dir, "pred_obs_lambda_mod_1.pdf"), width=10, height=10)
par(mfrow=c(2, 2))
for (l in 1:n_q) {
  plot(t(lambda_target)[, l], mean_lambda[, l],
       main=paste0("factor loadings lambda", l),
       xlab ="obs", ylab ="fitted")
  abline(a=0,b=1,col='red')
}
dev.off()

## Latent variables W
pdf(file=file.path(out_dir, "pred_obs_W_mod_1.pdf"), width=10, height=10)
par(mfrow=c(2, 2))
for (l in 1:n_q) {
  plot(W[,l],
       summary(mod_1$mcmc.latent[[paste0("lv_",l)]])[[1]][,"Mean"],
       main = paste0("Latent factor W_", l),
       xlab ="obs", ylab ="fitted")
  abline(a=0,b=1,col='red')
}
dev.off()

# =======================================
# Model 2 sorting species
# =======================================

# Sorting species
Y <- read.csv(file=file.path(out_dir, "Y.csv"))
Y_sort <- Y
Y_sort[, 1] <- Y[, 4]
Y_sort[, 2] <- Y[, 5]
Y_sort[, 3] <- Y[, 6]
Y_sort[, 4] <- Y[, 1]
Y_sort[, 5] <- Y[, 2]
Y_sort[, 6] <- Y[, 3]
Y <- Y_sort

## mod_2
mod_2 <- jSDM_binomial_probit(
  # Iteration
  burnin=1000,
  mcmc=1000,
  thin=1,
  # Response variable
  presence_data=Y,
  # Explanatory variables
  site_formula=~x1+x2,
  site_data=X,
  n_latent=3,
  site_effect="random",
  # Starting values
  alpha_start=0,
  beta_start=0,
  lambda_start=0,
  W_start=0,
  V_alpha=1,
  # Priors
  shape_Valpha=0.5,
  rate_Valpha=0.0005,
  mu_beta=0, V_beta=1,
  mu_lambda=0, V_lambda=1,
  seed=1234, verbose=1)

## beta_j
mean_beta <- matrix(0, n_species, ncol(X))
for (j in 1:n_species) {
  mean_beta[j,] <- apply(mod_2$mcmc.sp[[j]]
                         [,1:ncol(X)], 2, mean)  
}

## lambda_j
mean_lambda <- matrix(0, n_species, n_q)
for (j in 1:n_species) {
  mean_lambda[j,] <- apply(mod_2$mcmc.sp[[j]]
                           [,(ncol(X)+1):(ncol(X)+n_q)], 2, mean)
}

## Species effects beta
pdf(file=file.path(out_dir, "pred_obs_beta_mod_2.pdf"), width=5, height=5)
plot(t(beta_target), mean_beta,
     main="species effect beta",
     xlab ="obs", ylab ="fitted")
abline(a=0,b=1,col='red')
dev.off()

## Factor loadings lambda
pdf(file=file.path(out_dir, "pred_obs_lambda_mod_2.pdf"), width=10, height=10)
par(mfrow=c(2, 2))
for (l in 1:n_q) {
  plot(t(lambda_target)[, l], mean_lambda[, l],
       main=paste0("factor loadings lambda", l),
       xlab ="obs", ylab ="fitted")
  abline(a=0,b=1,col='red')
}
dev.off()

## Latent variables W
pdf(file=file.path(out_dir, "pred_obs_W_mod_2.pdf"), width=10, height=10)
par(mfrow=c(2, 2))
for (l in 1:n_q) {
  plot(W[,l],
       summary(mod_2$mcmc.latent[[paste0("lv_",l)]])[[1]][,"Mean"],
       main = paste0("Latent factor W_", l),
       xlab ="obs", ylab ="fitted")
  abline(a=0,b=1,col='red')
}
dev.off()

## End
