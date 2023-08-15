## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

## Libraries
library(jSDM)

## Simulating data
source("sim-data.R")
out_dir <- "outputs/R_sim_data/"
data <- sim_data_LVJSDM(n_species=30,
                        n_sites=100,
                        n_p=3,  # Number of explanatory variables (including intercept)
                        n_q=3,  # Number of latent variables
                        V_alpha_target=0.1,  # Variance of random site effects
                        neg_val_diag=-0.1, # Small negative values on the diagonal
                        seed=1234,
                        out_dir=out_dir)
X <- data$X
Y <- data$Y
W <- data$W
beta_target <- data$beta_target
lambda_target <- data$lambda_target
n_sites <- nrow(X)
n_q <- ncol(X)
n_species <- ncol(Y)

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
     main="Species effect beta",
     xlab ="obs", ylab ="fitted")
abline(a=0,b=1,col='red')
dev.off()

## Factor loadings lambda
pdf(file=file.path(out_dir, "pred_obs_lambda_mod_1.pdf"), width=10, height=10)
par(mfrow=c(2, 2))
for (l in 1:n_q) {
  plot(t(lambda_target)[, l], mean_lambda[, l],
       main=paste0("Factor loadings lambda", l),
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
Y <- read.csv(file=file.path(out_dir, "Y.csv"), header=TRUE, row.names=1)
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
     main="Species effect beta",
     xlab ="obs", ylab ="fitted")
abline(a=0,b=1,col='red')
dev.off()

## Factor loadings lambda
pdf(file=file.path(out_dir, "pred_obs_lambda_mod_2.pdf"), width=10, height=10)
par(mfrow=c(2, 2))
for (l in 1:n_q) {
  plot(t(lambda_target)[, l], mean_lambda[, l],
       main=paste0("Factor loadings lambda", l),
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
