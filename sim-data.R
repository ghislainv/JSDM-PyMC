## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

sim_data_LVJSDM <- function(n_species=30,
                            n_sites=100,
                            n_p=3,  # Number of explanatory variables (including intercept)
                            n_q=3,  # Number of latent variables
                            V_alpha_target=0.1,  # Variance of random site effects
                            neg_val_diag=-0.1, # Small negative values on the diagonal
                            seed=1234,
                            out_dir="outputs/R_sim_data/") {

  ## Output directory
  dir.create(out_dir)
  
  ## Set seed
  set.seed(seed)

  ## Explanatory variables X
  x <- matrix(rnorm(n_sites * (n_p - 1), 0, 1), nrow=n_sites, ncol=(n_p - 1))
  X <- cbind(rep(1, n_sites), x)
  colnames(X) <- c("Intercept", "x1", "x2")
  rownames(X) <- paste("site", 1:n_sites, sep="_")

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
  diag(lambda_target) <- rep(neg_val_diag, n_q)
  # Large positive values for some species
  for (i in 1:n_q) {
    lambda_target[i, n_q + i] <- n_q - i + 1
    if (i + 1 <= n_q) {
      lambda_target[(i + 1):n_q, n_q + i] <- 0
    }
  }

  ## Random site effect alpha
  alpha_target <- rnorm(n_sites, mean=0, sd=sqrt(V_alpha_target))

  ## Simulation of response data with probit link
  ## Latent variable Z
  e <- matrix(rnorm(n_species * n_sites, 0, 1), n_sites, n_species)
  Z <- alpha_target + X %*% beta_target + W %*% lambda_target + e
  # Presence-absence matrix Y
  Y <- matrix(0, n_sites, n_species)
  Y[Z > 0] <- 1
  Y <- data.frame(Y)
  colnames(Y) <- paste("sp", 1:n_species, sep="_")
  rownames(Y) <- paste("site", 1:n_sites, sep="_")

  ## Save data
  write.csv(X, file=file.path(out_dir, "X.csv"), row.names=TRUE)
  write.csv(Y, file=file.path(out_dir, "Y.csv"), row.names=TRUE)

  return(list(X=X, Y=Y, W=W,
              beta_target=beta_target,
              lambda_target=lambda_target))
  
}
  
## End
