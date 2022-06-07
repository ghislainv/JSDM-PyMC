import theano.tensor as tt

Lambda0 = tt.eye(n_species, n_q)
HALFNORMAL_SCALE = 1. / np.sqrt(1. - 2. / np.pi)

with pm.Model() as model:
    # Hyperpriors
    sigma_alpha = pm.HalfNormal("sigma_alpha", sigma=5.0)
    # Site random effect
    alpha = pm.Normal("alpha", mu=0, sigma=sigma_alpha, shape=n_sites)
    # Latent variables
    W = pm.Normal("W", mu=0, sigma=1, shape=(n_sites, n_q))
    # Factor loadings with constraints
    # Diagonal
    Lambda1 = tt.set_subtensor(Lambda0[np.arange(n_q), np.arange(n_q)],
        pm.HalfNormal("Lambda_diag", sigma=HALFNORMAL_SCALE, shape=n_q))
    # Inferior
    Lambda2 = tt.set_subtensor(Lambda1[1, 0],
        pm.Normal("Lambda_inf", mu=0, sigma=1))
    # Block
    Lambda = pm.Deterministic("Lambda", tt.set_subtensor(Lambda2[n_q:],
        pm.Normal("Lambda_block", mu=0, sigma=1, shape=(n_species-n_q, n_q))))
    # Species effects
    beta = pm.Normal("beta", mu=0, sigma=1, shape=(n_species, n_p))
    # Likelihood
    Xbeta = pm.math.dot(X, beta.transpose())
    Wlambda = pm.math.dot(W, Lambda.transpose()) 
    logit_theta = alpha[:, np.newaxis] + Xbeta + Wlambda
    obs = pm.Bernoulli("obs", logit_p=logit_theta, observed=Y)

# ===============================================================

N = n_sites * n_species
M = n_species
F = n_q
K = n_p
w_top = tt.eye(M, F)

HALFNORMAL_SCALE = 1. / np.sqrt(1. - 2. / np.pi)

with pm.Model() as rot_model:
    beta0 = pm.Normal("beta0", 0., 5., shape=M)
    beta = pm.Normal("beta", 0., 5., shape=(M, K))
    mu = beta0 + tt.dot(X, beta.T)
    
    sigma = pm.HalfNormal("sigma", 2.5)
    
    z = pm.Normal("z", 0., 1., shape=(N, F))

with rot_model:
    w_top_pos = tt.set_subtensor(
        w_top[F],
        pm.HalfNormal("w_pos", HALFNORMAL_SCALE, shape=F)
    )

with rot_model:
    w = pm.Deterministic("w", tt.set_subtensor(
        w_top_pos[F + 1:],
        pm.Normal("w_block", 0., 1., shape=(M - (F + 1), F))
    ))

w.eval()
w.tag.test_value
