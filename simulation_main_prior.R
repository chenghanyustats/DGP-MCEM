################################################################################
# DGP Simulation Study: Prior effects                                          #
# Cheng-Han Yu                                                                 #
################################################################################
## data referred to simulation_main.R

#############################################
## Algorithm implementation 
#############################################
library(parallel)
library(doParallel)
detectCores()
cl <- makeCluster(6, type = "FORK")
registerDoParallel(cl)
getDoParWorkers()
no_data <- 100

# ===========
## GPR no der
# ===========
# EB_gp <- foreach(k = 1:no_data) %dopar% {
#     res <- Rsolnp::solnp(pars = c(1, 1, 1), fun = log_mar_lik_gp,
#                          LB = c(0.0001, 0.0001, 0.0001), 
#                          UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001), 
#                          control = list(TOL = 1e-5, trace = 0),
#                          y = YY[[k]]$y, H0 = H0_diff[[k]])
#     res$par
# }

# ===========
## DGP-oracle
# ===========
# EB_dgp_oracle <- foreach(k = 1:no_data) %dopar% {
#     res <- Rsolnp::solnp(pars = c(1, 1, 1), fun = log_mar_lik_gp_der_oracle,
#                          LB = c(0.0001, 0.0001, 0.0001),
#                          UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
#                          control = list(TOL = 1e-5, trace = 0),
#                          y = YY[[k]]$y,
#                          x_vec = YY[[k]]$x, der_vec = cri_pts,
#                          H0 = H0_diff[[k]],
#                          is.sig.par = FALSE)
#     res$par}


# ===========
## DGP-single
# ===========
# no_data <- 5
system.time(StoEM_single_mixture_beta <- foreach(k = 1:no_data) %dopar% {
    stochastic_em_dgp(y = YY[[k]]$y, x = YY[[k]]$x, H0 = H0_diff[[k]], 
                      theta_init = c(1, 1, 1), epsilon = 1e-4, 
                      D = 100, a = 0, b = 2, n_sample = 1000, max_iter = 100,
                      lower = c(0.0001, 0.0001, 0.0001),
                      upper = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
                      shape1 = c(3, 10), shape2 = c(10, 3),
                      mixture_prob = c(0.5, 0.5),
                      ctrl = list(TOL = 1e-5, trace = 0),
                      is.sig.par = FALSE)})


system.time(StoEM_single_beta <- foreach(k = 1:no_data) %dopar% {
    stochastic_em_dgp(y = YY[[k]]$y, x = YY[[k]]$x, H0 = H0_diff[[k]], 
                      theta_init = c(1, 1, 1), epsilon = 1e-4, 
                      D = 100, a = 0, b = 2, n_sample = 1000, max_iter = 100,
                      lower = c(0.0001, 0.0001, 0.0001),
                      upper = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
                      shape1 = 2, shape2 = 2,
                      ctrl = list(TOL = 1e-5, trace = 0),
                      is.sig.par = FALSE)})


# =============
## DGP-multiple
# =============
# system.time(StoEM_multi <- foreach(k = 1:no_data) %dopar% {
#     stochastic_em_dgp_multi(y = YY[[k]]$y, x = YY[[k]]$x, H0 = H0_diff[[k]], 
#                             theta_init = c(1, 1, 1), epsilon = 1e-4, 
#                             D = 100, a_vec = c(0, 1), b_vec = c(1, 2),
#                             n_sample = 1000, max_iter = 100,
#                             lower = c(0.0001, 0.0001, 0.0001),
#                             upper = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
#                             ctrl = list(TOL = 1e-5, trace = 0),
#                             is.sig.par = FALSE)})
stopCluster(cl)




## prior of t is mixture of betas, e.g., p(t) = 0.5 beta(a1, b1) + 0.5 beta(a2, b2)
binary_idx <- rbinom(10000, 1, 0.5)
idx_1 <- which(binary_idx == 1)
idx_0 <- which(binary_idx == 0)
mixture_beta <- rbeta(10000, 2, 5)[which(binary_idx == 1)] + rbeta(10000, 5, 2)[which(binary_idx == 0)] 

beta_sample1_large <- rbeta(50000, 3, 10)
beta_sample2_large <- rbeta(50000, 10, 3)
beta_sample1_small <- rbeta(50000, 2, 6)
beta_sample2_small <- rbeta(50000, 6, 2)

# hist(beta_sample, breaks = 50)
# hist(beta_sample, breaks = 50)
mixture_beta_large <- c(beta_sample1_large, beta_sample2_large)
mixture_beta_small <- c(beta_sample1_small, beta_sample2_small)
# hist(mixture_beta, breaks = 50)

# rbeta(10000, 2, 5)

a <- 0
b <- 2

general_mixture_beta_large <- mixture_beta_large * (b - a) + a
general_mixture_beta_small <- mixture_beta_small * (b - a) + a
# hist(general_mixture_beta_small, breaks = 50)
d_large <- density(general_mixture_beta_large)
d_small <- density(general_mixture_beta_small)
d_beta <- density(rbeta(1000000, 2, 2) * (b - a) + a)
par(mfrow = c(1, 1))
plot(d_large, col = "red", lwd = 2, main = "Various priors", xlab = "x",
     xlim = c(0.06, 1.94))
lines(d_small, col = "blue", lwd = 2)
lines(d_beta, col = "purple", lwd = 2)
abline(h = 1 / 2, col = "green", lwd = 2)
legend("topright", c("0.6 * beta(20, 60; 0, 2) + 0.4 * beta(60, 20; 0, 2)",
                     "0.5 * beta(2, 6; 0, 2) + 0.5 * beta(6, 2; 0, 2)",
                     "Uniform(0, 2)", "beta(2, 2; 0, 2)"),
       col = c("red", "blue", "green", "purple"), lwd = c(2, 2, 2, 2), bty = "n")


par(mfrow = c(2, 3))
# d_beta <- density(rbeta(1000000, 2, 2) * (b - a) + a)
k = 8
plot(seq(0, 2, length.out = 100), rep(1 / 2, 100), col = "blue", type = "l", 
     lwd = 2, main = "uniform", xlab = "x", ylab = "", ylim = c(0, 1), 
     cex.main = 1.5)
plot(d_large, col = "red", lwd = 2, main = "0.5 * beta(3, 10) + 0.5 * beta(10, 3)",
     xlab = "x", ylab = "",ylim = c(0, 1), cex.main = 1.5,
     xlim = c(0.06, 1.94))
plot(d_beta, col = "green3", lwd = 2, xlim = c(0.06, 1.94), main = "beta(2, 2)",
     xlab = "x", ylab = "", ylim = c(0, 1), cex.main = 1.5)
hist_t(StoEM_single[[k]]$sample_t, 0, 2, ylim = c(0, 5), col = "blue", 
       den.line = FALSE, breaks = 60)
# title(list(paste("Dist. of t: Data", 1), cex = 1))
hist_t(StoEM_single_mixture_beta[[k]]$sample_t, breaks = 60, a = 0, b = 2, ylim = c(0, 5), 
       den.line = FALSE, col = "red")
hist_t(StoEM_single_beta[[k]]$sample_t, breaks = 60, a = 0, b = 2, ylim = c(0, 5),
       den.line = FALSE, col = "green3")


