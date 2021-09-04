################################################################################
# DGP-MCEM Demo                                                                #
# Cheng-Han Yu                                                                 #
################################################################################
# Load packages and R functions
pkg <- c("Matrix", "matrixcalc", "Rsolnp", "emulator", "R.utils")
lapply(pkg, require, character.only = TRUE)
sourceDirectory("./R")


# Load simulation data
load("./data/sim_data.Rdata", verbose = TRUE)


### set up
H0_diff <- lapply(YY, function(d) {
    outer(as.vector(d$x), as.vector(d$x),
          FUN = function(x1, x2) (x1 - x2))
})

idx <- seq(0, 2, length.out = 500)
x_test <- sort(seq(0, 2, length.out = 100))
n_test <- length(x_test)


## Take the first data and run the algorithm
fit_dgp <- mcem_dgp(y = YY[[1]]$y, x = YY[[1]]$x, 
                    H0 = H0_diff[[1]],
                    theta_init = c(1, 1), epsilon = 1e-3, 
                    D = 500, a = 0, b = 2, n_sample = 5000, 
                    max_iter = 100,
                    lower = c(0.001, 0.001),
                    upper = c(1/0.001, 1/0.001), init_t = 0.5, 
                    init_sig = 0.1,
                    shape1 = 1, shape2 = 1,
                    ctrl = list(TOL = 1e-5, trace = 0),
                    ga_shape = 1/2, ga_rate = 1/2,
                    a_h = 1, b_h = 1, 
                    mixture_prob = c(0.5, 0.5),
                    is.h.par = FALSE)


# Distribution of t
hist(fit_dgp$sample_t, breaks = 30, probability = TRUE)

# Distribution of sig
hist(fit_dgp$sample_sig, breaks = 30, probability = TRUE)


# hyperparameters
fit_dgp$thetas

tau <- fit_dgp$thetas[nrow(fit_dgp$thetas), 1]
h <- fit_dgp$thetas[nrow(fit_dgp$thetas), 2]

# sampling regression function f
sample_f <- get_pi_t_sig(yJ = c(YY[[1]]$y, rep(0, 1)),
                         x = YY[[1]]$x,
                         x_test = x_test,
                         idx_der = fit_dgp$sample_t,
                         sample_sig = fit_dgp$sample_sig,
                         tau = tau, h = h)


## plotting f
plot_pred_gp_f_y(x = YY[[1]]$x, y = YY[[1]]$y, idx = idx, x_test = x_test,
                 mu_test = sample_f$mu_test,
                 CI_Low_f = sample_f$ci_low,
                 CI_High_f = sample_f$ci_high,
                 xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                 ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                 plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
                 pred_lwd = 2, title = paste("DGP fitted f"), is.legend = FALSE)







