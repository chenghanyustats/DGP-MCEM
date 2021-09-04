################################################################################
# DGP Simulation Study: Compare Methods                                        #
# Cheng-Han Yu                                                                 #
################################################################################
rm(list = ls())


## load data sets
load("./Analysis/Data/SimCompareDataNew.Rdata", verbose = TRUE)
## change to load("./data/sim_data.rda", verbose = TRUE)

### set up
## small and large amp data have the same x
H0_diff_amp <- lapply(YY_small_amp, function(d) {
    outer(as.vector(d$x), as.vector(d$x),
          FUN = function(x1, x2) (x1 - x2))
})

idx <- seq(0, 2, length.out = 500)
x_test <- sort(seq(0, 2, length.out = 100))
n_test <- length(x_test)


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
EB_gp_tiny_amp <- foreach(k = 1:no_data) %dopar% {
    res <- Rsolnp::solnp(pars = c(1, 1, 1), fun = log_mar_lik_gp,
                         LB = c(0.0001, 0.0001, 0.0001),
                         UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
                         control = list(TOL = 1e-5, trace = 0),
                         y = YY_tiny_amp[[k]]$y, H0 = H0_diff_amp[[k]])
    res$par
}


EB_gp_small_amp <- foreach(k = 1:no_data) %dopar% {
    res <- Rsolnp::solnp(pars = c(1, 1, 1), fun = log_mar_lik_gp,
                         LB = c(0.0001, 0.0001, 0.0001),
                         UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
                         control = list(TOL = 1e-5, trace = 0),
                         y = YY_small_amp[[k]]$y, H0 = H0_diff_amp[[k]])
    res$par
}

EB_gp_large_amp <- foreach(k = 1:no_data) %dopar% {
    res <- Rsolnp::solnp(pars = c(1, 1, 1), fun = log_mar_lik_gp,
                         LB = c(0.0001, 0.0001, 0.0001),
                         UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
                         control = list(TOL = 1e-5, trace = 0),
                         y = YY_large_amp[[k]]$y, H0 = H0_diff_amp[[k]])
    res$par
}

# ===========
## DGP-oracle
# ===========
# EB_dgp_oracle_tiny_amp <- foreach(k = 1:no_data) %dopar% {
#     res <- Rsolnp::solnp(pars = c(1, 1, 1), fun = log_mar_lik_gp_der_oracle,
#                          LB = c(0.0001, 0.0001, 0.0001),
#                          UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
#                          control = list(TOL = 1e-5, trace = 0),
#                          y = YY_tiny_amp[[k]]$y,
#                          x_vec = YY_tiny_amp[[k]]$x, 
#                          der_vec = cri_pts_tiny_amp,
#                          H0 = H0_diff_amp[[k]],
#                          is.sig.par = FALSE)
#     res$par}
# 
# EB_dgp_oracle_small_amp <- foreach(k = 1:no_data) %dopar% {
#     res <- Rsolnp::solnp(pars = c(1, 1, 1), fun = log_mar_lik_gp_der_oracle,
#                          LB = c(0.0001, 0.0001, 0.0001),
#                          UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
#                          control = list(TOL = 1e-5, trace = 0),
#                          y = YY_small_amp[[k]]$y,
#                          x_vec = YY_small_amp[[k]]$x,
#                          der_vec = cri_pts_small_amp,
#                          H0 = H0_diff_amp[[k]],
#                          is.sig.par = FALSE)
#     res$par}


# EB_dgp_oracle_large_amp <- foreach(k = 1:no_data) %dopar% {
#     res <- Rsolnp::solnp(pars = c(1, 1, 1), fun = log_mar_lik_gp_der_oracle,
#                          LB = c(0.0001, 0.0001, 0.0001),
#                          UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
#                          control = list(TOL = 1e-5, trace = 0),
#                          y = YY_large_amp[[k]]$y,
#                          x_vec = YY_large_amp[[k]]$x,
#                          der_vec = cri_pts_large_amp,
#                          H0 = H0_diff_amp[[k]],
#                          is.sig.par = FALSE)
#     res$par}
# ===========
## DGP-single
# ===========
# no_data <- 5
system.time(StoEM_single_tiny_amp <- foreach(k = 1:no_data) %dopar% {
    stochastic_em_dgp(y = YY_tiny_amp[[k]]$y, x = YY_tiny_amp[[k]]$x,
                      H0 = H0_diff_amp[[k]],
                      theta_init = c(1, 1, 1), epsilon = 1e-4,
                      D = 100, a = 0, b = 2, n_sample = 1000, max_iter = 100,
                      lower = c(0.0001, 0.0001, 0.0001),
                      upper = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
                      shape1 = 1, shape2 = 1,
                      ctrl = list(TOL = 1e-5, trace = 0),
                      is.sig.par = FALSE)})

system.time(StoEM_single_small_amp <- foreach(k = 1:no_data) %dopar% {
    stochastic_em_dgp(y = YY_small_amp[[k]]$y, x = YY_small_amp[[k]]$x,
                      H0 = H0_diff_amp[[k]],
                      theta_init = c(1, 1, 1), epsilon = 1e-4,
                      D = 100, a = 0, b = 2, n_sample = 1000, max_iter = 100,
                      lower = c(0.0001, 0.0001, 0.0001),
                      upper = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
                      shape1 = 1, shape2 = 1,
                      ctrl = list(TOL = 1e-5, trace = 0),
                      is.sig.par = FALSE)})

system.time(StoEM_single_large_amp <- foreach(k = 1:no_data) %dopar% {
    stochastic_em_dgp(y = YY_large_amp[[k]]$y, x = YY_large_amp[[k]]$x,
                      H0 = H0_diff_amp[[k]],
                      theta_init = c(1, 1, 1), epsilon = 1e-4,
                      D = 100, a = 0, b = 2, n_sample = 1000, max_iter = 100,
                      lower = c(0.0001, 0.0001, 0.0001),
                      upper = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
                      shape1 = 1, shape2 = 1,
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

k <- 8
StoEM_multi_small_amp_8 <- stochastic_em_dgp_multi(y = YY_small_amp[[k]]$y, x = YY_small_amp[[k]]$x,
                                                 H0 = H0_diff_amp[[k]],
                            theta_init = c(1, 1, 1), epsilon = 1e-4,
                            D = 100, a_vec = c(0, 1), b_vec = c(1, 2),
                            n_sample = 1000, max_iter = 100,
                            lower = c(0.0001, 0.0001, 0.0001),
                            upper = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
                            ctrl = list(TOL = 1e-5, trace = 0),
                            is.sig.par = FALSE)




#############################################
## Extract sample of t and EB estimates from StoEM_single and StoEM_multiple
#############################################
sample_t_EM_single_lst_tiny_amp <- lapply(StoEM_single_tiny_amp, 
                                           function(x) { x$sample_t })
sample_t_EM_single_lst_small_amp <- lapply(StoEM_single_small_amp, 
                                           function(x) { x$sample_t })
sample_t_EM_single_lst_large_amp <- lapply(StoEM_single_large_amp,
                                           function(x) { x$sample_t })
# sample_t_EM_multi_lst <- lapply(StoEM_multi, function(x) { x$sample_t_mat })

EB_EM_single_lst_tiny_amp <- lapply(StoEM_single_tiny_amp, 
                                     function(x) { x$thetas[nrow(x$thetas), ] })
EB_EM_single_lst_small_amp <- lapply(StoEM_single_small_amp, 
                                     function(x) { x$thetas[nrow(x$thetas), ] })
EB_EM_single_lst_large_amp <- lapply(StoEM_single_large_amp,
                                     function(x) { x$thetas[nrow(x$thetas), ] })
# EB_EM_multi_lst <- lapply(StoEM_multi, function(x) { x$thetas[nrow(x$thetas), ] })


#############################################
## Predictive intervals from StoEM_single and StoEM_multi
#############################################

################################################################################
### Compare predictive f
################################################################################


## predictive f 
################################################################################
# ===========
## GPR no der
# ===========
pred_gp_lst_tiny_amp <- lapply(1:no_data, function(k) {
    get_pred_ci_gp(eb_par = EB_gp_tiny_amp[[k]], x = YY_tiny_amp[[k]]$x, 
                   x_test = x_test, y = YY_tiny_amp[[k]]$y)
})
pred_gp_lst_small_amp <- lapply(1:no_data, function(k) {
    get_pred_ci_gp(eb_par = EB_gp_small_amp[[k]], x = YY_small_amp[[k]]$x, 
                   x_test = x_test, y = YY_small_amp[[k]]$y)
})
pred_gp_lst_large_amp <- lapply(1:no_data, function(k) {
    get_pred_ci_gp(eb_par = EB_gp_large_amp[[k]], x = YY_large_amp[[k]]$x,
                   x_test = x_test, y = YY_large_amp[[k]]$y)
})
# ===========
## DGP-oracle
# ===========
# pred_dgp_oracle_lst_tiny_amp <- lapply(1:no_data, function(k) {
#     get_pred_ci_dgp_pt(yJ = c(YY_tiny_amp[[k]]$y, 0, 0),
#                        x = YY_tiny_amp[[k]]$x,
#                        x_test = x_test,
#                        idx_der = cri_pts_tiny_amp,
#                        sig = EB_dgp_oracle_tiny_amp[[k]][1],
#                        tau = EB_dgp_oracle_tiny_amp[[k]][2],
#                        h = EB_dgp_oracle_tiny_amp[[k]][3])
# })
# pred_dgp_oracle_lst_small_amp <- lapply(1:no_data, function(k) {
#     get_pred_ci_dgp_pt(yJ = c(YY_small_amp[[k]]$y, 0, 0),
#                        x = YY_small_amp[[k]]$x,
#                        x_test = x_test,
#                        idx_der = cri_pts_small_amp,
#                        sig = EB_dgp_oracle_small_amp[[k]][1],
#                        tau = EB_dgp_oracle_small_amp[[k]][2],
#                        h = EB_dgp_oracle_small_amp[[k]][3])
# })
# pred_dgp_oracle_lst_large_amp <- lapply(1:no_data, function(k) {
#     get_pred_ci_dgp_pt(yJ = c(YY_large_amp[[k]]$y, 0, 0),
#                        x = YY_large_amp[[k]]$x,
#                        x_test = x_test,
#                        idx_der = cri_pts_large_amp,
#                        sig = EB_dgp_oracle_large_amp[[k]][1],
#                        tau = EB_dgp_oracle_large_amp[[k]][2],
#                        h = EB_dgp_oracle_large_amp[[k]][3])
# })

# ===========
## DGP-multiple
# ===========
cl <- makeCluster(6, type = "FORK")
registerDoParallel(cl)

# system.time(pred_dgp_multi_lst <- foreach(k = 1:no_data) %dopar% {
#     get_pi_t(yJ = c(YY[[k]]$y, rep(0, 2)),
#              x = YY[[k]]$x, 
#              x_test = x_test, 
#              idx_der = sample_t_EM_multi_lst[[k]],
#              sig = EB_EM_multi_lst[[k]][1], 
#              tau = EB_EM_multi_lst[[k]][2],
#              h = EB_EM_multi_lst[[k]][3])
# })

# ===========
## DGP-single
# ===========
system.time(pred_dgp_single_lst_tiny_amp <- foreach(k = 1:no_data) %dopar% {
    get_pi_t(yJ = c(YY_tiny_amp[[k]]$y, rep(0, 1)),
             x = YY_tiny_amp[[k]]$x,
             x_test = x_test,
             idx_der = sample_t_EM_single_lst_tiny_amp[[k]],
             sig = EB_EM_single_lst_tiny_amp[[k]][1],
             tau = EB_EM_single_lst_tiny_amp[[k]][2],
             h = EB_EM_single_lst_tiny_amp[[k]][3])
})

system.time(pred_dgp_single_lst_small_amp <- foreach(k = 1:no_data) %dopar% {
    get_pi_t(yJ = c(YY_small_amp[[k]]$y, rep(0, 1)),
             x = YY_small_amp[[k]]$x,
             x_test = x_test,
             idx_der = sample_t_EM_single_lst_small_amp[[k]],
             sig = EB_EM_single_lst_small_amp[[k]][1],
             tau = EB_EM_single_lst_small_amp[[k]][2],
             h = EB_EM_single_lst_small_amp[[k]][3])
})
system.time(pred_dgp_single_lst_large_amp <- foreach(k = 1:no_data) %dopar% {
    get_pi_t(yJ = c(YY_large_amp[[k]]$y, rep(0, 1)),
             x = YY_large_amp[[k]]$x,
             x_test = x_test,
             idx_der = sample_t_EM_single_lst_large_amp[[k]],
             sig = EB_EM_single_lst_large_amp[[k]][1],
             tau = EB_EM_single_lst_large_amp[[k]][2],
             h = EB_EM_single_lst_large_amp[[k]][3])
})
stopCluster(cl)


save(EB_gp_tiny_amp, EB_gp_small_amp, EB_gp_large_amp,
     StoEM_single_tiny_amp, StoEM_single_small_amp, StoEM_single_large_amp,
     pred_dgp_single_lst_tiny_amp, pred_dgp_single_lst_small_amp, pred_dgp_single_lst_large_amp,
     pred_gp_lst_tiny_amp, pred_gp_lst_small_amp, pred_gp_lst_large_amp,
     file = "./simulation_amp.RData")

load("./simulation_amp.RData", verbose = TRUE)


pred_dgp_multi_lst_small_amp_8 <- get_pi_t(yJ = c(YY_small_amp[[k]]$y, rep(0, 2)),
         x = YY_small_amp[[k]]$x,
         x_test = x_test,
         idx_der = StoEM_multi_small_amp_8$sample_t_mat,
         sig = StoEM_multi_small_amp_8$thetas[5, 1],
         tau = StoEM_multi_small_amp_8$thetas[5, 2],
         h = StoEM_multi_small_amp_8$thetas[5, 3])

# =============================================================================
## Root Mean Square Error
# =============================================================================
true_fcn_val_tiny_amp <- regfcn_tiny_amp(x_test)
pred_mean_gp_tiny_amp <- lapply(pred_gp_lst_tiny_amp, function(x) x$mu_test)
pred_mean_dgp_oracle_tiny_amp <- lapply(pred_dgp_oracle_lst_tiny_amp, 
                                         function(x) x$mu_test)
# pred_mean_dgp_multi <- lapply(pred_dgp_multi_lst, function(x) x$mu_test)
pred_mean_dgp_single_tiny_amp <- lapply(pred_dgp_single_lst_tiny_amp, 
                                         function(x) x$mu_test)

rmse_gp_tiny_amp <- sapply(pred_mean_gp_tiny_amp, 
                            function(x) {rmse_f(x, true_fcn_val_tiny_amp)})
mean(rmse_gp_tiny_amp)

rmse_dgp_oracle_tiny_amp <- sapply(pred_mean_dgp_oracle_tiny_amp, function(x){
    rmse_f(x, true_fcn_val_tiny_amp)})
mean(rmse_dgp_oracle_tiny_amp)

# rmse_dgp_multi <- sapply(pred_mean_dgp_multi, function(x){
#     rmse_f(x, true_fcn_val)})
# mean(rmse_dgp_multi)

rmse_dgp_single_tiny_amp <- sapply(pred_mean_dgp_single_tiny_amp, function(x){
    rmse_f(x, true_fcn_val_tiny_amp)})
mean(rmse_dgp_single_tiny_amp)

# -----------------------------

true_fcn_val_small_amp <- regfcn_small_amp(x_test)
pred_mean_gp_small_amp <- lapply(pred_gp_lst_small_amp, function(x) x$mu_test)
pred_mean_dgp_oracle_small_amp <- lapply(pred_dgp_oracle_lst_small_amp, 
                                         function(x) x$mu_test)
# pred_mean_dgp_multi <- lapply(pred_dgp_multi_lst, function(x) x$mu_test)
pred_mean_dgp_single_small_amp <- lapply(pred_dgp_single_lst_small_amp, 
                                         function(x) x$mu_test)

rmse_gp_small_amp <- sapply(pred_mean_gp_small_amp, 
                            function(x) {rmse_f(x, true_fcn_val_small_amp)})
mean(rmse_gp_small_amp)

rmse_dgp_oracle_small_amp <- sapply(pred_mean_dgp_oracle_small_amp, function(x){
    rmse_f(x, true_fcn_val_small_amp)})
mean(rmse_dgp_oracle_small_amp)

# rmse_dgp_multi <- sapply(pred_mean_dgp_multi, function(x){
#     rmse_f(x, true_fcn_val)})
# mean(rmse_dgp_multi)

rmse_dgp_single_small_amp <- sapply(pred_mean_dgp_single_small_amp, function(x){
    rmse_f(x, true_fcn_val_small_amp)})
mean(rmse_dgp_single_small_amp)
# -----------------------------------------------------------
true_fcn_val_large_amp <- regfcn_large_amp(x_test)
pred_mean_gp_large_amp <- lapply(pred_gp_lst_large_amp, function(x) x$mu_test)
pred_mean_dgp_oracle_large_amp <- lapply(pred_dgp_oracle_lst_large_amp, 
                                         function(x) x$mu_test)
# pred_mean_dgp_multi <- lapply(pred_dgp_multi_lst, function(x) x$mu_test)
pred_mean_dgp_single_large_amp <- lapply(pred_dgp_single_lst_large_amp, 
                                         function(x) x$mu_test)

rmse_gp_large_amp <- sapply(pred_mean_gp_large_amp, 
                            function(x) {rmse_f(x, true_fcn_val_large_amp)})
mean(rmse_gp_large_amp)

rmse_dgp_oracle_large_amp <- sapply(pred_mean_dgp_oracle_large_amp, function(x){
    rmse_f(x, true_fcn_val_large_amp)})
mean(rmse_dgp_oracle_large_amp)

# rmse_dgp_multi <- sapply(pred_mean_dgp_multi, function(x){
#     rmse_f(x, true_fcn_val)})
# mean(rmse_dgp_multi)

rmse_dgp_single_large_amp <- sapply(pred_mean_dgp_single_large_amp, function(x){
    rmse_f(x, true_fcn_val_large_amp)})
mean(rmse_dgp_single_large_amp)


rmse_se_amp <- rbind(c(mean(rmse_gp_tiny_amp), 
                       mean(rmse_dgp_oracle_tiny_amp), 
                       mean(rmse_dgp_single_tiny_amp)),
                     c(mean(rmse_gp_small_amp), 
                          mean(rmse_dgp_oracle_small_amp), 
                          mean(rmse_dgp_single_small_amp)),
                        c(mean(rmse_gp), 
                          mean(rmse_dgp_oracle), 
                          mean(rmse_dgp_single)),
                        c(mean(rmse_gp_large_amp), 
                          mean(rmse_dgp_oracle_large_amp), 
                          mean(rmse_dgp_single_large_amp)))

rownames(rmse_se_amp) <- c("Tiny_amp", "Small_amp", "Medium_amp", "Large_amp")
colnames(rmse_se_amp) <- c("GPR", "Oracle", "DGP-single")
rmse_se_amp[-1, ]

################################################################################
### Compare distribution of t
################################################################################
# # ===============
# # GPR no der
# # ===============
# # generate M predictive curves and find the location of critical points 
# # for each f
# # (1) RMSE
# Deriv_1_arry <- array(0, dim = c(length(x_test) - 1, 1000, no_data))
# Deriv_2_arry <- array(0, dim = c(length(x_test) - 2, 1000, no_data))
# for (i in 1:no_data) {
#     Deriv_1_arry[, , i] <- apply(pred_gp_lst[[i]]$pred_f, 1, deriv_1st, 
#                                  x = x_test)
#     Deriv_2_arry[, , i] <- apply(pred_gp_lst[[i]]$pred_f, 1, deriv_2nd, 
#                                  x = x_test)
# }
# 
# zero_der_lst <- vector("list", length = no_data)
# for (k in 1:no_data) {
#     for (i in 1:1000) {
#         zero_der_lst[[k]][[i]] <- detect_zero_der(Deriv_1_arry[, i, k], x_test)
#     }
# }
# 
# zero_der_lst_1 <- lapply(zero_der_lst, function(x) unlist(x))
# 
# 
# # K-means clustering
# n_cl <- 1:4
# cl_kmeans_lst <- list()
# X1 <- unlist(zero_der_lst[[1]])
# for (k in n_cl) {
#     cl_kmeans_lst[[k]] <- stats::kmeans(X1, k, nstart = 50)
# }
# 
# par(mfrow = c(4, 1), mar = c(4, 1, 2, 1))
# 
# for (k in n_cl) {
#     plot(X1[cl_kmeans_lst[[k]]$cluster == 1], 
#          y = rep(1, length(X1[cl_kmeans_lst[[k]]$cluster == 1])), 
#          xlim = c(0, 2), col = 1, xlab = "t", ylab = "", axes = FALSE,
#          cex.lab = 1.5, cex.axis = 1.5)
#     axis(1, cex.axis = 1.5)
#     rug(X1, ticksize = 0.1, lwd = 0.2)
#     title(main = list(paste("k =", k), cex = 2))
#     for (j in 2:k) {
#         points(X1[cl_kmeans_lst[[k]]$cluster == j], 
#                y = rep(1, length(X1[cl_kmeans_lst[[k]]$cluster == j])), 
#                col = j)
#     }
#     points(cl_kmeans_lst[[k]]$centers, y = rep(1.1, k), col = 1:k, 
#            pch = 8, cex = 1)
# }
# par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
# wss <- unlist(lapply(cl_kmeans_lst, function(x) x$tot.withinss))
# plot(n_cl, wss, type = "b", col = "red",
#      xlab = "Number of clusters", ylab = "Within-cluster sum of squares",
#      cex.lab = 1.5)
# title(main = list("K-means scree plot", cex = 2))

# =============================================================================
# cl1_t_lst <- lapply(zero_der_lst_1, function(x){ x[x < 1] })
# 
# cl2_t_lst <- lapply(zero_der_lst_1_new, function(x){ x[x > 1] })
# 
# cl1_rmse_within_lst <- lapply(cl1_t_lst, function(x) {
#     sqrt(mean((x - cri_pts[1]) ^ 2))
# })
# cl2_rmse_within_lst <- lapply(cl2_t_lst, function(x) {
#     sqrt(mean((x - cri_pts[2]) ^ 2))
# })
# 
# cl1_sd_within_lst <- lapply(cl1_t_lst, function(x) sd(x))
# cl2_sd_within_lst <- lapply(cl2_t_lst, function(x) sd(x))
# 
# 
# cl1_mean_vec <- unlist(lapply(cl1_t_lst, mean))
# cl2_mean_vec <- unlist(lapply(cl2_t_lst, mean))
# 
# cl1_rmse_between <- mean((cl1_mean_vec - cri_pts[1])^2)
# cl2_rmse_between <- mean((cl2_mean_vec - cri_pts[2])^2)
# 
# cl1_sd_between <- sd(cl1_mean_vec)
# cl2_sd_between <- sd(cl2_mean_vec)
# 
# cl1_overall_mean <- mean(cl1_mean_vec)
# cl2_overall_mean <- mean(cl2_mean_vec)

# ===============
# DGP-mutiple
# ===============
# # (1) RMSE
# cl1_multi_t_lst <- lapply(sample_t_EM_multi_lst, function(x){ x[, 1] })
# cl2_multi_t_lst <- lapply(sample_t_EM_multi_lst, function(x){ x[, 2] })
# 
# cl1_rmse_within_multi_lst <- lapply(cl1_multi_t_lst, function(x) 
#     sqrt(mean((x - cri_pts[1]) ^ 2)))
# cl2_rmse_within_multi_lst <- lapply(cl2_multi_t_lst, function(x) 
#     sqrt(mean((x - cri_pts[2]) ^ 2)))
# 
# cl1_sd_within_multi_lst <- lapply(cl1_multi_t_lst, function(x) sd(x))
# cl2_sd_within_multi_lst <- lapply(cl2_multi_t_lst, function(x) sd(x))
# 
# cl1_mean_multi_vec <- unlist(lapply(cl1_multi_t_lst, mean))
# cl2_mean_multi_vec <- unlist(lapply(cl2_multi_t_lst, mean))
# 
# cl1_rmse_between_multi <- mean((cl1_mean_multi_vec - cri_pts[1]) ^ 2)
# cl2_rmse_between_multi <- mean((cl2_mean_multi_vec - cri_pts[2]) ^ 2)
# 
# cl1_sd_between_multi <- sd(cl1_mean_multi_vec)
# cl2_sd_between_multi <- sd(cl2_mean_multi_vec)
# 
# cl1_overall_mean_multi <- mean(cl1_mean_multi_vec)
# cl2_overall_mean_multi <- mean(cl2_mean_multi_vec)



# ===============
# DGP-single
# ===============
# # (1) RMSE
# cl1_single_t_lst <- lapply(sample_t_EM_single_lst, function(x){ x[x < 1] })
# cl2_single_t_lst <- lapply(sample_t_EM_single_lst, function(x){ x[x > 1] })
# 
# cl1_rmse_within_single_lst <- lapply(cl1_single_t_lst, function(x) { 
#     sqrt(mean((x - cri_pts[1]) ^ 2))
# })
# cl2_rmse_within_single_lst <- lapply(cl2_single_t_lst, function(x) {
#     sqrt(mean((x - cri_pts[2]) ^ 2))
# })
# 
# cl1_sd_within_single_lst <- lapply(cl1_single_t_lst, function(x) sd(x))
# cl2_sd_within_single_lst <- lapply(cl2_single_t_lst, function(x) sd(x))
# 
# cl1_mean_single_vec <- unlist(lapply(cl1_single_t_lst, mean))
# cl2_mean_single_vec <- unlist(lapply(cl2_single_t_lst, mean))
# 
# cl1_mode_single_vec <- unlist(lapply(cl1_single_t_lst, function(mm) {
#     density(mm)$x[which(density(mm)$y == max(density(mm)$y))]
# }))
# 
# cl2_mode_single_vec <- unlist(lapply(cl2_single_t_lst, function(mm) {
#     density(mm)$x[which(density(mm)$y == max(density(mm)$y))]
# }))
# 
# cl1_rmse_between_single <- mean((cl1_mean_single_vec - cri_pts[1]) ^ 2)
# cl2_rmse_between_single <- mean((cl2_mean_single_vec - cri_pts[2]) ^ 2)
# 
# cl1_rmse_between_single_mode <- mean((cl1_mode_single_vec - cri_pts[1]) ^ 2)
# cl2_rmse_between_single_mode <- mean((cl2_mode_single_vec - cri_pts[2]) ^ 2)
# 
# cl1_sd_between_single <- sd(cl1_mean_single_vec)
# cl2_sd_between_single <- sd(cl2_mean_single_vec)
# 
# cl1_overall_mean_single <- mean(cl1_mean_single_vec)
# cl2_overall_mean_single <- mean(cl2_mean_single_vec)


# ===============
# Plotting
# ===============
par(mfrow = c(2, 3))
for (k in 1:10) {
    hist_t(sample_t_EM_single_lst_tiny_amp[[k]], 0, 2, ylim = c(0, 10), 
           col = "blue", den.line = FALSE)
    title(list(paste("Dist. t DGP-single_tiny_amp Data", k), cex = 1.2))
    # hist_t(sample_t_EM_multi_lst_new[[k]], 0, 2, ylim = c(0, 10), col = "blue", 
    #        den.line = FALSE)
    # title(list(paste("Dist. of t of DGP-multiple: Data", k), cex = 2))
}


for (k in 1:10) {
    hist_t(sample_t_EM_single_lst_small_amp[[k]], 0, 2, ylim = c(0, 10), 
           col = "blue", den.line = FALSE)
    title(list(paste("Dist. t DGP-single_small_amp Data", k), cex = 1.2))
    # hist_t(sample_t_EM_multi_lst_new[[k]], 0, 2, ylim = c(0, 10), col = "blue", 
    #        den.line = FALSE)
    # title(list(paste("Dist. of t of DGP-multiple: Data", k), cex = 2))
}

par(mfrow = c(2, 3))
for (k in 1:10) {
    hist_t(sample_t_EM_single_lst_large_amp[[k]], 0, 2, ylim = c(0, 10),
           col = "blue", den.line = FALSE)
    title(list(paste("Dist. t DGP-single_large_amp Data", k), cex = 1.2))
    # hist_t(sample_t_EM_multi_lst_new[[k]], 0, 2, ylim = c(0, 10), col = "blue",
    #        den.line = FALSE)
    # title(list(paste("Dist. of t of DGP-multiple: Data", k), cex = 2))
}

# all_t <- unlist(sample_t_EM_single_lst)
# hist_t(all_t, 0, 2, ylim = c(0, 4), col = "blue", 
#        den.line = FALSE, cri_pt = cri_pts)
# title(list(paste("Dist. of t of all simulated data"), cex = 2))

# hist_t(sample_t_EM_multi_lst[[k]][, 1], 0, 2, ylim = c(0, 10), col = "blue", 
#        den.line = FALSE)
# hist_t(sample_t_EM_multi_lst[[k]][, 2], 0, 2, ylim = c(0, 10), col = "blue", 
#        den.line = FALSE)
# title(list(paste("Dist. of t of DGP-multiple: Data", k), cex = 2))


# =====================
### Predictive f
# =====================
par(mfcol = c(4, 5), mar = c(4, 4, 2, 1))
# pred_data1.png

for(k in 1:100) {
    # GPR no der 
    plot_pred_gp_f_y(x = YY_tiny_amp[[k]]$x, y = YY_tiny_amp[[k]]$y, 
                     idx = idx, x_test = x_test,
                     mu_test = pred_gp_lst_tiny_amp[[k]]$mu_test,
                     CI_Low_f = pred_gp_lst_tiny_amp[[k]]$ci_low,
                     CI_High_f = pred_gp_lst_tiny_amp[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, 
                     cri_pts = cri_pts_tiny_amp,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
                     pred_lwd = 2, title = paste("GPR_tiny_amp", k), 
                     is.legend = FALSE, true_fcn = regfcn_tiny_amp)
    plot_pred_gp_f_y(x = YY_tiny_amp[[k]]$x, y = YY_tiny_amp[[k]]$y, 
                     idx = idx, x_test = x_test,
                     mu_test = pred_dgp_single_lst_tiny_amp[[k]]$mu_test,
                     CI_Low_f = pred_dgp_single_lst_tiny_amp[[k]]$ci_low,
                     CI_High_f = pred_dgp_single_lst_tiny_amp[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, 
                     cri_pts = cri_pts_tiny_amp,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
                     pred_lwd = 2, title = paste("s-DGP_tiny_amp", k), 
                     is.legend = FALSE, true_fcn = regfcn_tiny_amp)
}


for(k in 1:100) {
    # GPR no der 
    plot_pred_gp_f_y(x = YY_small_amp[[k]]$x, y = YY_small_amp[[k]]$y, 
                     idx = idx, x_test = x_test,
                     mu_test = pred_gp_lst_small_amp[[k]]$mu_test,
                     CI_Low_f = pred_gp_lst_small_amp[[k]]$ci_low,
                     CI_High_f = pred_gp_lst_small_amp[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, 
                     cri_pts = cri_pts_small_amp,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
                     pred_lwd = 2, title = paste("GPR_small_amp", k), 
                     is.legend = FALSE, true_fcn = regfcn_small_amp)
    plot_pred_gp_f_y(x = YY_small_amp[[k]]$x, y = YY_small_amp[[k]]$y, 
                     idx = idx, x_test = x_test,
                     mu_test = pred_dgp_single_lst_small_amp[[k]]$mu_test,
                     CI_Low_f = pred_dgp_single_lst_small_amp[[k]]$ci_low,
                     CI_High_f = pred_dgp_single_lst_small_amp[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, 
                     cri_pts = cri_pts_small_amp,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
                     pred_lwd = 2, title = paste("s-DGP_s_amp", k), 
                     is.legend = FALSE, true_fcn = regfcn_small_amp)
}

# for(k in 1:100) {
#     # GPR no der 
#     plot_pred_gp_f_y(x = YY_large_amp[[k]]$x, y = YY_large_amp[[k]]$y, 
#                      idx = idx, x_test = x_test,
#                      mu_test = pred_gp_lst_large_amp[[k]]$mu_test,
#                      CI_Low_f = pred_gp_lst_large_amp[[k]]$ci_low,
#                      CI_High_f = pred_gp_lst_large_amp[[k]]$ci_high,
#                      xlim = c(0, 2), is.der.line = TRUE, 
#                      cri_pts = cri_pts_large_amp,
#                      ylim = c(0, 3), is.true.fcn = TRUE, cex = 1, 
#                      plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
#                      pred_lwd = 2, title = paste("GPR_large_amp", k), 
#                      is.legend = FALSE, true_fcn = regfcn_large_amp)
#     plot_pred_gp_f_y(x = YY_large_amp[[k]]$x, y = YY_large_amp[[k]]$y, 
#                      idx = idx, x_test = x_test,
#                      mu_test = pred_dgp_single_lst_large_amp[[k]]$mu_test,
#                      CI_Low_f = pred_dgp_single_lst_large_amp[[k]]$ci_low,
#                      CI_High_f = pred_dgp_single_lst_large_amp[[k]]$ci_high,
#                      xlim = c(0, 2), is.der.line = TRUE, 
#                      cri_pts = cri_pts_large_amp,
#                      ylim = c(0, 3), is.true.fcn = TRUE, cex = 1, 
#                      plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
#                      pred_lwd = 2, title = paste("s-DGP_l_amp", k), 
#                      is.legend = FALSE, true_fcn = regfcn_large_amp)
# }
# k = 8
par(mfrow = c(1,1))
for (k in 100) {
    plot_pred_gp_f_y(x = YY_small_amp[[k]]$x, y = YY_small_amp[[k]]$y, 
                     idx = idx, x_test = x_test,
                     mu_test = pred_gp_lst_small_amp[[k]]$mu_test,
                     CI_Low_f = pred_gp_lst_small_amp[[k]]$ci_low,
                     CI_High_f = pred_gp_lst_small_amp[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, 
                     cri_pts = cri_pts_small_amp,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 0, 1, 0.5), 
                     pred_lwd = 2, title = paste("GPR_small_amp", k), 
                     is.legend = FALSE, true_fcn = regfcn_small_amp)
}

plot_pred_gp_f_y(x = YY_small_amp[[k]]$x, y = YY_small_amp[[k]]$y, 
                 idx = idx, x_test = x_test,
                 mu_test = pred_gp_lst_small_amp[[k]]$mu_test,
                 CI_Low_f = pred_gp_lst_small_amp[[k]]$ci_low,
                 CI_High_f = pred_gp_lst_small_amp[[k]]$ci_high,
                 xlim = c(0, 2), is.der.line = TRUE, 
                 cri_pts = cri_pts_small_amp,
                 ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                 plot.type = "p", col.poly = rgb(0, 0, 1, 0.5), 
                 pred_lwd = 2, title = "Example of lower and upper bound", 
                 is.legend = FALSE, true_fcn = regfcn_small_amp)



gpr_no_cri_pt_idx_small <- c(8, 9, 19, 22, 27, 31, 32, 37, 38, 44,
                             62, 64, 65, 79, 80)

gpr_no_cri_pt_idx_tiny <- c(8, 9, 19, 22, 27, 31, 32, 37, 38, 44, 
                            62, 64, 65, 79, 80)

par(mfcol = c(2, 4), mar = c(4, 4, 2, 1))
plot_idx <- c(8, 32, 44, 65)
for (k in plot_idx) {
    matplot(x_test, t(pred_gp_lst_small_amp[[k]]$pred_f[1:10, ]), 
            type = "l", lwd = 2, main = paste("GPR_small_amp dataset", k), 
            xlab = "x", ylab = "y", ylim = c(0.3, 2.3))
    lines(x_test, pred_gp_lst_small_amp[[k]]$mu_test, 
          col = rgb(0, 0, 0, alpha = 0.6), lwd = 5)
    points(x = YY_small_amp[[k]]$x, y = YY_small_amp[[k]]$y, pch = 1, 
           col = rgb(0, 0, 0, alpha = 0.6), lwd = 2)
    abline(v = cri_pts_small_amp, lty = 2, cex = 0.1)
    
    matplot(x_test, pred_dgp_single_lst_small_amp[[k]]$pred_f[, 1:10], 
            type = "l", lwd = 2, main = paste("single-DGP_small_amp dataset", k), 
            xlab = "x", ylab = "y", ylim = c(0.3, 2.3))
    lines(x_test, pred_dgp_single_lst_small_amp[[k]]$mu_test, 
          col = rgb(0, 0, 0, alpha = 0.6), lwd = 5)
    points(x = YY_small_amp[[k]]$x, y = YY_small_amp[[k]]$y, pch = 1, 
           col = rgb(0, 0, 0, alpha = 0.6), lwd = 2)
    abline(v = cri_pts_small_amp, lty = 2, cex = 0.1)
}


k <- 8
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
# GPR no der 
plot_pred_gp_f_y(x = YY_small_amp[[k]]$x, y = YY_small_amp[[k]]$y, 
                 idx = idx, x_test = x_test,
                 mu_test = pred_gp_lst_small_amp[[k]]$mu_test,
                 CI_Low_f = pred_gp_lst_small_amp[[k]]$ci_low,
                 CI_High_f = pred_gp_lst_small_amp[[k]]$ci_high,
                 xlim = c(0, 2), is.der.line = TRUE, ylab = "f(x)",
                 cri_pts = cri_pts_small_amp, pred_col = "black",
                 ylim = c(0.3, 2.5), is.true.fcn = TRUE, cex = 0.5, 
                 plot.type = "p", col.poly = rgb(0, 0, 0, 0.2), 
                 pred_lwd = 2, title = paste("GPR"), 
                 is.legend = FALSE, true_fcn = regfcn_small_amp,
                 reg_fcn_col = "black", der_line_col = "black")


plot_pred_gp_f_y(x = YY_large_amp[[k]]$x, y = YY_large_amp[[k]]$y, 
                 idx = idx, x_test = x_test,
                 mu_test = pred_gp_lst_large_amp[[k]]$mu_test,
                 CI_Low_f = pred_gp_lst_large_amp[[k]]$ci_low,
                 CI_High_f = pred_gp_lst_large_amp[[k]]$ci_high,
                 xlim = c(0, 2), is.der.line = TRUE, 
                 cri_pts = cri_pts_large_amp,
                 ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                 plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
                 pred_lwd = 2, title = paste("GPR_large_amp"), 
                 is.legend = FALSE, true_fcn = regfcn_large_amp)

# DGP-oracle
plot_pred_gp_f_y(x = YY_small_amp[[k]]$x, y = YY_small_amp[[k]]$y, 
                 idx = idx, x_test = x_test,
                 mu_test = pred_dgp_oracle_lst_small_amp[[k]]$mu_test,
                 CI_Low_f = pred_dgp_oracle_lst_small_amp[[k]]$ci_low,
                 CI_High_f = pred_dgp_oracle_lst_small_amp[[k]]$ci_high,
                 xlim = c(0, 2), is.der.line = TRUE, 
                 cri_pts = cri_pts_small_amp,
                 ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                 plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
                 pred_lwd = 2, title = paste("oracle-DGP_small_amp"), 
                 is.legend = FALSE, legend.loc = "bottomleft",
                 true_fcn = regfcn_small_amp)
plot_pred_gp_f_y(x = YY_large_amp[[k]]$x, y = YY_large_amp[[k]]$y, 
                 idx = idx, x_test = x_test,
                 mu_test = pred_dgp_oracle_lst_large_amp[[k]]$mu_test,
                 CI_Low_f = pred_dgp_oracle_lst_large_amp[[k]]$ci_low,
                 CI_High_f = pred_dgp_oracle_lst_large_amp[[k]]$ci_high,
                 xlim = c(0, 2), is.der.line = TRUE, 
                 cri_pts = cri_pts_large_amp,
                 ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                 plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
                 pred_lwd = 2, title = paste("oracle-DGP_large_amp"), 
                 is.legend = FALSE, legend.loc = "bottomleft",
                 true_fcn = regfcn_large_amp)
# # DGP-multiple
# plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
#                  mu_test = pred_dgp_multi_lst[[k]]$mu_test,
#                  CI_Low_f = pred_dgp_multi_lst[[k]]$ci_low,
#                  CI_High_f = pred_dgp_multi_lst[[k]]$ci_high,
#                  xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
#                  ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
#                  plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
#                  pred_lwd = 2, title = paste("multiple-DGP"), 
#                  is.legend = FALSE)

# DGP-single
plot_pred_gp_f_y(x = YY_small_amp[[k]]$x, y = YY_small_amp[[k]]$y, 
                 idx = idx, x_test = x_test,
                 mu_test = pred_dgp_single_lst_small_amp[[k]]$mu_test,
                 CI_Low_f = pred_dgp_single_lst_small_amp[[k]]$ci_low,
                 CI_High_f = pred_dgp_single_lst_small_amp[[k]]$ci_high,
                 xlim = c(0, 2), is.der.line = TRUE, 
                 cri_pts = cri_pts_small_amp,
                 ylim = c(0.3, 2.5), is.true.fcn = TRUE, cex = 0.5, 
                 plot.type = "p", col.poly = rgb(0, 0, 0, 0.2), 
                 pred_lwd = 2, title = paste("single-DGP"), 
                 is.legend = FALSE, true_fcn = regfcn_small_amp,
                 reg_fcn_col = "black", der_line_col = "black", pred_col = "black")
plot_pred_gp_f_y(x = YY_large_amp[[k]]$x, y = YY_large_amp[[k]]$y, 
                 idx = idx, x_test = x_test,
                 mu_test = pred_dgp_single_lst_large_amp[[k]]$mu_test,
                 CI_Low_f = pred_dgp_single_lst_large_amp[[k]]$ci_low,
                 CI_High_f = pred_dgp_single_lst_large_amp[[k]]$ci_high,
                 xlim = c(0, 2), is.der.line = TRUE, 
                 cri_pts = cri_pts_large_amp,
                 ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                 plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
                 pred_lwd = 2, title = paste("single-DGP_large_amp"), 
                 is.legend = FALSE, true_fcn = regfcn_large_amp)


plot_pred_gp_f_y(x = YY_small_amp[[k]]$x, y = YY_small_amp[[k]]$y, 
                 idx = idx, x_test = x_test,
                 mu_test = pred_dgp_multi_lst_small_amp_8$mu_test,
                 CI_Low_f = pred_dgp_multi_lst_small_amp_8$ci_low,
                 CI_High_f = pred_dgp_multi_lst_small_amp_8$ci_high,
                 xlim = c(0, 2), is.der.line = TRUE, 
                 cri_pts = cri_pts_small_amp,
                 ylim = c(0.3, 2.5), is.true.fcn = TRUE, cex = 0.5, 
                 plot.type = "p", col.poly = rgb(0, 0, 0, 0.2), 
                 pred_lwd = 2, title = paste("multiple-DGP"), 
                 is.legend = FALSE, true_fcn = regfcn_small_amp,
                 reg_fcn_col = "black", der_line_col = "black", pred_col = "black")


# hist_t(sample_t_EM_multi_lst[[k]], 0, 2, ylim = c(0, 5), col = "blue", 
#        den.line = FALSE)
# title(list(paste("Distribution of t: multiple-DGP"), cex = 1.5))
# abline(v = cri_pts, col = "darkgreen", lty = 1, lwd = 1)
par(mfrow = c(1, 4))
for (k in plot_idx) {
    hist_t(sample_t_EM_single_lst_small_amp[[k]], 0, 2, ylim = c(0, 5), col = "blue", 
           den.line = FALSE)
    title(list(paste("Dist. of t: small_amp dataset", k), cex = 1.1))
    abline(v = cri_pts_small_amp, col = "darkgreen", lty = 1, lwd = 1)
}



hist_t(sample_t_EM_single_lst_large_amp[[k]], 0, 2, ylim = c(0, 5), col = "blue", 
       den.line = FALSE)
title(list(paste("Distribution of t: single-DGP_large_amp"), cex = 1.5))
abline(v = cri_pts_large_amp, col = "darkgreen", lty = 1, lwd = 1)

###### HPD interval ############################################################
# library(HDInterval)
# hdi(StoEM_one_lst[[100]]$sample_t, credMass = 0.95, allowSplit = TRUE)
# xxx <- density(StoEM_one_lst[[1]]$sample_t, from = 0, to = 2, bw = 0.01)
# plot(xxx)
# 
# CI_der_em_low <- matrix(0, 100, 2)
# CI_der_em_high <- matrix(0, 100, 2)
# for (k in 1:100) {
#     hdi_res <- hdi(density(StoEM_one_lst[[k]]$sample_t, from = 0, to = 2, 
#                            bw = 0.02), credMass = 0.95, allowSplit = TRUE)
#     CI_der_em_low[k, ] <- hdi_res[1, ]
#     CI_der_em_high[k, ] <- hdi_res[2, ]
# }
# 
# apply(CI_der_em_low, 2, mean)
# apply(CI_der_em_high, 2, mean)
# 
# CI_der_em <- hdi(density(StoEM_one_lst[[100]]$sample_t, 
#                          from = 0, to = 2, bw = 0.01), 
#                  credMass = 0.95, allowSplit = TRUE)
# 
# 
# hdi(density(StoEM_one_lst[[100]]$sample_t, from = 0, to = 2, bw = 0.01), 
#     credMass = 0.95, allowSplit = TRUE)[1, ]

##########################
hpd_interval_hist_lst_small_amp <- lapply(StoEM_single_small_amp, function(x) {
    get_hpd_interval_from_hist(samples = x$sample_t, breaks = 30, 
                               is.plot = FALSE)
})
hpd_interval_hist_lst_large_amp <- lapply(StoEM_single_large_amp, function(x) {
    get_hpd_interval_from_hist(samples = x$sample_t, breaks = 30, 
                               is.plot = FALSE)
})


par(mfcol = c(2, 3))
for (k in 1:3) {
    plot_pred_gp_f_y(x = YY_small_amp[[k]]$x, y = YY_small_amp[[k]]$y, 
                     idx = idx, x_test = x_test,
                     mu_test = pred_dgp_single_lst_small_amp[[k]]$mu_test,
                     CI_Low_f = pred_dgp_single_lst_small_amp[[k]]$ci_low,
                     CI_High_f = pred_dgp_single_lst_small_amp[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, 
                     cri_pts = cri_pts_small_amp,
                     ylim = c(0.3, 2.5), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 0, 1, 0.3), 
                     pred_lwd = 2, 
                     title = paste("Predictive f of DGP-single_small_amp: data", k), 
                     is.legend = TRUE, legend.loc = "bottomleft")
    hist_t(sample_t_EM_single_lst_small_amp[[k]], 0, 2, 
           ylim = c(0, 4), col = "blue",
           den.line = FALSE, cri_pt = cri_pts_small_amp)
    title(list(paste("Distribution of t_small_amp: data", k), cex = 1.8))
    segments(x0 = hpd_interval_hist_lst_small_amp[[k]]$ci_lower[[1]], 
             y0 = hpd_interval_hist_lst_small_amp[[k]]$den_value, 
             x1 = hpd_interval_hist_lst_small_amp[[k]]$ci_upper[[1]], 
             y1 = hpd_interval_hist_lst_small_amp[[k]]$den_value, 
             col = "red", lwd = 4)
    segments(x0 = hpd_interval_hist_lst_small_amp[[k]]$ci_lower[[2]], 
             y0 = hpd_interval_hist_lst_small_amp[[k]]$den_value, 
             x1 = hpd_interval_hist_lst_small_amp[[k]]$ci_upper[[2]], 
             y1 = hpd_interval_hist_lst_small_amp[[k]]$den_value, 
             col = "red", lwd = 4)
} 

par(mfcol = c(2, 3))
for (k in 1:3) {
    plot_pred_gp_f_y(x = YY_large_amp[[k]]$x, y = YY_large_amp[[k]]$y, 
                     idx = idx, x_test = x_test,
                     mu_test = pred_dgp_single_lst_large_amp[[k]]$mu_test,
                     CI_Low_f = pred_dgp_single_lst_large_amp[[k]]$ci_low,
                     CI_High_f = pred_dgp_single_lst_large_amp[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, 
                     cri_pts = cri_pts_large_amp,
                     ylim = c(0.3, 2.5), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 0, 1, 0.3), 
                     pred_lwd = 2, 
                     title = paste("Predictive f of DGP-single_large_amp: data", k), 
                     is.legend = TRUE, legend.loc = "bottomleft")
    hist_t(sample_t_EM_single_lst_large_amp[[k]], 0, 2, 
           ylim = c(0, 4), col = "blue",
           den.line = FALSE, cri_pt = cri_pts_large_amp)
    title(list(paste("Distribution of t_large_amp: data", k), cex = 1.8))
    segments(x0 = hpd_interval_hist_lst_large_amp[[k]]$ci_lower[[1]], 
             y0 = hpd_interval_hist_lst_large_amp[[k]]$den_value, 
             x1 = hpd_interval_hist_lst_large_amp[[k]]$ci_upper[[1]], 
             y1 = hpd_interval_hist_lst_large_amp[[k]]$den_value, 
             col = "red", lwd = 4)
    segments(x0 = hpd_interval_hist_lst_large_amp[[k]]$ci_lower[[2]], 
             y0 = hpd_interval_hist_lst_large_amp[[k]]$den_value, 
             x1 = hpd_interval_hist_lst_large_amp[[k]]$ci_upper[[2]], 
             y1 = hpd_interval_hist_lst_large_amp[[k]]$den_value, 
             col = "red", lwd = 4)
} 

# ## pred_dist_sim1.png
# par(mfrow = c(2, 2))
# for (k in c(1, 3)) {
#     plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
#                      mu_test = pred_EM_multi_lst[[k]]$mu_test,
#                      CI_Low_f = pred_EM_multi_lst[[k]]$ci_low,
#                      CI_High_f = pred_EM_multi_lst[[k]]$ci_high,
#                      xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
#                      ylim = c(0.3, 2.5), is.true.fcn = TRUE, cex = 1, 
#                      plot.type = "p", col.poly = rgb(0, 0, 1, 0.3), 
#                      pred_lwd = 2, 
#                      title = paste("Predictive f of DGP-multiple: data", k), 
#                      is.legend = FALSE, legend.loc = "bottomleft")
#     plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
#                      mu_test = pred_EM_single_lst[[k]]$mu_test,
#                      CI_Low_f = pred_EM_single_lst[[k]]$ci_low,
#                      CI_High_f = pred_EM_single_lst[[k]]$ci_high,
#                      xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
#                      ylim = c(0.3, 2.5), is.true.fcn = TRUE, cex = 1, 
#                      plot.type = "p", col.poly = rgb(0, 0, 1, 0.3), 
#                      pred_lwd = 2, 
#                      title = paste("Predictive f of DGP-single: data", k), 
#                      is.legend = FALSE, legend.loc = "bottomleft")
#     hist_t(sample_t_EM_multi_lst[[k]], 0, 2, ylim = c(0, 5), col = "blue",
#            den.line = FALSE, cri_pt = cri_pts)
#     title(list(paste("Dist. of t of DGP-multiple: data", k), cex = 1.6))
#     hist_t(sample_t_EM_single_lst[[k]], 0, 2, ylim = c(0, 5), col = "blue",
#            den.line = FALSE, cri_pt = cri_pts)
#     title(list(paste("Dist. of t of DGP-single: data", k), cex = 1.6))
# } 

# for (i in seq(0.01, 0.1, length.out = 10)) {
#     plot(density(unlist(sample_t_EM_single_lst), from = 0, to = 2, bw = i))
#     print(hdi(density(unlist(sample_t_EM_single_lst), from = 0, to = 2, bw = i), 
#               credMass = 0.95, allowSplit = TRUE))
# }
# 
# ## longer CI and larger bw
# bw_vec <- c("nrd0", "nrd", "ucv", "bcv", "SJ")
# par(mfrow = c(2, 3))
# for (i in 1:5) {
#     plot(density(unlist(sample_t_EM_single_lst), from = 0, to = 2, 
#                  bw = bw_vec[i]), main = paste("Method:", bw_vec[i]))
# }
# 
# ## shorter CI, smaller bw and close to the result from hist
# for (i in 1:5) {
#     hdi(density(StoEM_single_lst[[1]]$sample_t, from = 0, to = 2, 
#                 bw = bw_vec[i]), credMass = 0.95, allowSplit = TRUE)
# }


