################################################################################
# DGP Simulation Study: Compare Methods                                        #
# Cheng-Han Yu                                                                 #
################################################################################
pkg <- c("Matrix", "matrixcalc", "Rsolnp", "emulator", "R.utils")
lapply(pkg, require, character.only = TRUE)
sourceDirectory("./R")

## load data sets
# load("./Analysis/Data/SimCompareDataNew.Rdata", verbose = TRUE)
## change to load("./sim_data.Rdata", verbose = TRUE)
load("./data/sim_data.Rdata", verbose = TRUE)

### set up
H0_diff <- lapply(YY, function(d) {
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
EB_gp <- foreach(k = 1:no_data) %dopar% {
    res <- Rsolnp::solnp(pars = c(1, 1, 1), fun = log_mar_lik_gp,
                         LB = c(0.0001, 0.0001, 0.0001),
                         UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
                         control = list(TOL = 1e-5, trace = 0),
                         y = YY[[k]]$y, H0 = H0_diff[[k]])
    res$par
}

# ===========
## DGP-oracle
# ===========
EB_dgp_oracle <- foreach(k = 1:no_data) %dopar% {
    res <- Rsolnp::solnp(pars = c(1, 1, 1), fun = log_mar_lik_gp_der_oracle,
                         LB = c(0.0001, 0.0001, 0.0001),
                         UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
                         control = list(TOL = 1e-5, trace = 0),
                         y = YY[[k]]$y,
                         x_vec = YY[[k]]$x, der_vec = cri_pts,
                         H0 = H0_diff[[k]],
                         is.sig.par = FALSE)
    res$par}


# ===========
## DGP-single
# ===========
system.time(StoEM_single <- foreach(k = 1:100) %dopar% {
    mcem_dgp(y = YY[[k]]$y, x = YY[[k]]$x, 
             H0 = H0_diff[[k]],
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
             is.h.par = FALSE)})

# =============
## DGP-multiple
# =============
StoEM_multi <- foreach(k = 1:100) %dopar% {
    mcem_dgp_multi(y = YY[[k]]$y, x = YY[[k]]$x, 
                   H0 = H0_diff[[k]],
                   theta_init = c(1, 1), epsilon = 1e-3, 
                   D = 500, a_vec = c(0, 1), b_vec = c(1, 2), 
                   n_sample = 5000, 
                   max_iter = 100,
                   lower = c(0.001, 0.001),
                   upper = c(1/0.001, 1/0.001), 
                   init_t = c(0.5, 1.5), 
                   init_sig = 0.1,
                   shape1 = 1, shape2 = 1,
                   ctrl = list(TOL = 1e-5, trace = 0),
                   ga_shape = 1/2, ga_rate = 1/2,
                   a_h = 1, b_h = 1, is.h.par = FALSE)}


stopCluster(cl)




#############################################
## Extract sample of t and EB estimates from StoEM_single and StoEM_multiple
#############################################
sample_t_EM_single_lst <- lapply(StoEM_single, function(x) { x$sample_t })

sample_sig_EM_single_lst <- lapply(StoEM_single, function(x) { x$sample_sig })

sample_t_EM_multi_lst <- lapply(StoEM_multi, function(x) { x$sample_t_mat })

sample_sig_EM_multi_lst <- lapply(StoEM_multi, function(x) { x$sample_sig })

EB_EM_single_lst <- lapply(StoEM_single, function(x) { x$thetas[nrow(x$thetas), ] })

EB_EM_multi_lst <- lapply(StoEM_multi, function(x) { x$thetas[nrow(x$thetas), ] })





#############################################
## Predictive intervals 
#############################################

################################################################################
### Compare predictive f
################################################################################


## predictive f 
################################################################################
# ===========
## GPR no der
# ===========

pred_gp_lst <- lapply(1:no_data, function(k) {
    get_pred_ci_gp(eb_par = EB_gp[[k]], x = YY[[k]]$x, x_test = x_test, 
                   y = YY[[k]]$y)
})

# ===========
## DGP-oracle
# ===========
pred_dgp_oracle_lst <- lapply(1:no_data, function(k) {
    get_pred_ci_dgp_pt(yJ = c(YY[[k]]$y, 0, 0),
                       x = YY[[k]]$x, 
                       x_test = x_test, 
                       idx_der = cri_pts,
                       sig = EB_dgp_oracle[[k]][1], 
                       tau = EB_dgp_oracle[[k]][2],
                       h = EB_dgp_oracle[[k]][3])
})


# ===========
## DGP-multiple
# ===========

cl <- makeCluster(6, type = "FORK")
registerDoParallel(cl)

system.time(pred_dgp_multi_lst <- foreach(k = 1:no_data) %dopar% {
    get_pi_t_sig(yJ = c(YY[[k]]$y, rep(0, 2)),
                 x = YY[[k]]$x,
                 x_test = x_test,
                 idx_der = t(sample_t_EM_multi_lst[[k]]),
                 sample_sig = sample_sig_EM_multi_lst[[k]],
                 tau = EB_EM_multi_lst[[k]][1],
                 h = EB_EM_multi_lst[[k]][2])
})

# ===========
## DGP-single
# ===========
system.time(pred_dgp_single_lst <- foreach(k = 1:no_data) %dopar% {
    get_pi_t_sig(yJ = c(YY[[k]]$y, rep(0, 1)),
                 x = YY[[k]]$x,
                 x_test = x_test,
                 idx_der = sample_t_EM_single_lst[[k]],
                 sample_sig = sample_sig_EM_single_lst[[k]],
                 tau = EB_EM_single_lst[[k]][1],
                 h = EB_EM_single_lst[[k]][2])
})

stopCluster(cl)

save(EB_gp, EB_dgp_oracle,
     StoEM_multi, StoEM_single,
     pred_gp_lst, pred_dgp_oracle_lst,
     pred_dgp_multi_lst, pred_dgp_single_lst,
     file = "./simulation_main_result_sig.RData")


## load the results
load("./sim_result_algo.RData")
load("./sim_pred_gp_lst.RData")
load("./sim_pred_dgp_oracle_lst.RData")
load("./sim_pred_dgp_single_lst.RData")
load("./sim_pred_dgp_multiple_lst.RData")
# =============================================================================
## Root Mean Square Error
# =============================================================================
true_fcn_val <- regfcn(x_test)
pred_mean_gp <- lapply(pred_gp_lst, function(x) x$mu_test)
pred_mean_dgp_oracle <- lapply(pred_dgp_oracle_lst, function(x) x$mu_test)
pred_mean_dgp_multi <- lapply(pred_dgp_multi_lst, function(x) x$mu_test)
pred_mean_dgp_single <- lapply(pred_dgp_single_lst, function(x) x$mu_test)

rmse_gp <- sapply(pred_mean_gp, function(x) {rmse_f(x, true_fcn_val)})
mean(rmse_gp)

rmse_dgp_oracle <- sapply(pred_mean_dgp_oracle, function(x){
    rmse_f(x, true_fcn_val)})
mean(rmse_dgp_oracle)

rmse_dgp_multi <- sapply(pred_mean_dgp_multi, function(x){
    rmse_f(x, true_fcn_val)})
mean(rmse_dgp_multi)

rmse_dgp_single <- sapply(pred_mean_dgp_single, function(x){
    rmse_f(x, true_fcn_val)})
mean(rmse_dgp_single)

# ------
rmse_gp_each_f <- sapply(pred_gp_lst, function(x) {
    apply(t(x$pred_f), 2, function(s) sqrt(mean((s - true_fcn_val) ^ 2)))
})

mean(apply(rmse_gp_each_f, 2, mean))


rmse_dgp_oracle_each_f <- sapply(pred_dgp_oracle_lst, function(x) {
    apply(t(x$pred_f), 2, function(s) sqrt(mean((s - true_fcn_val) ^ 2)))
})

mean(apply(rmse_dgp_oracle_each_f, 2, mean))


rmse_dgp_single_each_f <- sapply(pred_dgp_single_lst, function(x) {
    apply(x$pred_f, 2, function(s) sqrt(mean((s - true_fcn_val) ^ 2)))
})

mean(apply(rmse_dgp_single_each_f, 2, mean))

# =============================================================================
## Root Mean Square Error a function of x and Plotting 
# =============================================================================
true_fcn_val <- regfcn(x_test)
# pred_mean_gp <- lapply(pred_no_der_lst_new, function(x) x$mu_test)
# pred_mean_dgp_oracle <- lapply(pred_known_a_lst_new, function(x) x$mu_test)
# pred_mean_dgp_multi <- lapply(pred_EM_multi_lst_new, function(x) x$mu_test)
# pred_mean_dgp_single <- lapply(pred_EM_one_lst_new, function(x) x$mu_test)

pred_mean_gp <- lapply(pred_gp_lst, function(x) x$mu_test)
pred_mean_dgp_oracle <- lapply(pred_dgp_oracle_lst, function(x) x$mu_test)
pred_mean_dgp_multi <- lapply(pred_dgp_multi_lst, function(x) x$mu_test)
pred_mean_dgp_single <- lapply(pred_dgp_single_lst, function(x) x$mu_test)



## GPR
diff_sq_gp <- sapply(pred_mean_gp, function(x) (x - true_fcn_val) ^ 2)
rmse_gp_x <- sqrt(apply(diff_sq_gp, 1, mean))


## Oracle-DGP
diff_sq_oracle <- sapply(pred_mean_dgp_oracle, function(x) (x - true_fcn_val) ^ 2)
rmse_dgp_oracle_x <- sqrt(apply(diff_sq_oracle, 1, mean))


## Multiple-DGP
diff_sq_multi <- sapply(pred_mean_dgp_multi, function(x) (x - true_fcn_val) ^ 2)
rmse_dgp_multi_x <- sqrt(apply(diff_sq_multi, 1, mean))

## Single-DGP
diff_sq_single <- sapply(pred_mean_dgp_single, function(x) (x - true_fcn_val) ^ 2)
rmse_dgp_single_x <- sqrt(apply(diff_sq_single, 1, mean))



## plotting
par(mar = c(4, 4, 2, 1))
par(mfrow = c(1, 1))
plot(x_test, rmse_gp_x, xlab = "x", ylab = "RMSE", main = "RMSE comparison", 
     col = 1, lty = 1, lwd = 2, type = "l", axes = F, ylim = c(0, 0.25))
lines(x_test, rmse_dgp_oracle_x, col = "purple", lty = 2, lwd = 2)
lines(x_test, rmse_dgp_multi_x, col = "red", lty = 3, lwd = 2)
lines(x_test, rmse_dgp_single_x, col = "blue", lty = 4, lwd = 2)
lines(x, rmse_song_x, col = "green3", lty = 5, lwd = 2)
legend("bottomright", bty = "n",
       legend = c("GPR", "Oracle-DGP", "Multiple-DGP", "Single-DGP", "NKS"),
       col = c("black", "purple", "red", "blue", "green3"),
       lty = c(1, 2, 3, 4, 5), 
       lwd = c(2, 2, 2, 2, 2))
axis(1)
axis(2, las = 1)







################################################################################
### Compare distribution of t
################################################################################
# ===============
# GPR no der
# ===============
# generate M predictive curves and find the location of critical points 
# for each f
# (1) RMSE
Deriv_1_arry <- array(0, dim = c(length(x_test) - 1, 1000, no_data))
Deriv_2_arry <- array(0, dim = c(length(x_test) - 2, 1000, no_data))
for (i in 1:no_data) {
    Deriv_1_arry[, , i] <- apply(pred_gp_lst[[i]]$pred_f, 1, deriv_1st, 
                                 x = x_test)
    Deriv_2_arry[, , i] <- apply(pred_gp_lst[[i]]$pred_f, 1, deriv_2nd, 
                                 x = x_test)
}

zero_der_lst <- vector("list", length = no_data)
for (k in 1:no_data) {
    for (i in 1:1000) {
        zero_der_lst[[k]][[i]] <- detect_zero_der(Deriv_1_arry[, i, k], x_test)
    }
}

zero_der_lst_1 <- lapply(zero_der_lst, function(x) unlist(x))


# K-means clustering
n_cl <- 1:4
cl_kmeans_lst <- list()
X1 <- unlist(zero_der_lst[[1]])
for (k in n_cl) {
    cl_kmeans_lst[[k]] <- stats::kmeans(X1, k, nstart = 50)
}

par(mfrow = c(4, 1), mar = c(4, 1, 2, 1))

for (k in n_cl) {
    plot(X1[cl_kmeans_lst[[k]]$cluster == 1], 
         y = rep(1, length(X1[cl_kmeans_lst[[k]]$cluster == 1])), 
         xlim = c(0, 2), col = 1, xlab = "t", ylab = "", axes = FALSE,
         cex.lab = 1.5, cex.axis = 1.5)
    axis(1, cex.axis = 1.5)
    rug(X1, ticksize = 0.1, lwd = 0.2)
    title(main = list(paste("k =", k), cex = 2))
    for (j in 2:k) {
        points(X1[cl_kmeans_lst[[k]]$cluster == j], 
               y = rep(1, length(X1[cl_kmeans_lst[[k]]$cluster == j])), 
               col = j)
    }
    points(cl_kmeans_lst[[k]]$centers, y = rep(1.1, k), col = 1:k, 
           pch = 8, cex = 1)
}
par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
wss <- unlist(lapply(cl_kmeans_lst, function(x) x$tot.withinss))
plot(n_cl, wss, type = "b", col = "red",
     xlab = "Number of clusters", ylab = "Within-cluster sum of squares",
     cex.lab = 1.5)
title(main = list("K-means scree plot", cex = 2))

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
# (1) RMSE
cl1_single_t_lst <- lapply(sample_t_EM_single_lst, function(x){ x[x < 1] })
cl2_single_t_lst <- lapply(sample_t_EM_single_lst, function(x){ x[x > 1] })

cl1_rmse_within_single_lst <- lapply(cl1_single_t_lst, function(x) { 
    sqrt(mean((x - cri_pts[1]) ^ 2))
})
cl2_rmse_within_single_lst <- lapply(cl2_single_t_lst, function(x) {
    sqrt(mean((x - cri_pts[2]) ^ 2))
})

cl1_sd_within_single_lst <- lapply(cl1_single_t_lst, function(x) sd(x))
cl2_sd_within_single_lst <- lapply(cl2_single_t_lst, function(x) sd(x))

cl1_mean_single_vec <- unlist(lapply(cl1_single_t_lst, mean))
cl2_mean_single_vec <- unlist(lapply(cl2_single_t_lst, mean))

cl1_mode_single_vec <- unlist(lapply(cl1_single_t_lst, function(mm) {
    density(mm)$x[which(density(mm)$y == max(density(mm)$y))]
}))

cl2_mode_single_vec <- unlist(lapply(cl2_single_t_lst, function(mm) {
    density(mm)$x[which(density(mm)$y == max(density(mm)$y))]
}))

(cl1_rmse_between_single <- sqrt(mean((cl1_mean_single_vec - cri_pts[1]) ^ 2)))
(cl2_rmse_between_single <- sqrt(mean((cl2_mean_single_vec - cri_pts[2]) ^ 2)))

(cl1_rmse_between_single_mode <- sqrt(mean((cl1_mode_single_vec - cri_pts[1]) ^ 2)))
(cl2_rmse_between_single_mode <- sqrt(mean((cl2_mode_single_vec - cri_pts[2]) ^ 2)))

cl1_sd_between_single <- sd(cl1_mean_single_vec)
cl2_sd_between_single <- sd(cl2_mean_single_vec)

cl1_overall_mean_single <- mean(cl1_mean_single_vec)
cl2_overall_mean_single <- mean(cl2_mean_single_vec)


# ===============
# Plotting
# ===============
par(mfrow = c(2, 3))
for (k in 1:100) {
    hist_t(sample_t_EM_single_lst[[k]], 0, 2, ylim = c(0, 10), col = "blue", 
           den.line = FALSE)
    title(list(paste("Dist. of t of DGP-single: Data", k), cex = 1.2))
}



# =====================
### Predictive f
# =====================
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
# pred_data1.png
k <- 2
# GPR no der 
plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                 mu_test = pred_gp_lst[[k]]$mu_test,
                 CI_Low_f = pred_gp_lst[[k]]$ci_low,
                 CI_High_f = pred_gp_lst[[k]]$ci_high,
                 xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                 ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                 plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
                 pred_lwd = 2, title = paste("GPR"), is.legend = FALSE)

# DGP-oracle
plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                 mu_test = pred_dgp_oracle_lst[[k]]$mu_test,
                 CI_Low_f = pred_dgp_oracle_lst[[k]]$ci_low,
                 CI_High_f = pred_dgp_oracle_lst[[k]]$ci_high,
                 xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                 ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                 plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
                 pred_lwd = 2, title = paste("oracle-DGP"), 
                 is.legend = FALSE, legend.loc = "bottomleft")

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
par(mfrow = c(6, 6), mar = c(4, 4, 2, 1))
for (k in 1:100) {
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_single_lst[[k]]$mu_test,
                     CI_Low_f = pred_dgp_single_lst[[k]]$ci_low,
                     CI_High_f = pred_dgp_single_lst[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = .6, 
                     plot.type = "p", col.poly = rgb(0, 0, 1, 0.4), 
                     pred_lwd = 2, title = paste("single-DGP-sig data", k), 
                     is.legend = FALSE)
    
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_EM_one_lst_new[[k]]$mu_test,
                     CI_Low_f = pred_EM_one_lst_new[[k]]$ci_low,
                     CI_High_f = pred_EM_one_lst_new[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = .6, 
                     plot.type = "p", col.poly = rgb(0, 0, 0.2, 0.4), 
                     pred_lwd = 2, 
                     title = paste("single-DGP data", k), 
                     is.legend = FALSE)
}


hist_t(sample_t_EM_single_lst[[k]], 0, 2, ylim = c(0, 5), col = "blue", 
       den.line = FALSE)
title(list(paste("Distribution of t: single-DGP"), cex = 1.5))
abline(v = cri_pts, col = "darkgreen", lty = 1, lwd = 1)

###### HPD interval ############################################################
library(HDInterval)
# hdi(StoEM_one_lst[[100]]$sample_t, credMass = 0.95, allowSplit = TRUE)
# xxx <- density(StoEM_one_lst[[1]]$sample_t, from = 0, to = 2, bw = 0.01)
# plot(xxx)

CI_der_em_low <- matrix(0, 100, 2)
CI_der_em_high <- matrix(0, 100, 2)

for (k in 1:100) {
    hdi_res <- hdi(density(StoEM_single[[k]]$sample_t, from = 0, to = 2, 
                           bw = 0.02), credMass = 0.95, allowSplit = TRUE)
    for (i in 1:nrow(hdi_res)) {
        if (((hdi_res[i, 1] - 0.1) < cri_pts[1]) & ((hdi_res[i, 2] + 0.1) > cri_pts[1])) {
            CI_der_em_low[k, 1] <- hdi_res[i, 1]
            CI_der_em_high[k, 1] <- hdi_res[i, 2]
            break
        }
    } 
    for (i in 1:nrow(hdi_res)) {
        if (((hdi_res[i, 1] - 0.1) < cri_pts[2]) & ((hdi_res[i, 2] + 0.1) > cri_pts[2])) {
            CI_der_em_low[k, 2] <- hdi_res[i, 1]
            CI_der_em_high[k, 2] <- hdi_res[i, 2]
            break
        }
    } 
}


apply(CI_der_em_low, 2, mean)
apply(CI_der_em_high, 2, mean)


##########################
hpd_interval_hist_lst <- lapply(StoEM_single, function(x) {
    get_hpd_interval_from_hist(samples = x$sample_t, breaks = 30, 
                               is.plot = FALSE)
})

mean(hpd_interval_hist_lst[[1]]$sample_cluster_lst[[1]])
mean(hpd_interval_hist_lst[[1]]$sample_cluster_lst[[2]])
hpd_interval_hist_lst[[1]]$ci_lower[[1]]
hpd_interval_hist_lst[[1]]$ci_lower[[2]]
hpd_interval_hist_lst[[1]]$ci_upper[[1]]
hpd_interval_hist_lst[[1]]$ci_upper[[2]]


par(mfcol = c(2, 3))
for (k in 1:3) {
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_single_lst[[k]]$mu_test,
                     CI_Low_f = pred_dgp_single_lst[[k]]$ci_low,
                     CI_High_f = pred_dgp_single_lst[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.5), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 0, 1, 0.3), 
                     pred_lwd = 2, 
                     title = paste("Predictive f of DGP-single: data", k), 
                     is.legend = TRUE, legend.loc = "bottomleft")
    hist_t(sample_t_EM_single_lst[[k]], 0, 2, ylim = c(0, 4), col = "blue",
           den.line = FALSE, cri_pt = cri_pts)
    title(list(paste("Distribution of t: data", k), cex = 1.8))
    segments(x0 = hpd_interval_hist_lst[[k]]$ci_lower[[1]], 
             y0 = hpd_interval_hist_lst[[k]]$den_value, 
             x1 = hpd_interval_hist_lst[[k]]$ci_upper[[1]], 
             y1 = hpd_interval_hist_lst[[k]]$den_value, col = "red", lwd = 4)
    segments(x0 = hpd_interval_hist_lst[[k]]$ci_lower[[2]], 
             y0 = hpd_interval_hist_lst[[k]]$den_value, 
             x1 = hpd_interval_hist_lst[[k]]$ci_upper[[2]], 
             y1 = hpd_interval_hist_lst[[k]]$den_value, col = "red", lwd = 4)
} 



par(mfrow = c(3, 2))
# pred_hist_black.png
for (k in 2) {
    # idx <- seq(min(YY[[k]]$x), max(YY[[k]]$x), length.out = 500)
    idx <- seq(0, 2, length.out = 500)
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_gp_lst[[k]]$mu_test,
                     CI_Low_f = pred_gp_lst[[k]]$ci_low,
                     CI_High_f = pred_gp_lst[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = .6, 
                     plot.type = "p", col.poly = rgb(0, 0, 0.2, 0.2), 
                     pred_lwd = 2,
                     reg_fcn_col = "black", pred_col = "black",
                     # title = paste("No der: data", k), 
                     title = paste("GPR"), 
                     der_line_col = "black",
                     is.legend = FALSE)
    
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_oracle_lst[[k]]$mu_test,
                     CI_Low_f = pred_dgp_oracle_lst[[k]]$ci_low,
                     CI_High_f = pred_dgp_oracle_lst[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = .6, 
                     plot.type = "p", col.poly = rgb(0, 0, 0.2, 0.2), 
                     pred_lwd = 2, 
                     reg_fcn_col = "black", pred_col = "black",
                     # title = paste("(a): data", k), 
                     title = paste("oracle-DGP"), 
                     der_line_col = "black",
                     is.legend = FALSE, legend.loc = "bottomleft")

    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_multi_lst[[k]]$mu_test,
                     CI_Low_f = pred_dgp_multi_lst[[k]]$ci_low,
                     CI_High_f = pred_dgp_multi_lst[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = .6, 
                     plot.type = "p", col.poly = rgb(0, 0, 0.2, 0.2), 
                     pred_lwd = 2, 
                     reg_fcn_col = "black", pred_col = "black",
                     # reg_fcn_col = "red", pred_col = rgb(1,0,0,0.8),
                     # title = paste("(b): data", k),
                     title = paste("multiple-DGP"), 
                     der_line_col = "black",
                     is.legend = FALSE)
    
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_single_lst[[k]]$mu_test,
                     CI_Low_f = pred_dgp_single_lst[[k]]$ci_low,
                     CI_High_f = pred_dgp_single_lst[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = .6, 
                     plot.type = "p", col.poly = rgb(0, 0, 0.2, 0.2), 
                     pred_lwd = 2, 
                     reg_fcn_col = "black", pred_col = "black",
                     # title = paste("(c): data", k), 
                     title = paste("single-DGP"), 
                     der_line_col = "black",
                     is.legend = FALSE)

    hist_t(sem_sig_test_mh_multi$sample_t_mat, 0, 2, ylim = c(0, 4), 
           col = "black", 
           den.line = FALSE)
    title(list(paste("Distribution of t: multiple-DGP"), cex = 1.5))
    abline(v = cri_pts, col = "black", lty = 1, lwd = 1)
    hist_t(sample_t_EM_single_lst[[k]], 0, 2, ylim = c(0, 4), col = "black", 
           den.line = FALSE)
    title(list(paste("Distribution of t: single-DGP"), cex = 1.5))
    abline(v = cri_pts, col = "black", lty = 1, lwd = 1)
}





## pred_dist_sim1.png
par(mfrow = c(2, 2))
for (k in c(1, 3)) {
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_EM_multi_lst[[k]]$mu_test,
                     CI_Low_f = pred_EM_multi_lst[[k]]$ci_low,
                     CI_High_f = pred_EM_multi_lst[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.5), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 0, 1, 0.3), 
                     pred_lwd = 2, 
                     title = paste("Predictive f of DGP-multiple: data", k), 
                     is.legend = FALSE, legend.loc = "bottomleft")
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_EM_single_lst[[k]]$mu_test,
                     CI_Low_f = pred_EM_single_lst[[k]]$ci_low,
                     CI_High_f = pred_EM_single_lst[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.5), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 0, 1, 0.3), 
                     pred_lwd = 2, 
                     title = paste("Predictive f of DGP-single: data", k), 
                     is.legend = FALSE, legend.loc = "bottomleft")
    hist_t(sample_t_EM_multi_lst[[k]], 0, 2, ylim = c(0, 5), col = "blue",
           den.line = FALSE, cri_pt = cri_pts)
    title(list(paste("Dist. of t of DGP-multiple: data", k), cex = 1.6))
    hist_t(sample_t_EM_single_lst[[k]], 0, 2, ylim = c(0, 5), col = "blue",
           den.line = FALSE, cri_pt = cri_pts)
    title(list(paste("Dist. of t of DGP-single: data", k), cex = 1.6))
} 

