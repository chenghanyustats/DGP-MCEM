################################################################################
# Reproduce Figures in the paper                                               #
# Cheng-Han Yu                                                                 #
################################################################################
pkg <- c("Matrix", "matrixcalc", "Rsolnp", "emulator", "R.utils")
lapply(pkg, require, character.only = TRUE)
sourceDirectory("./R")

# ==============================================================================
# Figure 1: Sample paths from a Gaussian process with (right) and without (left) 
# first order derivative information. A squared exponential kernel with 
# parameters $\tau=1$ and $h=1$ was used. Input values are indicated with 
# plus signs in the plots, and stationary points $t_1 = -4, t_2 = 0,$ and 
# $t_3 = 4$ are indicated with purple vertical lines in the right-hand plot.
# ==============================================================================
# ==============================================================================
# noise_free_path.png
# ===================
n_path <- 5
idx_obs <- sort(runif(5, -5, 5))
idx_der <- seq(-4, 4, length = 3)
idx_test <- seq(-5, 5, length = 101)
y_der <- rep(0, length(idx_der))
tau <- 1
phi <- 1
h <- sqrt(1 / phi)
sig <- 0
set.seed(10000)
# --------------
K_X <- compute_cov_1d(idx1 = idx_obs, tau = tau, h = h) + 
    diag(sig ^ 2, length(idx_obs))
K_XtestX <- compute_cov_1d(idx1 = idx_test, idx2 = idx_obs, tau = tau, h = h)
K_XtestXtest <- compute_cov_1d(idx1 = idx_test, idx2 = idx_test, tau = tau, 
                               h = h)
y_obs <- mvnfast::rmvn(1, rep(0, length(idx_obs)), K_X)

predict_mean_var <- cond_norm_mean_var(mean_vec_1 = 0, mean_vec_2 = 0, 
                                       obs_2 = y_obs[1, ],
                                       cov_mat_1 = K_XtestXtest, 
                                       cov_mat_2 = K_X, cov_mat_12 = K_XtestX)

predict_val <- mvnfast::rmvn(n_path, predict_mean_var$condMean,
                             predict_mean_var$condVar)
# --------------
y_joint <- c(y_obs[1, ], y_der)
K_joint <- compute_joint_cov(idx_obs = idx_obs, idx_der = idx_der, sig = sig, 
                             tau = tau, h = h, w = 0)
Kffnew <- compute_cov_1d(idx1 = idx_obs, idx2 = idx_test, tau = tau, h = h)
Kdfnew <- computeCovDer1(idx1 = idx_der, idx2 = idx_test, tau = tau, h = h) 
Kbar <- rbind(Kffnew, Kdfnew)
Kfnewfnew <- compute_cov_1d(idx1 = idx_test, idx2 = idx_test, tau = tau, h = h)

predict_mean_var_der <- cond_norm_mean_var(mean_vec_1 = 0, mean_vec_2 = 0,
                                           obs_2 = y_joint, 
                                           cov_mat_1 = Kfnewfnew,
                                           cov_mat_2 = K_joint,
                                           cov_mat_12 = t(Kbar))

predict_val_der <- mvnfast::rmvn(n_path, predict_mean_var_der$condMean, 
                                 predict_mean_var_der$condVar)
# --------------
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
gp_plot(idx = idx_obs, idx_test = idx_test, obs = y_obs[1, ], 
        predict_val = predict_val, ylim = c(-3, 3), is.multiplepath = TRUE,
        n_path = 5, is.title = FALSE, is.ci = FALSE)
title(substitute(paste("Paths without derivative: (", tau, ", ", h, ") = (", 
                       v2, ", ", v3, ")"), list(v2 = tau, v3 = phi)),
      cex.main = 1)

gp_plot(idx = idx_obs, idx_test = idx_test, obs = y_obs[1, ], 
        predict_val = predict_val_der, idx_der = idx_der, 
        ylim = c(-3, 3), is.multiplepath = TRUE, is.ci = FALSE,
        n_path = 5, is.title = FALSE, is.derline = TRUE)
title(substitute(
    paste("Paths with 1st derivative: (", tau, ", ", h, ") = (", v2, ", ", 
          v3, ")"), list(v2 = tau, v3 = phi)), cex.main = 1)

# ==============================================================================
# Figure 2: Simulated data. Root mean squared errors (RMSE) between the true 
# regression function and the estimated curves, as a function of $x$, averaged 
# across the 100 simulated datasets, for standard GPR, three different DGP models 
# (oracle, multiple and single), as described in the text, and the NKS method 
# of Song et al. (2006).The GPR and DGP methods use the test input $x_i^*$ 
# while NKS uses the original input $x_i$.
# ==============================================================================
# library(Matrix)
# library(matrixcalc)
# library(Rsolnp)
# library(emulator)

## load data sets
load("./data/sim_data.RData", verbose = TRUE)
# load("./simulation_main_result_sig.RData", verbose = TRUE)
load("./sim_result_algo.RData")
load("./sim_pred_gp_lst.RData")
load("./sim_pred_dgp_oracle_lst.RData")
load("./sim_pred_dgp_single_lst.RData")
load("./sim_pred_dgp_multiple_lst.RData")

x_test <- sort(seq(0, 2, length.out = 100))

true_fcn_val <- regfcn(x_test)
pred_mean_gp <- lapply(pred_no_der_lst_new, function(x) x$mu_test)
pred_mean_dgp_oracle <- lapply(pred_known_a_lst_new, function(x) x$mu_test)
pred_mean_dgp_multi <- lapply(pred_EM_multi_lst_new, function(x) x$mu_test)
pred_mean_dgp_single <- lapply(pred_EM_one_lst_new, function(x) x$mu_test)

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



# ==============================================================================
# Figure 3:Simulated data. {\it Top two rows:} Estimated curves for one simulated 
# dataset under different models: A standard GPR model and three different DGP 
# models (oracle, multiple and single), as described in the text. The true 
# regression function is shown as a solid line and the estimated curves as 
# dashed lines. Dots indicate input values. {\it Bottom row:} Posterior 
# distributions of $t$, for single and multiple DGPs, with vertical lines 
# indicating the locations of the true stationary points.
# ==============================================================================
# ==============================================================================

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



# ==============================================================================
# Figure 4: ERP Data Analysis: {\it Top left:} Illustration of amplitude and 
# latency of characteristic components of an ERP signal. {\it Top right}: 
#     Subject-level ERP waveforms, averaged over all trial conditions. 
# The continuous and dashed thick curves are the benchmark ERP averaged across 
# young and older subjects, respectively. The 0ms time, which corresponds to time 
# point 100, is the start of the onset of sound. Points 101 - 350 represent the 
# time the stimulus is played. The N100 component of interest is the amplitude of 
# the dip characterizing the signal in the time window between the vertical dashed
# lines. {\it Bottom plots}: ERP waveforms averaged across older and young subjects,
# for the two (voiced vs unvoiced) lexical biasing trial conditions of interest.
# ==============================================================================
# ==============================================================================
## amp_lat.png
par(mfrow = c(1, 1))
par(mar = c(2,2,2,1))
plot(seq(1, 2.5*pi, length.out = 1000), lwd = 4,
     cos(seq(-0, 2.5*pi, length.out = 1000) + pi / 2) , 
     type = "l", xlab= "", ylab = "",
     frame.plot = FALSE, ylim = c(-1.1, 1.1),
     col.axis = "white", col.ticks='white',
     cex.axis = 1, 
     axes = FALSE)
title(ylab = "mu V", 
      xlab = "Time (msec)", 
      line = 0.5, 
      cex.lab = 1.5, 
      family="Calibri Light")
abline(h = 0)
abline(v = 2)
arrows(x0 = 2, y0 = 1.02, x1 = 5, y1 = 1.02, lwd = 2, length = 0.1)
arrows(x0 = 5, y0 = 1.02, x1 = 2, y1 = 1.02, lwd = 2, length = 0.1)
arrows(x0 = 5.1, y0 = 0, x1 = 5.1, y1 = 0.95, lwd = 2, length = 0.1)
arrows(x0 = 5.1, y0 = 0.95, x1 = 5.1, y1 = 0, lwd = 2, length = 0.1)
text(3.5, 0.9, "latency", cex = 1.5)
text(6, 0.5, "amplitude", cex = 1.5)

# ------------------------------------------------------------------------------
# load erp data
load("./data/subj_erp_lst_group_mean.RData", verbose = TRUE)
load("./data/subj_y.RData", verbose = TRUE)
subj_erp_lst_group_mean_ts <- lapply(subj_erp_lst_group_mean, 
                                     function(x) apply(x, 1, mean)[11:360])
subj_erp_mat_group_mean_ts <- sapply(subj_erp_lst_group_mean, 
                                     function(x) apply(x, 1, mean)[11:360])
idx_11_good <- c(1, 5, 7, 9, 11, 12, 13, 14, 15, 17, 20)
subj_erp_group_mean_all_11 <- apply(subj_erp_mat_group_mean_ts[, idx_11_good], 
                                    1, mean)

load("./data/subj_erp_lst_group_mean_old_bias_only.RData", verbose = TRUE)
# subj_erp_lst_group_mean_old <- lapply(subj_data_lst_old,
#                                       function(x) apply(x, c(2, 3), mean))
subj_erp_lst_group_mean_old_ts <- lapply(subj_erp_lst_group_mean_old, 
                                         function(x) apply(x, 1, mean)[11:360])
subj_erp_mat_group_mean_old_ts <- sapply(subj_erp_lst_group_mean_old, 
                                     function(x) apply(x, 1, mean)[11:360])
subj_erp_group_mean_all_old <- apply(subj_erp_mat_group_mean_old_ts, 1, mean)
x_all <- seq(0.002, 0.7, by = 0.002)
# ------------------------------------------------------------------------------

### erp_data_young_old11.png
plot(x_all * 1000 - 200, subj_erp_lst_group_mean_ts[[idx_11_good[1]]], type = "l", 
     lty = 1, lwd = 0.5, ylim = c(-8, 3), ylab = "Voltage (mu V)",
     xlab = "Time (msec)", col = "black")
title(main = "All trials (Young and Older)", cex.main = 1.5)
for (i in 2:length(idx_11_good)) {
    lines(x_all * 1000 - 200, subj_erp_lst_group_mean_ts[[idx_11_good[i]]], 
          col = "black", lwd = 0.5)
}
lines(x_all * 1000 - 200, subj_erp_group_mean_all_11, lwd = 4, col = "black")
for (i in 1:11) {
    lines(x_all * 1000 - 200, subj_erp_lst_group_mean_old_ts[[i]], col = "black", 
          lwd = 0.7, lty = 2)
}
lines(x_all * 1000 - 200, subj_erp_group_mean_all_old, lwd = 4, col = "black",
      lty = 2)
# abline(h = 0)
abline(v = c(60, 130), lty = 2, lwd = 1.5)


# erp_bias_cond_old11.png
# midVOT_idx <- data_old_elec[1, 5, ] %in% c(20, 25, 30)
# data_old_elec_midVOT <- data_old_elec[, , midVOT_idx]
# get_bias_cond <- function(data, bias_cond = c(0, 1)) {
#     data_lst <- list()
#     for (i in 1:length(bias_cond)) {
#         idx_bias <- which(data[1, 6, ] == bias_cond[i])
#         data_i <- data[, , idx_bias]
#         data_lst[[i]] <- data_i
#     }
#     return(data_lst)
# }
# data_old_elec_bias <- get_bias_cond(data_old_elec_midVOT, bias_cond = c(-1, 1))
# data_old_bias_lst <- list()
# for (i in 1:length(bias_cond_vec)) {
#     print(i)
#     data_bias_i <- data_old_elec_bias[[i]]
#     
#     data_bias_i_avg_elect <- apply(data_bias_i, c(2, 3), mean)
#     data_bias_i_avg_elect_trl_ts <- apply(data_bias_i_avg_elect, 1, mean)[11:360]
#     data_old_bias_lst[[i]] <- data_bias_i_avg_elect_trl_ts
# }

load("./data/data_old_bias_lst.RData", verbose = TRUE)

plot_idx <- 125:225
plot((x_all[plot_idx] * 1000 - 200), data_old_bias_lst[[2]][plot_idx], 
     type = "l", lty = 1, ylim = c(-2.8, 1.4),
     ylab = "Voltage (mu V)", xlab = "Time (msec)", lwd = 3)
# title(main = ("ERPs by bias cond (avg elect, trl and subj)"))
title(main = ("Biasing conditions (Older)"), cex.main = 1.5)
lines((x_all[plot_idx] * 1000 - 200), data_old_bias_lst[[1]][plot_idx]+0.04, 
      col = 1, lwd = 3, lty = 2)
abline(h = 0)
# abline(v = c(60, 130), lwd = 1, lty = 1, col = "darkgreen")
legend("bottomright", c("Unvoiced Biasing", "Voiced Biasing"), lwd = 3, 
       lty = 1:2, bty = "n")

# erp_bias_cond_young11.png
load("./data/data_bias.RData", verbose = TRUE)

data_bias_mean_lst_11 <- lapply(data_bias_mat_lst, 
                                function(x) apply(x[, idx_11_good], 1, mean))
x_all_msec_50 <- sec_to_msec(x_all[plot_idx], pre_stimulus = 200)
plot(x_all_msec_50, data_bias_mean_lst_11[[2]][plot_idx], type = "l", lty = 1, 
     ylim = c(-2.8, 1.4), ylab = "Voltage (mu V)", xlab = "Time (msec)", lwd = 3)
title(main = ("Biasing conditions (Young)"), cex.main = 1.5)
lines(x_all_msec_50, data_bias_mean_lst_11[[1]][plot_idx] + 0.04, col = 1, 
      lwd = 3, lty = 2)
abline(h = 0)
legend("bottomleft", c("Unvoiced Biasing", "Voiced Biasing"), lwd = 3, 
       lty = 1:2, bty = "n")


# ==============================================================================
# Figure 5: ERP Data Analysis:  {\it Left column:} 95\% HPD regions of the 
# posterior distributions of $t$ for the older group, for all subjects and for 
# the voiced and unvoiced biasing trial conditions. {\it Right column:} 95\% HPDs 
# regions of the posterior distributions of $t$ for the young subjects group. In 
# all plots, the dash vertical lines indicate the posterior means obtained by 
# fitting a Gaussian mixture to the posterior samples of $t$,  averaged across 
# subjects, and the 95\% CI calculated as (mean $\pm$ 1.96 std) and shown as 
# shaded areas.
# ==============================================================================
# ==============================================================================
par(mfrow = c(1, 1))
par(mar = c(4, 4, 2, 1))

x_msec <- seq(50, 250, by = 2)

load("./erp_res_mcem.RData", verbose = TRUE)


# ------------------------------------------------------------------------------
mcmc_output_uniform_mcem <- res_uniform_mcem$mcmc_output
draws_msec_uniform_mcem <- apply(mcmc_output_uniform_mcem$draws, 2, function(a) {
    t_sam <- a * (max(x_all) - min(x_all)) + min(x_all)
    t_sam_msec <- sec_to_msec(t_sam, 200)
    return(t_sam_msec)})
mu_lst_young <- vector("list", 11)
var_lst_young <- vector("list", 11)
proportion_lst_young <- vector("list", 11)
library(Rmixmod)
for (s in 1:11) {
    data <- draws_msec_uniform_mcem[, s]
    out_test <- mixmodCluster(data, nbCluster = 2, criterion = "ICL")
    mu_lst_young[[s]] <- out_test@results[[1]]@parameters@mean
    var_lst_young[[s]] <- out_test@results[[1]]@parameters@variance
    proportion_lst_young[[s]] <- out_test@results[[1]]@parameters@proportions
}
# t1_mean_vec_young <- sapply(mu_lst_young, function(x) min(x))
# t2_mean_vec_young <- sapply(mu_lst_young, function(x) max(x))
# 
# t_mean_mat_young <- cbind(t1_mean_vec_young, t2_mean_vec_young)

order_prop <- lapply(proportion_lst_young, function(x) order(x, decreasing = TRUE))
t1_t2_mean_vec_young <- sapply(1:11, function(x) 
    mu_lst_young[[x]][order_prop[[x]]][1:2])
order_mean <- apply(t1_t2_mean_vec_young, 2, order)
t1_mean_vec_young <- apply(t1_t2_mean_vec_young, 2, min)
t2_mean_vec_young <- apply(t1_t2_mean_vec_young, 2, max)
t_mean_mat_young <- cbind(t1_mean_vec_young, t2_mean_vec_young)

r1_r2_mean_young <- apply(t_mean_mat_young, 2, mean)
r1_r2_sd_young <- apply(t_mean_mat_young, 2, sd)
names(r1_r2_mean_young) <- names(r1_r2_sd_young) <- c("r1", "r2")
hpdi_young_lst <- vector("list", 11)
for (k in 1:11) {
    hpdi_young <- get_hpd_interval_from_hist(samples = draws_msec_uniform_mcem[, k],
                                             breaks = 70,
                                             is.plot = FALSE)
    hpdi_young_lst[[k]][[1]] <- unlist(hpdi_young$ci_lower)
    hpdi_young_lst[[k]][[2]] <- unlist(hpdi_young$ci_upper)
}

plot_subj_hpd_bar(x_msec = x_msec, hpdi_lst = hpdi_young_lst, 
                  r1_r2_mean = r1_r2_mean_young, r1_r2_sd = r1_r2_sd_young, 
                  title = "Young")

# ------------------------------------------------------------------------------
mcmc_output_uniform_voiced_11 <- res_uniform_voiced_11$mcmc_output
draws_msec_uniform_voiced_11 <- apply(mcmc_output_uniform_voiced_11$draws, 2, function(a) {
    t_sam <- a * (max(x_all) - min(x_all)) + min(x_all)
    t_sam_msec <- sec_to_msec(t_sam, 200)
    return(t_sam_msec)})

mcmc_output_uniform_unvoiced_11 <- res_uniform_unvoiced_11$mcmc_output
draws_msec_uniform_unvoiced_11 <- apply(mcmc_output_uniform_unvoiced_11$draws, 2, function(a) {
    t_sam <- a * (max(x_all) - min(x_all)) + min(x_all)
    t_sam_msec <- sec_to_msec(t_sam, 200)
    return(t_sam_msec)})

mu_lst_voiced_11 <- vector("list", 11)
var_lst_voiced_11 <- vector("list", 11)
proportion_lst_voiced_11 <- vector("list", 11)
for (s in 1:11) {
    data <- draws_msec_uniform_voiced_11[, s]
    out_test <- mixmodCluster(data, nbCluster = 2, criterion = "ICL")
    mu_lst_voiced_11[[s]] <- out_test@results[[1]]@parameters@mean
    var_lst_voiced_11[[s]] <- out_test@results[[1]]@parameters@variance
}

t1_mean_vec_voiced_11 <- sapply(mu_lst_voiced_11, function(x) min(x))
t2_mean_vec_voiced_11 <- sapply(mu_lst_voiced_11, function(x) max(x))
t_mean_mat_voiced_11 <- cbind(t1_mean_vec_voiced_11, t2_mean_vec_voiced_11)
r1_r2_mean_voiced_11 <- apply(t_mean_mat_voiced_11, 2, mean)
r1_r2_sd_voiced_11 <- apply(t_mean_mat_voiced_11, 2, sd)
names(r1_r2_mean_voiced_11) <- names(r1_r2_sd_voiced_11) <- c("r1", "r2")

mu_lst_unvoiced_11 <- vector("list", 11)
var_lst_unvoiced_11 <- vector("list", 11)
proportion_lst_unvoiced_11 <- vector("list", 11)
par(mfrow = c(4, 3))
for (s in 1:11) {
    data <- draws_msec_uniform_unvoiced_11[, s]
    out_test <- mixmodCluster(data, nbCluster = 2, criterion = "ICL")
    mu_lst_unvoiced_11[[s]] <- out_test@results[[1]]@parameters@mean
    var_lst_unvoiced_11[[s]] <- out_test@results[[1]]@parameters@variance
    proportion_lst_unvoiced_11[[s]] <- out_test@results[[1]]@parameters@proportions
}

t1_mean_vec_unvoiced_11 <- sapply(mu_lst_unvoiced_11, function(x) min(x))
t2_mean_vec_unvoiced_11 <- sapply(mu_lst_unvoiced_11, function(x) max(x))
t_mean_mat_unvoiced_11 <- cbind(t1_mean_vec_unvoiced_11, t2_mean_vec_unvoiced_11)
r1_r2_mean_unvoiced_11 <- apply(t_mean_mat_unvoiced_11, 2, mean)
r1_r2_sd_unvoiced_11 <- apply(t_mean_mat_unvoiced_11, 2, sd)
names(r1_r2_mean_unvoiced_11) <- names(r1_r2_sd_unvoiced_11) <- c("r1", "r2")

hpdi_voiced_lst_11 <- vector("list", 11)
hpdi_unvoiced_lst_11 <- vector("list", 11)

for (k in 1:11) {
    hpdi_voiced_11 <- get_hpd_interval_from_hist(samples = draws_msec_uniform_voiced_11[, k],
                                                 breaks = 70,
                                                 is.plot = FALSE)
    hpdi_unvoiced_11 <- get_hpd_interval_from_hist(samples = draws_msec_uniform_unvoiced_11[, k],
                                                   breaks = 70,
                                                   is.plot = FALSE)
    hpdi_voiced_lst_11[[k]][[1]] <- unlist(hpdi_voiced_11$ci_lower)
    hpdi_voiced_lst_11[[k]][[2]] <- unlist(hpdi_voiced_11$ci_upper)
    hpdi_unvoiced_lst_11[[k]][[1]] <- unlist(hpdi_unvoiced_11$ci_lower)
    hpdi_unvoiced_lst_11[[k]][[2]] <- unlist(hpdi_unvoiced_11$ci_upper)
}

plot_subj_hpd_bar(x_msec = x_msec, hpdi_lst = hpdi_voiced_lst_11, 
                  r1_r2_mean = r1_r2_mean_voiced_11, 
                  r1_r2_sd = r1_r2_sd_voiced_11, 
                  title = "Voiced (Young)")

plot_subj_hpd_bar(x_msec = x_msec, hpdi_lst = hpdi_unvoiced_lst_11, 
                  r1_r2_mean = r1_r2_mean_unvoiced_11, 
                  r1_r2_sd = r1_r2_sd_unvoiced_11, 
                  title = "Unvoiced (Young)")
# ------------------------------------------------------------------------------

mcmc_output_uniform_old <- res_uniform_old$mcmc_output
draws_msec_uniform_old <- apply(mcmc_output_uniform_old$draws, 2, function(a) {
    t_sam <- a * (max(x_all) - min(x_all)) + min(x_all)
    t_sam_msec <- sec_to_msec(t_sam, 200)
    return(t_sam_msec)})


mu_lst_elderly <- vector("list", 11)
var_lst_elderly <- vector("list", 11)
proportion_lst_elderly <- vector("list", 11)
par(mfrow = c(4, 3))
for (s in 1:11) {
    data <- draws_msec_uniform_old[, s]
    out_test <- mixmodCluster(data, nbCluster = 2, )
    mu_lst_elderly[[s]] <- out_test@results[[1]]@parameters@mean
    var_lst_elderly[[s]] <- out_test@results[[1]]@parameters@variance
    proportion_lst_elderly[[s]] <- out_test@results[[1]]@parameters@proportions
}

t1_mean_vec_elderly <- sapply(mu_lst_elderly, function(x) min(x))
t2_mean_vec_elderly <- sapply(mu_lst_elderly, function(x) max(x))
t_mean_mat_elderly <- cbind(t1_mean_vec_elderly, t2_mean_vec_elderly)
r1_r2_mean_elderly <- apply(t_mean_mat_elderly, 2, mean)
r1_r2_sd_elderly <- apply(t_mean_mat_elderly, 2, sd)
names(r1_r2_mean_elderly) <- names(r1_r2_sd_elderly) <- c("r1", "r2")

hpdi_elderly_lst <- vector("list", 11)

for (k in 1:11) {
    hpdi_elderly <- get_hpd_interval_from_hist(samples = draws_msec_uniform_old[, k],
                                               breaks = 70,
                                               is.plot = FALSE)
    hpdi_elderly_lst[[k]][[1]] <- unlist(hpdi_elderly$ci_lower)
    hpdi_elderly_lst[[k]][[2]] <- unlist(hpdi_elderly$ci_upper)
}

plot_subj_hpd_bar(x_msec = x_msec, hpdi_lst = hpdi_elderly_lst, 
                  r1_r2_mean = r1_r2_mean_elderly, r1_r2_sd = r1_r2_sd_elderly, 
                  title = "Older")



# ------------------------------------------------------------------------------

mcmc_output_uniform_voiced_old <- res_uniform_voiced_old$mcmc_output
draws_msec_uniform_voiced_old <- apply(mcmc_output_uniform_voiced_old$draws, 2, function(a) {
    t_sam <- a * (max(x_all) - min(x_all)) + min(x_all)
    t_sam_msec <- sec_to_msec(t_sam, 200)
    return(t_sam_msec)})

mcmc_output_uniform_unvoiced_old <- res_uniform_unvoiced_old$mcmc_output
draws_msec_uniform_unvoiced_old <- apply(mcmc_output_uniform_unvoiced_old$draws, 2, function(a) {
    t_sam <- a * (max(x_all) - min(x_all)) + min(x_all)
    t_sam_msec <- sec_to_msec(t_sam, 200)
    return(t_sam_msec)})

mu_lst_voiced_old <- vector("list", 11)
var_lst_voiced_old <- vector("list", 11)
proportion_lst_voiced_old <- vector("list", 11)
par(mfrow = c(4, 3))
par(mar = c(4, 4, 2, 1))
for (s in 1:11) {
    data <- draws_msec_uniform_voiced_old[, s]
    out_test <- mixmodCluster(data, nbCluster = 2, criterion = "ICL")
    mu_lst_voiced_old[[s]] <- out_test@results[[1]]@parameters@mean
    var_lst_voiced_old[[s]] <- out_test@results[[1]]@parameters@variance
    proportion_lst_voiced_old[[s]] <- out_test@results[[1]]@parameters@proportions
}

t1_mean_vec_voiced_old <- sapply(mu_lst_voiced_old, function(x) min(x))
t2_mean_vec_voiced_old <- sapply(mu_lst_voiced_old, function(x) max(x))
t_mean_mat_voiced_old <- cbind(t1_mean_vec_voiced_old, t2_mean_vec_voiced_old)
r1_r2_mean_voiced_old <- apply(t_mean_mat_voiced_old, 2, mean)
r1_r2_sd_voiced_old <- apply(t_mean_mat_voiced_old, 2, sd)
names(r1_r2_mean_voiced_old) <- names(r1_r2_sd_voiced_old) <- c("r1", "r2")

# -----------------------
#### old-unvoiced
# -----------------------
mu_lst_unvoiced_old <- vector("list", 11)
var_lst_unvoiced_old <- vector("list", 11)
proportion_lst_unvoiced_old <- vector("list", 11)
par(mfrow = c(4, 3))
for (s in 1:11) {
    data <- draws_msec_uniform_unvoiced_old[, s]
    out_test <- mixmodCluster(data, nbCluster = 2, criterion = "ICL")
    mu_lst_unvoiced_old[[s]] <- out_test@results[[1]]@parameters@mean
    var_lst_unvoiced_old[[s]] <- out_test@results[[1]]@parameters@variance
    proportion_lst_unvoiced_old[[s]] <- out_test@results[[1]]@parameters@proportions
}

t1_mean_vec_unvoiced_old <- sapply(mu_lst_unvoiced_old, function(x) min(x))
t2_mean_vec_unvoiced_old <- sapply(mu_lst_unvoiced_old, function(x) max(x))
t_mean_mat_unvoiced_old <- cbind(t1_mean_vec_unvoiced_old, t2_mean_vec_unvoiced_old)
r1_r2_mean_unvoiced_old <- apply(t_mean_mat_unvoiced_old, 2, mean)
r1_r2_sd_unvoiced_old <- apply(t_mean_mat_unvoiced_old, 2, sd)
names(r1_r2_mean_unvoiced_old) <- names(r1_r2_sd_unvoiced_old) <- c("r1", "r2")


hpdi_voiced_lst_old <- vector("list", 11)
hpdi_unvoiced_lst_old <- vector("list", 11)

for (k in 1:11) {
    hpdi_voiced_old <- get_hpd_interval_from_hist(samples = draws_msec_uniform_voiced_old[, k],
                                                  breaks = 70,
                                                  is.plot = FALSE)
    hpdi_unvoiced_old <- get_hpd_interval_from_hist(samples = draws_msec_uniform_unvoiced_old[, k],
                                                    breaks = 70,
                                                    is.plot = FALSE)
    hpdi_voiced_lst_old[[k]][[1]] <- unlist(hpdi_voiced_old$ci_lower)
    hpdi_voiced_lst_old[[k]][[2]] <- unlist(hpdi_voiced_old$ci_upper)
    hpdi_unvoiced_lst_old[[k]][[1]] <- unlist(hpdi_unvoiced_old$ci_lower)
    hpdi_unvoiced_lst_old[[k]][[2]] <- unlist(hpdi_unvoiced_old$ci_upper)
}


plot_subj_hpd_bar(x_msec = x_msec, hpdi_lst = hpdi_voiced_old, 
                  r1_r2_mean = r1_r2_mean_voiced_old, 
                  r1_r2_sd = r1_r2_sd_voiced_old, 
                  title = "Voiced (Older)")

plot_subj_hpd_bar(x_msec = x_msec, hpdi_lst = hpdi_unvoiced_old, 
                  r1_r2_mean = r1_r2_mean_unvoiced_old, 
                  r1_r2_sd = r1_r2_sd_unvoiced_old, 
                  title = "Unvoiced (Older)")










































