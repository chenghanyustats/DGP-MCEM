pkg <- c("Matrix", "matrixcalc", "Rsolnp", "emulator", "R.utils")
lapply(pkg, require, character.only = TRUE)
sourceDirectory("./R")
load("./data/subj_y.RData", verbose = TRUE)
load("./data/y_bias.RData", verbose = TRUE)
load("./data/subj_y_old.RData", verbose = TRUE)
load("./data/y_bias_old_11_bias_only.RData", verbose = TRUE)
sub_data_idx <- 125:225
x_s_all <- (x_all - min(x_all)) / (max(x_all) - min(x_all))   ## standardize to (0, 1)
x <- x_all[sub_data_idx]
x_s <- x_s_all[sub_data_idx]
x_msec <- sec_to_msec(x, pre_stimulus = 200)
x_test <- sort(seq(min(x), max(x), length.out = 150))
x_test_s <- sort(seq(min(x_s), max(x_s), length.out = 150))

x_test_msec <- sec_to_msec(x_test, pre_stimulus = 200)
idx_11_good <- c(1, 5, 7, 9, 11, 12, 13, 14, 15, 17, 20)
shape1 <- shape2 <- 3
H0_mcmc <- outer(as.vector(x_s), as.vector(x_s), FUN = function(x1, x2) (x1 - x2))
# -----------------------
#### young
# -----------------------
## implementation of the algorithm
run_idx <- 11
listname <- c(paste("t.", 1:run_idx, sep = ""))
keep.lst <- sapply(rep(0, run_idx), as.list)
tune.lst <- sapply(rep(.1, run_idx), as.list)
names(keep.lst) <- listname
names(tune.lst) <- listname
name.par <- c(paste0("t", 1:run_idx), "sig")
start.lst <- list(t = rep(0.5, run_idx), sig = 0.5)
ga_shape = 200
ga_rate = 270


system.time(
    res_uniform_mcem <- sem_dgp_subj_sample_sig_mh(multi_y = matrix(subj_y[, idx_11_good],
                                                                    ncol = run_idx),
                                                   x = x_s, H0 = H0_mcmc, 
                                                   theta_init = c(1, 0.1), 
                                                   epsilon = 1e-4, D = 500, 
                                                   max_iter = 100,
                                                   a = min(x_s), b = max(x_s),
                                                   lower = c(0.0001, 0.0001),
                                                   upper = c(1/0.0001, 1/0.0001), 
                                                   shape1 = shape1, shape2 = shape2,
                                                   ctrl = list(TOL = 1e-5, trace = 0),
                                                   ga_shape = ga_shape, 
                                                   ga_rate = ga_rate,
                                                   n_mcmc = 11000, start.lst = start.lst,
                                                   burn_e = 1000, burn_final = 1000,
                                                   thin_e = 10, thin_final = 10,
                                                   name.par = name.par, adapt = TRUE,
                                                   tune.lst = tune.lst, keep.lst = keep.lst, 
                                                   tune.len = 20,
                                                   target.accept.rate = 0.35,
                                                   proposal_type = "uniform",
                                                   is.h.par = is.h.par, 
                                                   a_h = a_h, b_h = b_h))

eb_theta_uniform_mcem <- res_uniform_mcem$theta_mat[nrow(res_uniform_mcem$theta_mat), ]

mcmc_output_uniform_mcem <- res_uniform_mcem$mcmc_output
draws_msec_uniform_mcem <- apply(mcmc_output_uniform_mcem$draws, 2, function(a) {
    t_sam <- a * (max(x_all) - min(x_all)) + min(x_all)
    t_sam_msec <- sec_to_msec(t_sam, 200)
    return(t_sam_msec)})

mu_lst_young <- vector("list", 11)
var_lst_young <- vector("list", 11)
proportion_lst_young <- vector("list", 11)
par(mfrow = c(4, 3))
for (s in 1:11) {
    data <- draws_msec_uniform_mcem[, s]
    out_test <- mixmodCluster(data, nbCluster = 2, criterion = "ICL")
    mu_lst_young[[s]] <- out_test@results[[1]]@parameters@mean
    var_lst_young[[s]] <- out_test@results[[1]]@parameters@variance
    proportion_lst_young[[s]] <- out_test@results[[1]]@parameters@proportions
    hist_mix_gaussian(data = data, mu_vec = mu_lst_young[[s]], 
                      var_vec = var_lst_young[[s]], breaks = 70,
                      proportion_vec = proportion_lst_young[[s]], 
                      title = paste("Hist with GMM - young subj", s),
                      col = "blue")
}

t1_mean_vec_young <- sapply(mu_lst_young, function(x) min(x))
t2_mean_vec_young <- sapply(mu_lst_young, function(x) max(x))

t_mean_mat_young <- cbind(t1_mean_vec_young, t2_mean_vec_young)
min_idx_young <- sapply(mu_lst_young, function(x) which.min(x))
t1_var_vec_young <- sapply(1:11, function(x) var_lst_young[[x]][[min_idx_young[x]]])
t2_var_vec_young <- sapply(1:11, function(x) var_lst_young[[x]][[3-min_idx_young[x]]])

r1_r2_mean_young <- apply(t_mean_mat_young, 2, mean)
r1_r2_sd_young <- apply(t_mean_mat_young, 2, sd)
names(r1_r2_mean_young) <- names(r1_r2_sd_young) <- c("r1", "r2")

hpdi_young_lst <- vector("list", 11)
# hpdi_elderly_lst <- vector("list", 11)

for (k in 1:11) {
    hpdi_young <- get_hpd_interval_from_hist(samples = draws_msec_uniform[, k],
                                             breaks = 70,
                                             is.plot = FALSE)
    hpdi_young_lst[[k]][[1]] <- unlist(hpdi_young$ci_lower)
    hpdi_young_lst[[k]][[2]] <- unlist(hpdi_young$ci_upper)
}

first_der_all_11 <- deriv_1st(x_s, subj_erp_group_mean_all_11[sub_data_idx])
(zero_der_11 <- detect_zero_der(der_vec = first_der_all_11, 
                                x_test = seq(0.002, 0.7, by = 0.002)[sub_data_idx]))

par(mfrow = c(1, 1))
par(mar = c(4, 4, 2, 1))
plot(rep(10, 11), seq(11), type = "l", xlab = "Time (msec)", yaxt='n',
     ylab = "Subject", cex.axis = 1, tcl = -0.4, xlim = c(min(x_msec), max(x_msec)), 
     main = paste("Young"))
axis(2, at = 1:11, las = 2, tcl = -0.3, cex.axis = 0.9)
for (k in 1:11) {
    n_clu <- length(hpdi_young_lst[[k]][[1]])
    for (j in 1:n_clu) {
        segments(x0 = hpdi_young_lst[[k]][[1]][j], 
                 y0 = k, 
                 x1 = hpdi_young_lst[[k]][[2]][j], 
                 y1 = k, col = "black", lwd = 4)
    }
}
abline(v = r1_r2_mean_young, col = "black", lwd = 2, lty = 2)
r1_intvl_young <- c(r1_r2_mean_young[1] - 1.96 * r1_r2_sd_young[1],
                    r1_r2_mean_young[1] + 1.96 * r1_r2_sd_young[1])
r2_intvl_young <- c(r1_r2_mean_young[2] - 1.96 * r1_r2_sd_young[2],
                    r1_r2_mean_young[2] + 1.96 * r1_r2_sd_young[2])
polygon(c(seq(r1_intvl_young[1], r1_intvl_young[2], length = 100), 
          rev(seq(r1_intvl_young[1], r1_intvl_young[2], length = 100))), 
        c(rep(21, 100), rep(0, 100)), 
        col = rgb(0, 0, 0, 0.1), border = NA)
polygon(c(seq(r2_intvl_young[1], r2_intvl_young[2], length = 100), 
          rev(seq(r2_intvl_young[1], r2_intvl_young[2], length = 100))), 
        c(rep(21, 100), rep(0, 100)), 
        col = rgb(0, 0, 0, 0.1), border = NA)


# ------------------
## young-11 voiced
# -------------------
y_voiced_sd <- apply(y_voiced[, idx_11_good], 1, sd)
y_unvoiced_sd <- apply(y_unvoiced[, idx_11_good], 1, sd)
y_voiced_var <- y_voiced_sd ^ 2
y_unvoiced_var <- y_unvoiced_sd ^ 2
mean(y_voiced_var)  ## 1.70486  1.389491 1.344074
mean(y_unvoiced_var)  ## 1.941219  1.485492
var(y_voiced_var)
var(y_unvoiced_var)


## implementation of the algorithm
a_h = 1
b_h = 1
shape1 = 3
shape2 = 4
is.sig.par = TRUE
is.h.par = FALSE
run_idx <- 11
listname <- c(paste("t.", 1:run_idx, sep = ""))
keep.lst <- sapply(rep(0, run_idx), as.list)
tune.lst <- sapply(rep(.1, run_idx), as.list)
names(keep.lst) <- listname
names(tune.lst) <- listname
name.par <- c(paste0("t", 1:run_idx), "sig")
start.lst <- list(t = rep(0.5, run_idx), sig = 0.5)
ga_shape = 100
ga_rate = 140


system.time(
    res_uniform_voiced_11 <- sem_dgp_subj_sample_sig_mh(multi_y = matrix(y_voiced[, idx_11_good], 
                                                                              ncol = run_idx), 
                                                        x = x_s, H0 = H0_mcmc, 
                                                        theta_init = c(1, 0.1), 
                                                        epsilon = 1e-4, D = 500, max_iter = 100,
                                                        a = min(x_s), b = max(x_s),
                                                        lower = c(0.0001, 0.0001),
                                                        upper = c(1/0.0001, 1/0.0001), 
                                                        shape1 = 3, shape2 = 3,
                                                        ctrl = list(TOL = 1e-5, trace = 0),
                                                        ga_shape = ga_shape, ga_rate = ga_rate,
                                                        n_mcmc = 11000, start.lst = start.lst,
                                                        burn_e = 100, burn_final = 1000,
                                                        thin_e = 10, thin_final = 10,
                                                        name.par = name.par, adapt = TRUE,
                                                        tune.lst = tune.lst, keep.lst = keep.lst, 
                                                        tune.len = 20,
                                                        target.accept.rate = 0.35,
                                                        proposal_type = "uniform",
                                                        is.h.par = FALSE, 
                                                        a_h = 10, b_h = 10))
# ------------------
## young-11 unvoiced
# -------------------
system.time(
    res_uniform_unvoiced_11 <- sem_dgp_subj_sample_sig_mh(multi_y = matrix(y_unvoiced[, idx_11_good], 
                                                                           ncol = run_idx), 
                                                          x = x_s, H0 = H0_mcmc, 
                                                          theta_init = c(1, 0.1), 
                                                          epsilon = 1e-4, D = 500, max_iter = 100,
                                                          a = min(x_s), b = max(x_s),
                                                          lower = c(0.0001, 0.0001),
                                                          upper = c(1/0.0001, 1/0.0001), 
                                                          shape1 = 3, shape2 = 3,
                                                          ctrl = list(TOL = 1e-5, trace = 0),
                                                          ga_shape = 100, ga_rate = 140,
                                                          n_mcmc = 11000, start.lst = start.lst,
                                                          burn_e = 100, burn_final = 1000,
                                                          thin_e = 10, thin_final = 10,
                                                          name.par = name.par, adapt = TRUE,
                                                          tune.lst = tune.lst, keep.lst = keep.lst, 
                                                          tune.len = 20,
                                                          target.accept.rate = 0.35,
                                                          proposal_type = "uniform",
                                                          is.h.par = FALSE, 
                                                          a_h = 10, b_h = 10))

eb_theta_uniform_voiced_11 <- res_uniform_voiced_11$theta_mat[nrow(res_uniform_voiced_11$theta_mat), ]
mcmc_output_uniform_voiced_11 <- res_uniform_voiced_11$mcmc_output
draws_msec_uniform_voiced_11 <- apply(mcmc_output_uniform_voiced_11$draws, 2, function(a) {
    t_sam <- a * (max(x_all) - min(x_all)) + min(x_all)
    t_sam_msec <- sec_to_msec(t_sam, 200)
    return(t_sam_msec)})


eb_theta_uniform_unvoiced_11 <- res_uniform_unvoiced_11$theta_mat[nrow(res_uniform_unvoiced_11$theta_mat), ]
mcmc_output_uniform_unvoiced_11 <- res_uniform_unvoiced_11$mcmc_output
draws_msec_uniform_unvoiced_11 <- apply(mcmc_output_uniform_unvoiced_11$draws, 2, function(a) {
    t_sam <- a * (max(x_all) - min(x_all)) + min(x_all)
    t_sam_msec <- sec_to_msec(t_sam, 200)
    return(t_sam_msec)})


mu_lst_voiced_11 <- vector("list", 11)
var_lst_voiced_11 <- vector("list", 11)
proportion_lst_voiced_11 <- vector("list", 11)
par(mfrow = c(4, 3))
par(mar = c(4, 4, 2, 1))
for (s in 1:11) {
    data <- draws_msec_uniform_voiced_11[, s]
    out_test <- mixmodCluster(data, nbCluster = 2, criterion = "ICL")
    mu_lst_voiced_11[[s]] <- out_test@results[[1]]@parameters@mean
    var_lst_voiced_11[[s]] <- out_test@results[[1]]@parameters@variance
    proportion_lst_voiced_11[[s]] <- out_test@results[[1]]@parameters@proportions
    hist_mix_gaussian(data = data, mu_vec = mu_lst_voiced_11[[s]], 
                      var_vec = var_lst_voiced_11[[s]], breaks = 60,
                      proportion_vec = proportion_lst_voiced_11[[s]], 
                      title = paste("Hist with GMM - young11-voiced subj", s),
                      col = "blue")
}

t1_mean_vec_voiced_11 <- sapply(mu_lst_voiced_11, function(x) min(x))
t2_mean_vec_voiced_11 <- sapply(mu_lst_voiced_11, function(x) max(x))

t_mean_mat_voiced_11 <- cbind(t1_mean_vec_voiced_11, t2_mean_vec_voiced_11)
min_idx_voiced_11 <- sapply(mu_lst_voiced_11, function(x) which.min(x))
t1_var_vec_voiced_11 <- sapply(1:6, function(x) var_lst_voiced_11[[x]][[min_idx_voiced_11[x]]])
t2_var_vec_voiced_11 <- sapply(1:6, function(x) var_lst_voiced_11[[x]][[3-min_idx_voiced_11[x]]])

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
    hist_mix_gaussian(data = data, mu_vec = mu_lst_unvoiced_11[[s]], 
                      var_vec = var_lst_unvoiced_11[[s]], breaks = 60,
                      proportion_vec = proportion_lst_unvoiced_11[[s]], 
                      title = paste("Hist with GMM - young11-unvoiced subj", s),
                      col = "blue")
}

t1_mean_vec_unvoiced_11 <- sapply(mu_lst_unvoiced_11, function(x) min(x))
t2_mean_vec_unvoiced_11 <- sapply(mu_lst_unvoiced_11, function(x) max(x))

t_mean_mat_unvoiced_11 <- cbind(t1_mean_vec_unvoiced_11, t2_mean_vec_unvoiced_11)
min_idx_unvoiced_11 <- sapply(mu_lst_unvoiced_11, function(x) which.min(x))
t1_var_vec_unvoiced_11 <- sapply(1:6, function(x) var_lst_unvoiced_11[[x]][[min_idx_unvoiced_11[x]]])
t2_var_vec_unvoiced_11 <- sapply(1:6, function(x) var_lst_unvoiced_11[[x]][[3-min_idx_unvoiced_11[x]]])

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


par(mfrow = c(1, 1))
plot(rep(10, 11), seq(11), type = "l", xlab = "Time (msec)", yaxt='n',
     ylab = "Subject", cex.axis = 1, tcl = -0.4, xlim = c(min(x_msec), max(x_msec)), 
     main = paste("Voiced (Young)"))
axis(2, at = 1:11, las = 2, tcl = -0.3, cex.axis = 0.9)
for (k in 1:11) {
    n_clu <- length(hpdi_voiced_lst_11[[k]][[1]])
    for (j in 1:n_clu) {
        segments(x0 = hpdi_voiced_lst_11[[k]][[1]][j], 
                 y0 = k, 
                 x1 = hpdi_voiced_lst_11[[k]][[2]][j], 
                 y1 = k, col = "black", lwd = 4)
    }
}

abline(v = r1_r2_mean_voiced_11, col = "black", lwd = 2, lty = 2)
r1_intvl_voiced_11 <- c(r1_r2_mean_voiced_11[1] - 1.96 * r1_r2_sd_voiced_11[1],
                        r1_r2_mean_voiced_11[1] + 1.96 * r1_r2_sd_voiced_11[1])
r2_intvl_voiced_11 <- c(r1_r2_mean_voiced_11[2] - 1.96 * r1_r2_sd_voiced_11[2],
                        r1_r2_mean_voiced_11[2] + 1.96 * r1_r2_sd_voiced_11[2])
polygon(c(seq(r1_intvl_voiced_11[1], r1_intvl_voiced_11[2], length = 100), 
          rev(seq(r1_intvl_voiced_11[1], r1_intvl_voiced_11[2], length = 100))), 
        c(rep(21, 100), rep(0, 100)), 
        col = rgb(0, 0, 0, 0.1), border = NA)
polygon(c(seq(r2_intvl_voiced_11[1], r2_intvl_voiced_11[2], length = 100), 
          rev(seq(r2_intvl_voiced_11[1], r2_intvl_voiced_11[2], length = 100))), 
        c(rep(21, 100), rep(0, 100)), 
        col = rgb(0, 0, 0, 0.1), border = NA)


plot(rep(10, 11), seq(11), type = "l", xlab = "Time (msec)", yaxt='n',
     ylab = "Subject", cex.axis = 1, tcl = -0.4, xlim = c(min(x_msec), max(x_msec)), 
     main = paste("Unvoiced (Young)"))
axis(2, at = 1:11, las = 2, tcl = -0.3, cex.axis = 0.9)
for (k in 1:11) {
    n_clu <- length(hpdi_unvoiced_lst_11[[k]][[1]])
    for (j in 1:n_clu) {
        segments(x0 = hpdi_unvoiced_lst_11[[k]][[1]][j], 
                 y0 = k, 
                 x1 = hpdi_unvoiced_lst_11[[k]][[2]][j], 
                 y1 = k, col = "black", lwd = 4)
    }
}

abline(v = r1_r2_mean_unvoiced_11, col = "black", lwd = 2, lty = 2)
r1_intvl_unvoiced_11 <- c(r1_r2_mean_unvoiced_11[1] - 1.96 * r1_r2_sd_unvoiced_11[1],
                          r1_r2_mean_unvoiced_11[1] + 1.96 * r1_r2_sd_unvoiced_11[1])
r2_intvl_unvoiced_11 <- c(r1_r2_mean_unvoiced_11[2] - 1.96 * r1_r2_sd_unvoiced_11[2],
                          r1_r2_mean_unvoiced_11[2] + 1.96 * r1_r2_sd_unvoiced_11[2])
polygon(c(seq(r1_intvl_unvoiced_11[1], r1_intvl_unvoiced_11[2], length = 100), 
          rev(seq(r1_intvl_unvoiced_11[1], r1_intvl_unvoiced_11[2], length = 100))), 
        c(rep(21, 100), rep(0, 100)), 
        col = rgb(0, 0, 0, 0.1), border = NA)
polygon(c(seq(r2_intvl_unvoiced_11[1], r2_intvl_unvoiced_11[2], length = 100), 
          rev(seq(r2_intvl_unvoiced_11[1], r2_intvl_unvoiced_11[2], length = 100))), 
        c(rep(21, 100), rep(0, 100)), 
        col = rgb(0, 0, 0, 0.1), border = NA)







################################################################################
## Old people data
################################################################################
# subj_y_old <- subj_erp_mat_group_mean_old_ts[sub_data_idx, ]  ## data_aged.R

(subj_sd_old <- apply(subj_y_old, 1, sd))
(subj_var_old <- subj_sd_old ^ 2)

mean(subj_var_old)  
var(subj_var_old)

## implementation of the algorithm

run_idx <- 11
listname <- c(paste("t.", 1:run_idx, sep = ""))
keep.lst <- sapply(rep(0, run_idx), as.list)
tune.lst <- sapply(rep(.1, run_idx), as.list)
names(keep.lst) <- listname
names(tune.lst) <- listname
name.par <- c(paste0("t", 1:run_idx), "sig")
start.lst <- list(t = rep(0.5, run_idx), sig = 0.5)



ga_shape = 60
ga_rate = 141
a_h = 100
b_h = 100
shape1 = 2
shape2 = 3.5
is.sig.par = TRUE
is.h.par = FALSE
system.time(
    res_uniform_old <- sem_dgp_subj_sample_sig_mh(multi_y = matrix(subj_y_old,
                                                                   ncol = run_idx),
                                                       x = x_s, H0 = H0_mcmc, 
                                                       theta_init = c(1, 1), 
                                                       epsilon = 1e-4, D = 500, max_iter = 100,
                                                       a = min(x_s), b = max(x_s),
                                                       lower = c(0.0001, 0.0001),
                                                       upper = c(1/0.0001, 1/0.0001), 
                                                       shape1 = shape1, shape2 = shape2,
                                                       ctrl = list(TOL = 1e-5, trace = 0),
                                                       # ga_shape = 8.25, ga_rate = 9.0625,
                                                       # ga_shape = 3.538462, ga_rate = 2.53,
                                                       ga_shape = ga_shape, 
                                                       ga_rate = ga_rate,
                                                       n_mcmc = 11000, start.lst = start.lst,
                                                       burn_e = 1000, burn_final = 1000,
                                                       thin_e = 10, thin_final = 10,
                                                       name.par = name.par, adapt = TRUE,
                                                       tune.lst = tune.lst, keep.lst = keep.lst, 
                                                       tune.len = 20,
                                                       target.accept.rate = 0.35,
                                                       proposal_type = "uniform",
                                                       is.h.par = is.h.par, 
                                                       a_h = a_h, b_h = b_h))
# res_uniform_old <- res_uniform_old_test

eb_theta_uniform_old <- res_uniform_old$theta_mat[nrow(res_uniform_old$theta_mat), ]

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
    hist_mix_gaussian(data = data, mu_vec = mu_lst_elderly[[s]], 
                      var_vec = var_lst_elderly[[s]], 
                      proportion_vec = proportion_lst_elderly[[s]], 
                      title = paste("Hist with GMM - old subj", s),
                      col = "blue")
}

################
par(mfrow = c(2, 1))
s <- 5
data <- draws_msec_uniform_old[, s]
hist_mix_gaussian(data = data, mu_vec = mu_lst_elderly[[s]], 
                  var_vec = var_lst_elderly[[s]], 
                  proportion_vec = proportion_lst_elderly[[s]], 
                  title = "Single-DGP-Older t", xlab = "Time (msec)",
                  col = "gray", breaks = 70, ylim = c(0, 0.07))

hist(mcmc_output_uniform_old[[s]]$sample_sig, freq = FALSE, border = "white",
     main = "Single-DGP-Older sigma", xlab = "sigma", breaks = 25)


# par(mfrow = c(1, 1))
s <- 5
data <- as.vector(draws_msec_old[[s]])
hist_mix_gaussian(data = data, mu_vec = mu_lst_elderly[[s]], 
                  var_vec = var_lst_elderly[[s]], 
                  proportion_vec = proportion_lst_elderly[[s]], 
                  title = "Multiple-DGP-Older t", xlab = "Time (msec)",
                  col = "gray", breaks = 70, ylim = c(0, 0.07))

hist(multi_mcem_lst_3_3_old[[s]]$sample_sig, freq = FALSE, border = "white",
     main = "Multiple-DGP-Older sigma", xlab = "sigma", breaks = 25)


s <- 8
par(mfrow = c(1, 1))
data <- as.vector(draws_msec[[s]])
hist_mix_gaussian(data = data, mu_vec = mu_lst_young[[s]], 
                  var_vec = var_lst_young[[s]],
                  proportion_vec = proportion_lst_young[[s]], 
                  title = "Multiple-DGP-Young t",
                  col = "gray", breaks = 70)

hist(multi_mcem_lst_3_3_young[[s]]$sample_sig, freq = FALSE, border = "white",
     main = "Multiple-DGP-Young sigma", xlab = "sigma", breaks = 30)

################







t1_mean_vec_elderly <- sapply(mu_lst_elderly, function(x) min(x))
t2_mean_vec_elderly <- sapply(mu_lst_elderly, function(x) max(x))

t_mean_mat_elderly <- cbind(t1_mean_vec_elderly, t2_mean_vec_elderly)
min_idx <- sapply(mu_lst_elderly, function(x) which.min(x))
t1_var_vec_elderly <- sapply(1:6, function(x) var_lst_elderly[[x]][[min_idx[x]]])
t2_var_vec_elderly <- sapply(1:6, function(x) var_lst_elderly[[x]][[3-min_idx[x]]])

r1_r2_mean_elderly <- apply(t_mean_mat_elderly, 2, mean)
r1_r2_sd_elderly <- apply(t_mean_mat_elderly, 2, sd)
names(r1_r2_mean_elderly) <- names(r1_r2_sd_elderly) <- c("r1", "r2")


##########################################################
hpdi_elderly_lst <- vector("list", 11)

for (k in 1:11) {
    hpdi_elderly <- get_hpd_interval_from_hist(samples = draws_msec_uniform_old[, k],
                                                breaks = 70,
                                                is.plot = FALSE)
    hpdi_elderly_lst[[k]][[1]] <- unlist(hpdi_elderly$ci_lower)
    hpdi_elderly_lst[[k]][[2]] <- unlist(hpdi_elderly$ci_upper)
}



par(mfrow = c(1, 1))
par(mar=c(4, 4, 2, 1))
plot(rep(10, 11), seq(11), type = "l", xlab = "Time (msec)", yaxt='n',
     ylab = "Subject", cex.axis = 1, tcl = -0.4, xlim = c(min(x_msec), max(x_msec)), 
     main = paste("Older"))
axis(2, at = 1:11, las = 2, tcl = -0.3, cex.axis = 0.9)
for (k in 1:11) {
    n_clu <- length(hpdi_elderly_lst[[k]][[1]])
    for (j in 1:n_clu) {
        segments(x0 = hpdi_elderly_lst[[k]][[1]][j], 
                 y0 = k, 
                 x1 = hpdi_elderly_lst[[k]][[2]][j], 
                 y1 = k, col = "black", lwd = 4)
    }
}
# abline(v = zero_der_old[[1]] * 1000 - 200, col = "black", lwd = 2, lty = 1)
abline(v = r1_r2_mean_elderly, col = "black", lwd = 2, lty = 2)
r1_intvl_elderly <- c(r1_r2_mean_elderly[1] - 1.96 * r1_r2_sd_elderly[1],
                      r1_r2_mean_elderly[1] + 1.96 * r1_r2_sd_elderly[1])
r2_intvl_elderly <- c(r1_r2_mean_elderly[2] - 1.96 * r1_r2_sd_elderly[2],
                      r1_r2_mean_elderly[2] + 1.96 * r1_r2_sd_elderly[2])
polygon(c(seq(r1_intvl_elderly[1], r1_intvl_elderly[2], length = 100), 
          rev(seq(r1_intvl_elderly[1], r1_intvl_elderly[2], length = 100))), 
        c(rep(21, 100), rep(0, 100)), 
        col = rgb(0, 0, 0, 0.1), border = NA)
polygon(c(seq(r2_intvl_elderly[1], r2_intvl_elderly[2], length = 100), 
          rev(seq(r2_intvl_elderly[1], r2_intvl_elderly[2], length = 100))), 
        c(rep(21, 100), rep(0, 100)), 
        col = rgb(0, 0, 0, 0.1), border = NA)


# ============================
#### Voiced and Unvoiced
# ============================
# load("y_bias_old_11_bias_only.RData", verbose = TRUE)
y_voiced_sd_old <- apply(y_voiced_old, 1, sd)
y_unvoiced_sd_old <- apply(y_unvoiced_old, 1, sd)
y_voiced_var_old <- y_voiced_sd_old ^ 2
y_unvoiced_var_old <- y_unvoiced_sd_old ^ 2
mean(y_voiced_var_old)  #2.739144 2.958615
mean(y_unvoiced_var_old)  #0.7453595 2.820785
var(y_voiced_var_old)
var(y_unvoiced_var_old)

run_idx <- 11
listname <- c(paste("t.", 1:run_idx, sep = ""))
keep.lst <- sapply(rep(0, run_idx), as.list)
tune.lst <- sapply(rep(.1, run_idx), as.list)
names(keep.lst) <- listname
names(tune.lst) <- listname
name.par <- c(paste0("t", 1:run_idx), "sig")
start.lst <- list(t = rep(0.5, run_idx), sig = 0.5)

ga_shape = 100
ga_rate = 300
a_h = 1
b_h = 1
shape1 = 3
shape2 = 3.5
is.sig.par = TRUE
is.h.par = FALSE
system.time(
    res_uniform_voiced_old <- sem_dgp_subj_sample_sig_mh(multi_y = matrix(y_voiced_old[, 1:run_idx], 
                                                                               ncol = run_idx), 
                                                         x = x_s, H0 = H0_mcmc, 
                                                         theta_init = c(1, 0.1), 
                                                         epsilon = 1e-4, D = 500, max_iter = 100,
                                                         a = min(x_s), b = max(x_s),
                                                         lower = c(0.0001, 0.0001),
                                                         upper = c(1/0.0001, 1/0.0001), 
                                                         shape1 = 3, shape2 = 3,
                                                         ctrl = list(TOL = 1e-5, trace = 0),
                                                         ga_shape = ga_shape, ga_rate = ga_rate,
                                                         n_mcmc = 11000, start.lst = start.lst,
                                                         burn_e = 100, burn_final = 1000,
                                                         thin_e = 10, thin_final = 10,
                                                         name.par = name.par, adapt = TRUE,
                                                         tune.lst = tune.lst, keep.lst = keep.lst, 
                                                         tune.len = 20,
                                                         target.accept.rate = 0.35,
                                                         proposal_type = "uniform",
                                                         is.h.par = FALSE, 
                                                         a_h = 10, b_h = 10))

system.time(
    res_uniform_unvoiced_old <- sem_dgp_subj_sample_sig_mh(multi_y = matrix(y_unvoiced_old[, 1:run_idx], 
                                                                                 ncol = run_idx), 
                                                           x = x_s, H0 = H0_mcmc, 
                                                           theta_init = c(1, 0.1), 
                                                           epsilon = 1e-4, D = 500, max_iter = 100,
                                                           a = min(x_s), b = max(x_s),
                                                           lower = c(0.0001, 0.0001),
                                                           upper = c(1/0.0001, 1/0.0001), 
                                                           shape1 = 3, shape2 = 3,
                                                           ctrl = list(TOL = 1e-5, trace = 0),
                                                           ga_shape = 100, ga_rate = 300,
                                                           n_mcmc = 11000, start.lst = start.lst,
                                                           burn_e = 100, burn_final = 1000,
                                                           thin_e = 10, thin_final = 10,
                                                           name.par = name.par, adapt = TRUE,
                                                           tune.lst = tune.lst, keep.lst = keep.lst, 
                                                           tune.len = 20,
                                                           target.accept.rate = 0.35,
                                                           proposal_type = "uniform",
                                                           is.h.par = FALSE, 
                                                           a_h = 10, b_h = 10))

# first_der_all_voiced_old <- deriv_1st(x_s, data_bias_mean_lst_old[[1]][125:225])
# first_der_all_unvoiced_old <- deriv_1st(x_s, data_bias_mean_lst_old[[2]][125:225])
# (zero_der_voiced_old <- detect_zero_der(der_vec = first_der_all_voiced_old, 
#                                         x_test = seq(0.002, 0.7, by = 0.002)[sub_data_idx]))
# (zero_der_unvoiced_old <- detect_zero_der(der_vec = first_der_all_unvoiced_old, 
#                                           x_test = seq(0.002, 0.7, by = 0.002)[sub_data_idx]))

eb_theta_uniform_voiced_old <- res_uniform_voiced_old$theta_mat[nrow(res_uniform_voiced_old$theta_mat), ]

mcmc_output_uniform_voiced_old <- res_uniform_voiced_old$mcmc_output
draws_msec_uniform_voiced_old <- apply(mcmc_output_uniform_voiced_old$draws, 2, function(a) {
    t_sam <- a * (max(x_all) - min(x_all)) + min(x_all)
    t_sam_msec <- sec_to_msec(t_sam, 200)
    return(t_sam_msec)})

# par(mfrow = c(4, 3), mar = c(4, 4, 2, 1))
# for (i in 1:run_idx) {
#     # plot(draws_msec_uniform_voiced[, i], type = "l", main = paste(name.par[i], "Unif_voiced"))
#     hist(draws_msec_uniform_voiced_old[, i], breaks = 80, xlim = c(50, 250))
#     abline(v = zero_der_voiced_old[[1]][c(1, 2)] * 1000 - 200, col = "red")
# }
# hist(draws_msec_uniform_voiced, breaks = 150, xlim = c(50, 250), freq = FALSE)
# abline(v = zero_der_voiced[[1]] * 1000 - 200, col = "red")


eb_theta_uniform_unvoiced_old <- res_uniform_unvoiced_old$theta_mat[nrow(res_uniform_unvoiced_old$theta_mat), ]
mcmc_output_uniform_unvoiced_old <- res_uniform_unvoiced_old$mcmc_output
draws_msec_uniform_unvoiced_old <- apply(mcmc_output_uniform_unvoiced_old$draws, 2, function(a) {
    t_sam <- a * (max(x_all) - min(x_all)) + min(x_all)
    t_sam_msec <- sec_to_msec(t_sam, 200)
    return(t_sam_msec)})

# par(mfrow = c(4, 3), mar = c(4, 4, 2, 1))
# for (i in 1:run_idx) {
#     # plot(draws_msec_uniform_unvoiced[, i], type = "l", main = paste(name.par[i], "Unif_unvoiced"))
#     hist(draws_msec_uniform_unvoiced_old[, i], breaks = 80, xlim = c(50, 250))
#     abline(v = zero_der_unvoiced_old[[1]][c(1, 2)] * 1000 - 200, col = "red")
# }

# -----------------------
#### old-voiced
# -----------------------
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
    hist_mix_gaussian(data = data, mu_vec = mu_lst_voiced_old[[s]], 
                      var_vec = var_lst_voiced_old[[s]], breaks = 60,
                      proportion_vec = proportion_lst_voiced_old[[s]], 
                      title = paste("Hist with GMM - old-voiced subj", s),
                      col = "blue")
}

t1_mean_vec_voiced_old <- sapply(mu_lst_voiced_old, function(x) min(x))
t2_mean_vec_voiced_old <- sapply(mu_lst_voiced_old, function(x) max(x))

t_mean_mat_voiced_old <- cbind(t1_mean_vec_voiced_old, t2_mean_vec_voiced_old)
min_idx_voiced_old <- sapply(mu_lst_voiced_old, function(x) which.min(x))
t1_var_vec_voiced_old <- sapply(1:11, function(x) var_lst_voiced_old[[x]][[min_idx_voiced_old[x]]])
t2_var_vec_voiced_old <- sapply(1:11, function(x) var_lst_voiced_old[[x]][[3-min_idx_voiced_old[x]]])

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
    hist_mix_gaussian(data = data, mu_vec = mu_lst_unvoiced_old[[s]], 
                      var_vec = var_lst_unvoiced_old[[s]], breaks = 60,
                      proportion_vec = proportion_lst_unvoiced_old[[s]], 
                      title = paste("Hist with GMM - old-unvoiced subj", s),
                      col = "blue")
}

t1_mean_vec_unvoiced_old <- sapply(mu_lst_unvoiced_old, function(x) min(x))
t2_mean_vec_unvoiced_old <- sapply(mu_lst_unvoiced_old, function(x) max(x))

t_mean_mat_unvoiced_old <- cbind(t1_mean_vec_unvoiced_old, t2_mean_vec_unvoiced_old)
min_idx_unvoiced_old <- sapply(mu_lst_unvoiced_old, function(x) which.min(x))
t1_var_vec_unvoiced_old <- sapply(1:11, function(x) var_lst_unvoiced_old[[x]][[min_idx_unvoiced_old[x]]])
t2_var_vec_unvoiced_old <- sapply(1:11, function(x) var_lst_unvoiced_old[[x]][[3-min_idx_unvoiced_old[x]]])

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

par(mfrow = c(1, 1))
plot(rep(10, 11), seq(11), type = "l", xlab = "Time (msec)", yaxt='n',
     ylab = "Subject", cex.axis = 1, tcl = -0.4, xlim = c(min(x_msec), max(x_msec)), 
     main = paste("Voiced (Older)"))
axis(2, at = 1:11, las = 2, tcl = -0.3, cex.axis = 0.9)
for (k in 1:11) {
    n_clu <- length(hpdi_voiced_lst_old[[k]][[1]])
    for (j in 1:n_clu) {
        segments(x0 = hpdi_voiced_lst_old[[k]][[1]][j], 
                 y0 = k, 
                 x1 = hpdi_voiced_lst_old[[k]][[2]][j], 
                 y1 = k, col = "black", lwd = 4)
    }
}
# abline(v = zero_der_voiced_old[[1]][c(1, 2)] * 1000 - 200, col = "black", lwd = 2)
abline(v = r1_r2_mean_voiced_old, col = "black", lwd = 2, lty = 2)
r1_intvl_voiced_old <- c(r1_r2_mean_voiced_old[1] - 1.96 * r1_r2_sd_voiced_old[1],
                         r1_r2_mean_voiced_old[1] + 1.96 * r1_r2_sd_voiced_old[1])
r2_intvl_voiced_old <- c(r1_r2_mean_voiced_old[2] - 1.96 * r1_r2_sd_voiced_old[2],
                         r1_r2_mean_voiced_old[2] + 1.96 * r1_r2_sd_voiced_old[2])
polygon(c(seq(r1_intvl_voiced_old[1], r1_intvl_voiced_old[2], length = 100), 
          rev(seq(r1_intvl_voiced_old[1], r1_intvl_voiced_old[2], length = 100))), 
        c(rep(21, 100), rep(0, 100)), 
        col = rgb(0, 0, 0, 0.1), border = NA)
polygon(c(seq(r2_intvl_voiced_old[1], r2_intvl_voiced_old[2], length = 100), 
          rev(seq(r2_intvl_voiced_old[1], r2_intvl_voiced_old[2], length = 100))), 
        c(rep(21, 100), rep(0, 100)), 
        col = rgb(0, 0, 0, 0.1), border = NA)


plot(rep(10, 11), seq(11), type = "l", xlab = "Time (msec)", yaxt='n',
     ylab = "Subject", cex.axis = 1, tcl = -0.4, xlim = c(min(x_msec), max(x_msec)), 
     main = paste("Unvoiced (Older)"))
axis(2, at = 1:11, las = 2, tcl = -0.3, cex.axis = 0.9)
for (k in 1:11) {
    n_clu <- length(hpdi_unvoiced_lst_old[[k]][[1]])
    for (j in 1:n_clu) {
        segments(x0 = hpdi_unvoiced_lst_old[[k]][[1]][j], 
                 y0 = k, 
                 x1 = hpdi_unvoiced_lst_old[[k]][[2]][j], 
                 y1 = k, col = "black", lwd = 4)
    }
}
# abline(v = zero_der_unvoiced_old[[1]][c(1, 2)] * 1000 - 200, col = "black", lwd = 2)
abline(v = r1_r2_mean_unvoiced_old, col = "black", lwd = 2, lty = 2)
r1_intvl_unvoiced_old <- c(r1_r2_mean_unvoiced_old[1] - 1.96 * r1_r2_sd_unvoiced_old[1],
                           r1_r2_mean_unvoiced_old[1] + 1.96 * r1_r2_sd_unvoiced_old[1])
r2_intvl_unvoiced_old <- c(r1_r2_mean_unvoiced_old[2] - 1.96 * r1_r2_sd_unvoiced_old[2],
                           r1_r2_mean_unvoiced_old[2] + 1.96 * r1_r2_sd_unvoiced_old[2])
polygon(c(seq(r1_intvl_unvoiced_old[1], r1_intvl_unvoiced_old[2], length = 100), 
          rev(seq(r1_intvl_unvoiced_old[1], r1_intvl_unvoiced_old[2], length = 100))), 
        c(rep(21, 100), rep(0, 100)), 
        col = rgb(0, 0, 0, 0.1), border = NA)
polygon(c(seq(r2_intvl_unvoiced_old[1], r2_intvl_unvoiced_old[2], length = 100), 
          rev(seq(r2_intvl_unvoiced_old[1], r2_intvl_unvoiced_old[2], length = 100))), 
        c(rep(21, 100), rep(0, 100)), 
        col = rgb(0, 0, 0, 0.1), border = NA)



