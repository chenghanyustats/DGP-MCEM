################################################################################
# real data analysis
# multiple MCEM
################################################################################
pkg <- c("Matrix", "matrixcalc", "Rsolnp", "emulator", "R.utils")
lapply(pkg, require, character.only = TRUE)
sourceDirectory("./R")


multi_mcem_lst_3_3_young <- list()
multi_mcem_lst_3_3_old <- list()
idx_11_good <- c(1, 5, 7, 9, 11, 12, 13, 14, 15, 17, 20)
subj_y_good <- subj_y[, idx_11_good]

system.time(for(k in 1:11) {
    print(paste("k=", k))
    multi_mcem_lst_3_3_young[[k]] <- mcem_dgp_multi(subj_y_good[, k], x = x_s, H0 = H0_mcmc,
                                              theta_init = c(1, 0.1), epsilon = 1e-4,
                                              D = 100, a_vec = c(min(x_s), median(x_s)),
                                              b_vec = c(median(x_s), max(x_s)),
                                              n_sample = 5000, max_iter = 100,
                                              lower = c(0.0001, 0.0001),
                                              upper = c(1/0.0001, 1/0.0001),
                                              init_t = c(0.4, 0.6), 
                                              init_sig = 1,
                                              ctrl = list(TOL = 1e-5, trace = 0),
                                              ga_shape = 200, ga_rate = 270,
                                              shape1 = c(3, 3), shape2 = c(3, 3))
})

multi_mcem_lst_3_3_young[[9]] <- mcem_dgp_multi(subj_y_good[, 9], x = x_s, H0 = H0_mcmc,
                                                theta_init = c(1, 0.1), epsilon = 1e-4,
                                                D = 100, a_vec = c(min(x_s), median(x_s)),
                                                b_vec = c(median(x_s), max(x_s)),
                                                n_sample = 5000, max_iter = 100,
                                                lower = c(0.0001, 0.0001),
                                                upper = c(1/0.0001, 1/0.0001),
                                                init_t = c(0.4, 0.6), 
                                                init_sig = 1,
                                                ctrl = list(TOL = 1e-5, trace = 0),
                                                ga_shape = 200, ga_rate = 270,
                                                shape1 = c(3, 3), shape2 = c(3, 3.5))

hist(multi_mcem_lst_3_3_young[[9]]$sample_t_mat, 70)
hist(multi_mcem_lst_3_3_young[[1]]$sample_sig, 70)



y = subj_y_good[, 9]
x = x_s
H0 = H0_mcmc
theta_init = c(1, 0.1)
epsilon = 1e-4
D = 100
a_vec = c(min(x_s), median(x_s))
b_vec = c(median(x_s), max(x_s))
n_sample = 5000
max_iter = 100
lower = c(0.0001, 0.0001)
upper = c(1/0.0001, 1/0.0001)
init_t = c(0.4, 0.6)
init_sig = 1
ctrl = list(TOL = 1e-5, trace = 0)
ga_shape = 200
ga_rate = 270
shape1 = c(3, 3)
shape2 = c(3, 3)




system.time(for(k in 1:11) {
    print(paste("k=", k))
    multi_mcem_lst_3_3_old[[k]] <- mcem_dgp_multi(subj_y_old[, k], x = x_s, H0 = H0_mcmc,
                                              theta_init = c(1, 0.1), epsilon = 1e-4,
                                              D = 100, a_vec = c(min(x_s), median(x_s)),
                                              b_vec = c(median(x_s), max(x_s)),
                                              n_sample = 5000, max_iter = 100,
                                              lower = c(0.0001, 0.0001),
                                              upper = c(1/0.0001, 1/0.0001),
                                              init_t = c(0.4, 0.6), 
                                              init_sig = 1,
                                              ctrl = list(TOL = 1e-5, trace = 0),
                                              ga_shape = 60, ga_rate = 140,
                                              shape1 = c(3, 3), shape2 = c(3, 3))
})


multi_mcem_lst_3_3_old[[9]] <- mcem_dgp_multi(subj_y_good[, 9], x = x_s, H0 = H0_mcmc,
                                 theta_init = c(1, 0.1), epsilon = 1e-4,
                                 D = 100, a_vec = c(min(x_s), median(x_s)),
                                 b_vec = c(median(x_s), max(x_s)),
                                 n_sample = 2000, max_iter = 100,
                                 lower = c(0.0001, 0.0001),
                                 upper = c(1/0.0001, 1/0.0001),
                                 init_t = c(0.4, 0.6),
                                 init_sig = 1,
                                 ctrl = list(TOL = 1e-5, trace = 0),
                                 ga_shape = 200, ga_rate = 270,
                                 mixture_prob = c(0.5, 0.5),
                                 shape1 = c(3, 3), shape2 = c(3, 5))

hist(multi_mcem_3_3$sample_t_mat, breaks = 170)
hist(multi_mcem_3_3$sample_sig, breaks = 170)
 
# 
mcem_3_3 <- mcem_dgp(subj_y_good[, 9], x = x_s, H0 = H0_mcmc,
                     theta_init = c(1, .1), epsilon = 1e-4,
                     D = 100, n_sample = 5000, a = min(x_s), b = max(x_s),
                     max_iter = 100,
                     lower = c(0.001, 0.001),
                     upper = c(1/0.001, 1/0.001), init_t = 0.5,
                     init_sig = 1,
                     shape1 = 3, shape2 = 3,
                     ctrl = list(TOL = 1e-5, trace = 0),
                     ga_shape = 200, ga_rate = 270,
                     a_h = 1, b_h = 1, mixture_prob = c(0.5, 0.5),
                     is.h.par = FALSE)

hist(mcem_3_3$sample_t, breaks = 70)
# hist(multi_mcem_3_3$sample_t_mat)
# hist(multi_mcem_3_3$sample_sig)
# 
# hist(multi_StoEM_lst_3_3[[1]]$sample_t_mat)


# -----------------------
#### young
# -----------------------

eb_theta_multi_mcem <- 
    multi_mcem_lst_3_3_young[[1]]$theta_mat[nrow(multi_mcem_lst_3_3_young[[1]]$theta_mat), ]

mcmc_output_uniform_mcem <- res_uniform_mcem$mcmc_output
draws_msec_uniform_mcem <- apply(mcmc_output_uniform_mcem$draws, 2, function(a) {
    t_sam <- a * (max(x_all) - min(x_all)) + min(x_all)
    t_sam_msec <- sec_to_msec(t_sam, 200)
    return(t_sam_msec)})



multi_t_mcem_lst_3_3 <- lapply(multi_mcem_lst_3_3_young, function(x) {
    x$sample_t_mat
})

multi_EB_mcem_lst_3_3 <- lapply(multi_mcem_lst_3_3_young, function(x) {
    x$thetas[nrow(x$thetas), ]
})



mu_lst_young <- vector("list", 11)
var_lst_young <- vector("list", 11)
proportion_lst_young <- vector("list", 11)

draws_msec <- lapply(multi_t_mcem_lst_3_3, function(a) {
    t_sam <- a * (max(x_all) - min(x_all)) + min(x_all)
    t_sam_msec <- sec_to_msec(t_sam, 200)
    return(t_sam_msec)
})


par(mfrow = c(4, 3))
for (s in 1:11) {
    # data <- draws_msec_uniform_mcem[, s]
    data <- as.vector(draws_msec[[s]])
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

order_prop <- lapply(proportion_lst_young, function(x) order(x, decreasing = TRUE))
t1_t2_mean_vec_young <- sapply(1:11, function(x) mu_lst_young[[x]][order_prop[[x]]][1:2])

order_mean <- apply(t1_t2_mean_vec_young, 2, order)


t1_mean_vec_young <- apply(t1_t2_mean_vec_young, 2, min)
t2_mean_vec_young <- apply(t1_t2_mean_vec_young, 2, max)

t_mean_mat_young <- cbind(t1_mean_vec_young, t2_mean_vec_young)

# min_idx_young <- sapply(mu_lst_young, function(x) which.min(x))

t1_t2_var_vec_young <- sapply(1:11, function(x) var_lst_young[[x]][order_prop[[x]]][1:2])

t1_t2_var_vec_young <- sapply(1:11, function(x) t1_t2_var_vec_young[, x][order_mean[, x]])
t1_var_vec_young <- t1_t2_var_vec_young[1, ]
t2_var_vec_young <- t1_t2_var_vec_young[2, ]

(r1_r2_mean_young <- apply(t_mean_mat_young, 2, mean))
(r1_r2_sd_young <- apply(t_mean_mat_young, 2, sd))
names(r1_r2_mean_young) <- names(r1_r2_sd_young) <- c("r1", "r2")

hpdi_young_lst <- vector("list", 11)
# hpdi_elderly_lst <- vector("list", 11)

for (k in 1:11) {
    hpdi_young <- get_hpd_interval_from_hist(samples = as.vector(draws_msec[[k]]),
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


# -----------------------
#### old
# -----------------------


# eb_theta_uniform_old <- res_uniform_old$theta_mat[nrow(res_uniform_old$theta_mat), ]
# 
# mcmc_output_uniform_old <- res_uniform_old$mcmc_output
# draws_msec_uniform_old <- apply(mcmc_output_uniform_old$draws, 2, function(a) {
#     t_sam <- a * (max(x_all) - min(x_all)) + min(x_all)
#     t_sam_msec <- sec_to_msec(t_sam, 200)
#     return(t_sam_msec)})





multi_t_mcem_lst_3_3_old <- lapply(multi_mcem_lst_3_3_old, function(x) {
    x$sample_t_mat
})

multi_EB_mcem_lst_3_3_old <- lapply(multi_mcem_lst_3_3_old, function(x) {
    x$thetas[nrow(x$thetas), ]
})

draws_msec_old <- lapply(multi_t_mcem_lst_3_3_old, function(a) {
    t_sam <- a * (max(x_all) - min(x_all)) + min(x_all)
    t_sam_msec <- sec_to_msec(t_sam, 200)
    return(t_sam_msec)
})



mu_lst_elderly <- vector("list", 11)
var_lst_elderly <- vector("list", 11)
proportion_lst_elderly <- vector("list", 11)
par(mfrow = c(4, 3))
for (s in 1:11) {
    # data <- draws_msec_uniform_old[, s]
    data <- as.vector(draws_msec_old[[s]])
    out_test <- mixmodCluster(data, nbCluster = 2, criterion = "ICL")
    mu_lst_elderly[[s]] <- out_test@results[[1]]@parameters@mean
    var_lst_elderly[[s]] <- out_test@results[[1]]@parameters@variance
    proportion_lst_elderly[[s]] <- out_test@results[[1]]@parameters@proportions
    hist_mix_gaussian(data = data, mu_vec = mu_lst_elderly[[s]], 
                      var_vec = var_lst_elderly[[s]], 
                      proportion_vec = proportion_lst_elderly[[s]], 
                      title = paste("Hist with GMM - old subj", s),
                      col = "blue")
}




hist_mix_gaussian <- function(data, mu_vec, var_vec, proportion_vec,
         title = "", col = "navy", breaks = 40, xlab = "Time (msec)", 
         is.curve = FALSE, ylim = c(0, 0.05)) {
    ### one dimensional histogram with Gaussian density fitting
    hist(data, freq = FALSE, xlim = c(50, 250), breaks = breaks, ylim = ylim,
         col = col, border = "white", main = "", xlab = xlab)
    title(main = title)
    if (is.curve) {
        for (i in 1:length(mu_vec)) {
            curve(proportion_vec[i] * dnorm(x, mean = mu_vec[i], 
                                            sd = sqrt(var_vec[[i]])), 
                  col = i + 1, lwd = 3, add = TRUE, yaxt = "n")
        }
    }
}
par(mfrow = c(1, 1))
s <- 5
data <- as.vector(draws_msec_old[[s]])
hist_mix_gaussian(data = data, mu_vec = mu_lst_elderly[[s]], 
                  var_vec = var_lst_elderly[[s]], 
                  proportion_vec = proportion_lst_elderly[[s]], 
                  title = "Multiple-DGP-Older t", xlab = "Time (msec)",
                  col = "gray", breaks = 70)

mcem_3_3_old <- mcem_dgp(subj_y_old[, s], x = x_s, H0 = H0_mcmc,
                         theta_init = c(1, .1), epsilon = 1e-4,
                         D = 100, n_sample = 5000, a = min(x_s), b = max(x_s),
                         max_iter = 100,
                         lower = c(0.001, 0.001),
                         upper = c(1/0.001, 1/0.001), init_t = 0.5,
                         init_sig = 1,
                         shape1 = 3, shape2 = 3,
                         ctrl = list(TOL = 1e-5, trace = 0),
                         ga_shape = 60, ga_rate = 140,
                         a_h = 1, b_h = 1, mixture_prob = c(0.5, 0.5),
                         is.h.par = FALSE)

t_sam <- mcem_3_3_old$sample_t * (max(x_all) - min(x_all)) + min(x_all)
t_sam_msec <- sec_to_msec(t_sam, 200)
hist(t_sam_msec, border = "white", xlim = c(50, 250), las = 1,
     main = "Single-DGP-Older t", xlab = "Time (msec)", breaks = 50)
hist(data, border = "white", xlim = c(50, 250), las = 1,
     main = "Multiple-DGP-Older t", xlab = "Time (msec)", breaks = 60)


hist(mcem_3_3_old$sample_sig, border = "white",las = 1,
     main = "Single-DGP-Older sigma", xlab = "sigma", breaks = 25)
hist(multi_mcem_lst_3_3_old[[s]]$sample_sig, border = "white",las = 1,
     main = "Multiple-DGP-Older sigma", xlab = "sigma", breaks = 25)




par(mfrow = c(1, 1))

# hist_mix_gaussian(data = data, mu_vec = mu_lst_young[[s]], 
#                   var_vec = var_lst_young[[s]],
#                   proportion_vec = proportion_lst_young[[s]], 
#                   title = "Multiple-DGP-Young t",
#                   col = "gray", breaks = 70)
# 
# hist(multi_mcem_lst_3_3_young[[s]]$sample_sig, freq = FALSE, border = "white",
#      main = "Multiple-DGP-Young sigma", xlab = "sigma", breaks = 30)

mcem_3_3_young <- mcem_dgp(subj_y_good[, s], x = x_s, H0 = H0_mcmc,
                         theta_init = c(1, .1), epsilon = 1e-4,
                         D = 100, n_sample = 5000, a = min(x_s), b = max(x_s),
                         max_iter = 100,
                         lower = c(0.001, 0.001),
                         upper = c(1/0.001, 1/0.001), init_t = 0.5,
                         init_sig = 1,
                         shape1 = 3, shape2 = 3,
                         ctrl = list(TOL = 1e-5, trace = 0),
                         ga_shape = 200, ga_rate = 270,
                         a_h = 1, b_h = 1, mixture_prob = c(0.5, 0.5),
                         is.h.par = FALSE)

t_sam <- mcem_3_3_young$sample_t * (max(x_all) - min(x_all)) + min(x_all)
t_sam_msec <- sec_to_msec(t_sam, 200)
hist(t_sam_msec, border = "white", xlim = c(50, 250), las = 1,
     main = "Single-DGP-Young t", xlab = "Time (msec)", breaks = 50)
s <- 8
data <- as.vector(draws_msec[[s]])
hist(data, border = "white", xlim = c(50, 250), las = 1,
     main = "Multiple-DGP-Young t", xlab = "Time (msec)", breaks = 60)


hist(mcem_3_3_young$sample_sig, border = "white", las = 1,
     main = "Single-DGP-Young sigma", xlab = "sigma", breaks = 25)
hist(multi_mcem_lst_3_3_young[[s]]$sample_sig, border = "white", las = 1,
     main = "Multiple-DGP-Young sigma", xlab = "sigma", breaks = 25)


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
    hpdi_elderly <- get_hpd_interval_from_hist(samples = as.vector(draws_msec_old[[k]]),
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



