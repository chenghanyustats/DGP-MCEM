### supplement material
## Effect of sample size and SNR
pkg <- c("Matrix", "matrixcalc", "Rsolnp", "emulator", "R.utils")
lapply(pkg, require, character.only = TRUE)
sourceDirectory("./R")

# ===========================
## sig 0.25; n = 50 (normal amp)
# ===========================
gp_sig <- sapply(EB_gp, function(x) x[1])
mean(gp_sig)
sd(gp_sig)
dgp_single_sig <- sapply(sample_sig_EM_single_lst, function(x) mean(x))
mean(dgp_single_sig)
dgp_multi_sig <- sapply(sample_sig_EM_multi_lst, function(x) mean(x))
mean(dgp_multi_sig)


par(mfrow = c(5, 5))
for (i in 1:100) {
    hist(sample_sig_EM_single_lst[[i]], main = i, xlim = c(0.15, 0.6))
}


# ===========================
## sig 0.7; n = 50
# ===========================
load("./data/SimSig07Data.Rdata", verbose = TRUE)

H0_lst_sig07 <- lapply(YY_sig07, function(d) {
    outer(as.vector(d$x), as.vector(d$x), 
          FUN = function(x1, x2) (x1 - x2))
})


library(parallel)
library(doParallel)
detectCores()
cl <- makeCluster(6, type = "FORK")
registerDoParallel(cl)
getDoParWorkers()


# ===========
## GPR no der
# ===========
# EB_gp_lst_sig07 <- foreach(k = 1:no_data) %dopar% {
#     res <- Rsolnp::solnp(pars = c(1, 1, 1), fun = log_mar_lik_gp,
#                          LB = c(0.0001, 0.0001, 0.0001),
#                          UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
#                          control = list(TOL = 1e-5, trace = 0),
#                          y = YY_sig07[[k]]$y, H0 = H0_lst_sig07[[k]])
#     res$par
# }

# system.time(StoEM_single_lst_sig07 <- foreach(k = 1:100) %dopar% {
#     mcem_dgp(y = YY_sig07[[k]]$y, x = YY_sig07[[k]]$x, 
#              H0 = H0_lst_sig07[[k]],
#              theta_init = c(1, 1), epsilon = 1e-3, 
#              D = 500, a = 0, b = 2, n_sample = 5000, 
#              max_iter = 100,
#              lower = c(0.001, 0.001),
#              upper = c(1/0.001, 1/0.001), init_t = 0.5, 
#              init_sig = 1,
#              shape1 = 1, shape2 = 1,
#              ctrl = list(TOL = 1e-5, trace = 0),
#              ga_shape = 1/2, ga_rate = 1/2,
#              a_h = 1, b_h = 1, 
#              mixture_prob = c(0.5, 0.5),
#              is.h.par = FALSE)})

system.time(StoEM_multi_lst_sig07 <- foreach(k = 1:100) %dopar% {
    mcem_dgp_multi(y = YY_sig07[[k]]$y, x = YY_sig07[[k]]$x,
                   H0 = H0_lst_sig07[[k]],
                   theta_init = c(1, 1), epsilon = 1e-3,
                   D = 500, a_vec = c(0, 1), b_vec = c(1, 2),
                   n_sample = 5000,
                   max_iter = 100,
                   lower = c(0.001, 0.001),
                   upper = c(1/0.001, 1/0.001),
                   init_t = c(0.5, 1.5),
                   init_sig = 1,
                   shape1 = 1, shape2 = 1,
                   ctrl = list(TOL = 1e-5, trace = 0),
                   ga_shape = 1/2, ga_rate = 1/2,
                   a_h = 1, b_h = 1, is.h.par = FALSE)})

stopCluster(cl)
# 
# k <- 35
# 
# test <- mcem_dgp(y = YYn_lst$n10[[k]]$y, x = YYn_lst$n10[[k]]$x,
#          H0 = H0_lst_10[[k]],
#          theta_init = c(1, 1), epsilon = 1e-3,
#          D = 500, a = 0, b = 2, n_sample = 5000,
#          max_iter = 100,
#          lower = c(0.001, 0.001),
#          upper = c(1/0.001, 1/0.001), init_t = 0.5,
#          init_sig = .1,
#          shape1 = 1, shape2 = 1,
#          ctrl = list(TOL = 1e-5, trace = 0),
#          ga_shape = 1/2, ga_rate = 1/2,
#          a_h = 1, b_h = 1,
#          mixture_prob = c(0.5, 0.5),
#          is.h.par = FALSE)
# # 
# test <- mcem_dgp(y = YY_sig07[[k]]$y, x = YY_sig07[[k]]$x, 
#                  H0 = H0_lst_sig07[[k]],
#                  theta_init = c(1, 1), epsilon = 1e-3, 
#                  D = 500, a = 0, b = 2, n_sample = 5000, 
#                  max_iter = 100,
#                  lower = c(0.001, 0.001),
#                  upper = c(1/0.001, 1/0.001), init_t = 0.5, 
#                  init_sig = 1,
#                  shape1 = 1, shape2 = 1,
#                  ctrl = list(TOL = 1e-5, trace = 0),
#                  ga_shape = 1/2, ga_rate = 1/2,
#                  a_h = 1, b_h = 1, 
#                  mixture_prob = c(0.5, 0.5),
#                  is.h.par = FALSE)

find_mode_from_hist <- function(sample) {
    den <- density(sample)
    den$x[which.max(den$y)]
}

gp_sig07 <- sapply(EB_gp_lst_sig07, function(x) x[1])
mean(gp_sig07)
find_mode_from_hist(gp_sig07)
sd(gp_sig07)

sample_sig_EM_single_lst_sig07 <- lapply(StoEM_single_lst_sig07, function(x) { x$sample_sig })
dgp_single_sig07 <- sapply(sample_sig_EM_single_lst_sig07, function(x) mean(x))
mean(dgp_single_sig07)
sd(dgp_single_sig07)

dgp_single_sig07_mode <- sapply(sample_sig_EM_single_lst_sig07, function(x) find_mode_from_hist(x))
find_mode_from_hist(dgp_single_sig07_mode)
sd(dgp_single_sig07_mode)


sample_sig_EM_multi_lst_sig07 <- lapply(StoEM_multi_lst_sig07, function(x) { x$sample_sig })
dgp_multi_sig07 <- sapply(sample_sig_EM_multi_lst_sig07, function(x) mean(x))
mean(dgp_multi_sig07)
sd(dgp_multi_sig07)

dgp_multi_sig07_mode <- sapply(sample_sig_EM_multi_lst_sig07, function(x) find_mode_from_hist(x))
find_mode_from_hist(dgp_multi_sig07_mode)
sd(dgp_multi_sig07_mode)



# ===========================
## sig 0.25; n = 10
# ===========================
load("./data/SimDataSizeNew.Rdata", verbose = TRUE)
str(YYn_lst)
idx <- seq(0, 2, length.out = 500)
x_test <- sort(seq(0, 2, length.out = 100))

H0_lst_10 <- lapply(YYn_lst$n10, function(d) {
    outer(as.vector(d$x), as.vector(d$x), 
          FUN = function(x1, x2) (x1 - x2))
})



library(parallel)
library(doParallel)
detectCores()
cl <- makeCluster(6, type = "FORK")
registerDoParallel(cl)
getDoParWorkers()


# ===========
## GPR no der
# ===========
# EB_gp_lst_n10 <- foreach(k = 1:no_data) %dopar% {
#     res <- Rsolnp::solnp(pars = c(1, 1, 1), fun = log_mar_lik_gp,
#                          LB = c(0.0001, 0.0001, 0.0001),
#                          UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
#                          control = list(TOL = 1e-5, trace = 0),
#                          y = YYn_lst$n10[[k]]$y, H0 = H0_lst_10[[k]])
#     res$par
# }

# system.time(StoEM_single_lst_n10 <- foreach(k = 1:100) %dopar% {
#     mcem_dgp(y = YYn_lst$n10[[k]]$y, x = YYn_lst$n10[[k]]$x,
#              H0 = H0_lst_10[[k]],
#              theta_init = c(1, 1), epsilon = 1e-3,
#              D = 500, a = 0, b = 2, n_sample = 5000,
#              max_iter = 100,
#              lower = c(0.001, 0.001),
#              upper = c(1/0.001, 1/0.001), init_t = 0.5,
#              init_sig = 0.1,
#              shape1 = 1, shape2 = 1,
#              ctrl = list(TOL = 1e-5, trace = 0),
#              ga_shape = 1/2, ga_rate = 1/2,
#              a_h = 1, b_h = 1,
#              mixture_prob = c(0.5, 0.5),
#              is.h.par = FALSE)})

system.time(StoEM_multi_lst_n10 <- foreach(k = 1:100) %dopar% {
    mcem_dgp_multi(y = YYn_lst$n10[[k]]$y, x = YYn_lst$n10[[k]]$x,
                   H0 = H0_lst_10[[k]],
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
                   a_h = 1, b_h = 1, is.h.par = FALSE)})

stopCluster(cl)

gp_sig_n10 <- sapply(EB_gp_lst_n10, function(x) x[1])
mean(gp_sig_n10)
sd(gp_sig_n10)
find_mode_from_hist(gp_sig_n10)
sd(gp_sig_n10)



sample_sig_EM_single_lst_n10 <- lapply(StoEM_single_lst_n10, function(x) { x$sample_sig })
dgp_single_n10 <- sapply(sample_sig_EM_single_lst_n10, function(x) mean(x))
mean(dgp_single_n10)
sd(dgp_single_n10)

par(mfrow = c(5, 5)) 
for (i in 1:100) {
    hist(sample_sig_EM_single_lst_n10[[i]], main = i, xlim = c(0, 2))
}

dgp_single_n10_mode <- sapply(sample_sig_EM_single_lst_n10, function(x) find_mode_from_hist(x))
find_mode_from_hist(dgp_single_n10_mode)
sd(dgp_single_n10_mode)



sample_sig_EM_multi_lst_n10 <- lapply(StoEM_multi_lst_n10, function(x) { x$sample_sig })
dgp_multi_n10 <- sapply(sample_sig_EM_multi_lst_n10, function(x) mean(x))
mean(dgp_multi_n10)
sd(dgp_multi_n10)

dgp_multi_n10_mode <- sapply(sample_sig_EM_multi_lst_n10, function(x) find_mode_from_hist(x))
find_mode_from_hist(dgp_multi_n10_mode)
sd(dgp_multi_n10_mode)

# ===========================
## sig 0.25; n = 50 (small amp)
# ===========================
load("./Analysis/Data/SimCompareDataNew.Rdata", verbose = TRUE)


## small and large amp data have the same x
H0_diff_amp <- lapply(YY_small_amp, function(d) {
    outer(as.vector(d$x), as.vector(d$x),
          FUN = function(x1, x2) (x1 - x2))
})



library(parallel)
library(doParallel)
detectCores()
cl <- makeCluster(6, type = "FORK")
registerDoParallel(cl)
getDoParWorkers()


# ===========
## GPR no der
# ===========
EB_gp_lst_amp <- foreach(k = 1:no_data) %dopar% {
    res <- Rsolnp::solnp(pars = c(1, 1, 1), fun = log_mar_lik_gp,
                         LB = c(0.0001, 0.0001, 0.0001),
                         UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
                         control = list(TOL = 1e-5, trace = 0),
                         y = YY_small_amp[[k]]$y, H0 = H0_diff_amp[[k]])
    res$par
}

system.time(StoEM_single_lst_amp <- foreach(k = 1:100) %dopar% {
    mcem_dgp(y = YY_small_amp[[k]]$y, x = YY_small_amp[[k]]$x,
             H0 = H0_diff_amp[[k]],
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

system.time(StoEM_multi_lst_amp <- foreach(k = 1:100) %dopar% {
    mcem_dgp_multi(y = YY_small_amp[[k]]$y, x = YY_small_amp[[k]]$x,
                   H0 = H0_diff_amp[[k]],
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
                   a_h = 1, b_h = 1, is.h.par = FALSE)})

stopCluster(cl)

gp_sig_amp <- sapply(EB_gp_lst_amp, function(x) x[1])
mean(gp_sig_amp)
find_mode_from_hist(gp_sig_amp)
sd(gp_sig_amp)


sample_sig_EM_single_lst_amp <- lapply(StoEM_single_lst_amp, function(x) { x$sample_sig })
dgp_single_amp <- sapply(sample_sig_EM_single_lst_amp, function(x) mean(x))
mean(dgp_single_amp)
sd(dgp_single_amp)

dgp_single_amp_mode <- sapply(sample_sig_EM_single_lst_amp, function(x) find_mode_from_hist(x))
find_mode_from_hist(dgp_single_amp_mode)
sd(dgp_single_amp_mode)


sample_sig_EM_multi_lst_amp <- lapply(StoEM_multi_lst_amp, function(x) { x$sample_sig })
dgp_multi_amp <- sapply(sample_sig_EM_multi_lst_amp, function(x) mean(x))
mean(dgp_multi_amp)
sd(dgp_multi_amp)

dgp_multi_amp_mode <- sapply(sample_sig_EM_multi_lst_amp, function(x) find_mode_from_hist(x))
find_mode_from_hist(dgp_multi_amp_mode)
sd(dgp_multi_amp_mode)





