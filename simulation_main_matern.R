################################################################################
# DGP Simulation Study: Compare Methods                                        #
# Cheng-Han Yu                                                                 #
################################################################################
rm(list = ls())


## load data sets
load("./Analysis/Data/SimCompareDataNew.Rdata", verbose = TRUE)
## change to load("./data/sim_data.rda", verbose = TRUE)

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

# # ===========
# ## GPR no der
# # ===========
# EB_gp_matern_3_2 <- foreach(k = 1:no_data) %dopar% {
#     res <- Rsolnp::solnp(pars = c(1, 1, 1), fun = log_mar_lik_gp_matern,
#                          LB = c(0.0001, 0.0001, 0.0001), 
#                          UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001), 
#                          control = list(TOL = 1e-5, trace = 0),
#                          y = YY[[k]]$y, H0 = H0_diff[[k]], nu = 1.5)
#     res$par
# }
# 
# EB_gp_matern_5_2 <- foreach(k = 1:no_data) %dopar% {
#     res <- Rsolnp::solnp(pars = c(1, 1, 1), fun = log_mar_lik_gp_matern,
#                          LB = c(0.0001, 0.0001, 0.0001), 
#                          UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001), 
#                          control = list(TOL = 1e-5, trace = 0),
#                          y = YY[[k]]$y, H0 = H0_diff[[k]], nu = 2.5)
#     res$par
# }
# 
# # ===========
# ## DGP-oracle
# # ===========
# EB_dgp_oracle_3_2 <- foreach(k = 1:no_data) %dopar% {
#     res <- Rsolnp::solnp(pars = c(1, 1, 1), 
#                          fun = log_mar_lik_gp_der_oracle_matern,
#                          LB = c(0.0001, 0.0001, 0.0001),
#                          UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
#                          control = list(TOL = 1e-5, trace = 0),
#                          y = YY[[k]]$y,
#                          x_vec = YY[[k]]$x, der_vec = cri_pts,
#                          H0 = H0_diff[[k]], nu = 1.5, 
#                          is.sig.par = FALSE)
#     res$par}
# 
# EB_dgp_oracle_5_2 <- foreach(k = 1:no_data) %dopar% {
#     res <- Rsolnp::solnp(pars = c(1, 1, 1), 
#                          fun = log_mar_lik_gp_der_oracle_matern,
#                          LB = c(0.0001, 0.0001, 0.0001),
#                          UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
#                          control = list(TOL = 1e-5, trace = 0),
#                          y = YY[[k]]$y,
#                          x_vec = YY[[k]]$x, der_vec = cri_pts,
#                          H0 = H0_diff[[k]], nu = 2.5, 
#                          is.sig.par = FALSE)
#     res$par}
# 
# # ===========
# ## DGP-single
# # ===========
# no_data <- 100
# system.time(StoEM_single_matern_3_2 <- foreach(k = 1:no_data) %dopar% {
#     stochastic_em_dgp_matern(y = YY[[k]]$y, x = YY[[k]]$x, H0 = H0_diff[[k]], 
#                       theta_init = c(1, 1, 1), epsilon = 1e-4, 
#                       D = 150, a = 0, b = 2, n_sample = 2000, max_iter = 100,
#                       lower = c(0.0001, 0.0001, 0.0001),
#                       upper = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
#                       shape1 = 1, shape2 = 1, nu = 1.5,
#                       ctrl = list(TOL = 1e-5, trace = 0),
#                       is.sig.par = FALSE,
#                       is.h.par = FALSE)})
# 
# system.time(StoEM_single_matern_5_2 <- foreach(k = 1:no_data) %dopar% {
#     stochastic_em_dgp_matern(y = YY[[k]]$y, x = YY[[k]]$x, H0 = H0_diff[[k]], 
#                              theta_init = c(1, 1, 1), epsilon = 1e-4, 
#                              D = 150, a = 0, b = 2, n_sample = 2000, max_iter = 100,
#                              lower = c(0.0001, 0.0001, 0.0001),
#                              upper = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
#                              shape1 = 1, shape2 = 1, nu = 2.5,
#                              ctrl = list(TOL = 1e-5, trace = 0),
#                              is.sig.par = FALSE,
#                              is.h.par = FALSE)})
# ===========
## DGP-multiple
# ===========
no_data <- 100

# stochastic_em_dgp_multi_matern <- function(y, x, H0, theta_init = c(1, 1, 1), epsilon = 1e-4, 
#                                            D = 100, a_vec = c(0, 1), b_vec = c(1, 2), 
#                                            n_sample = 4000, max_iter = 100,
#                                            lower = c(0.001, 0.001, 0.001),
#                                            upper = c(1/0.001, 1/0.001, 1/0.001),
#                                            ctrl = list(TOL = 1e-5, trace = 0),
#                                            ga_shape = 5, ga_rate = 5, nu = 1.5,
#                                            is.sig.par = TRUE, tune = 1)
    
system.time(StoEM_multi_matern_3_2 <- foreach(k = 1:no_data) %dopar% {
    stochastic_em_dgp_multi_matern(y = YY[[k]]$y, x = YY[[k]]$x, H0 = H0_diff[[k]], 
                                   theta_init = c(1, 1, 1), epsilon = 1e-4, 
                                   D = 150, a_vec = c(0, 1), b_vec = c(1, 2), 
                                   n_sample = 2000, max_iter = 100,
                                   lower = c(0.0001, 0.0001, 0.0001),
                                   upper = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
                                   nu = 1.5,
                                   ctrl = list(TOL = 1e-5, trace = 0),
                                   is.sig.par = FALSE)})

system.time(StoEM_multi_matern_5_2 <- foreach(k = 1:no_data) %dopar% {
    stochastic_em_dgp_multi_matern(y = YY[[k]]$y, x = YY[[k]]$x, H0 = H0_diff[[k]], 
                                   theta_init = c(1, 1, 1), epsilon = 1e-4, 
                                   D = 150, a_vec = c(0, 1), b_vec = c(1, 2), 
                                   n_sample = 2000, max_iter = 100,
                                   lower = c(0.0001, 0.0001, 0.0001),
                                   upper = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
                                   nu = 2.5,
                                   ctrl = list(TOL = 1e-5, trace = 0),
                                   is.sig.par = FALSE)})
stopCluster(cl)
# y = YY[[k]]$y
# x = YY[[k]]$x
# H0 = H0_diff[[k]]


# system.time(multi_dgp_matern_3_2_test <- stochastic_em_dgp_multi_matern(y = YY[[k]]$y, x = YY[[k]]$x, H0 = H0_diff[[k]], 
#                                                                         theta_init = c(1, 1, 1), epsilon = 1e-4, 
#                                                                         D = 150, a_vec = c(0, 1), b_vec = c(1, 2), 
#                                                                         n_sample = 2000, max_iter = 100,
#                                                                         lower = c(0.0001, 0.0001, 0.0001),
#                                                                         upper = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
#                                                                         nu = 1.5,
#                                                                         ctrl = list(TOL = 1e-5, trace = 0),
#                                                                         is.sig.par = FALSE))



#############################################
## Extract sample of t and EB estimates from StoEM_single and StoEM_multiple
#############################################
sample_t_EM_single_lst_matern_3_2 <- lapply(StoEM_single_matern_3_2, 
                                            function(x) { x$sample_t })
sample_t_EM_single_lst_matern_5_2 <- lapply(StoEM_single_matern_5_2, 
                                            function(x) { x$sample_t })
sample_t_EM_multi_lst_matern_3_2 <- lapply(StoEM_multi_matern_3_2, 
                                           function(x) { x$sample_t_mat })
sample_t_EM_multi_lst_matern_5_2 <- lapply(StoEM_multi_matern_5_2, 
                                           function(x) { x$sample_t_mat })

EB_EM_single_lst_matern_3_2 <- lapply(StoEM_single_matern_3_2, 
                                  function(x) { x$thetas[nrow(x$thetas), ] })
EB_EM_single_lst_matern_5_2 <- lapply(StoEM_single_matern_5_2, 
                                  function(x) { x$thetas[nrow(x$thetas), ] })
EB_EM_multi_lst_matern_3_2 <- lapply(StoEM_multi_matern_3_2, 
                                      function(x) { x$thetas[nrow(x$thetas), ] })
EB_EM_multi_lst_matern_5_2 <- lapply(StoEM_multi_matern_5_2, 
                                      function(x) { x$thetas[nrow(x$thetas), ] })


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
pred_gp_lst_matern_3_2 <- lapply(1:no_data, function(k) {
    get_pred_ci_gp_matern(eb_par = EB_gp_matern_3_2[[k]], x = YY[[k]]$x, x_test = x_test, 
                   y = YY[[k]]$y, nu = 1.5)
})


pred_gp_lst_matern_5_2 <- lapply(1:no_data, function(k) {
    get_pred_ci_gp_matern(eb_par = EB_gp_matern_5_2[[k]], x = YY[[k]]$x, x_test = x_test, 
                   y = YY[[k]]$y, nu = 2.5)
})
# # ===========
# ## DGP-oracle
# # ===========
pred_dgp_oracle_lst_3_2 <- lapply(1:no_data, function(k) {
    get_pred_ci_dgp_pt_matern(yJ = c(YY[[k]]$y, 0, 0),
                       x = YY[[k]]$x,
                       x_test = x_test,
                       idx_der = cri_pts,
                       sig = EB_dgp_oracle_3_2[[k]][1],
                       tau = EB_dgp_oracle_3_2[[k]][2],
                       l = EB_dgp_oracle_3_2[[k]][3],
                       nu = 1.5)
})


pred_dgp_oracle_lst_5_2 <- lapply(1:no_data, function(k) {
    get_pred_ci_dgp_pt_matern(yJ = c(YY[[k]]$y, 0, 0),
                       x = YY[[k]]$x,
                       x_test = x_test,
                       idx_der = cri_pts,
                       sig = EB_dgp_oracle_5_2[[k]][1],
                       tau = EB_dgp_oracle_5_2[[k]][2],
                       l = EB_dgp_oracle_5_2[[k]][3],
                       nu = 2.5)
})

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

system.time(pred_dgp_multi_lst_matern_3_2 <- foreach(k = 1:no_data) %dopar% {
    get_pi_t_matern(yJ = c(YY[[k]]$y, rep(0, 2)),
                    x = YY[[k]]$x, 
                    x_test = x_test, 
                    idx_der = sample_t_EM_multi_lst_matern_3_2[[k]],
                    sig = EB_EM_multi_lst_matern_3_2[[k]][1], 
                    tau = EB_EM_multi_lst_matern_3_2[[k]][2],
                    l = EB_EM_multi_lst_matern_3_2[[k]][3],
                    nu = 1.5)
})

system.time(pred_dgp_multi_lst_matern_5_2 <- foreach(k = 1:no_data) %dopar% {
    get_pi_t_matern(yJ = c(YY[[k]]$y, rep(0, 2)),
                    x = YY[[k]]$x, 
                    x_test = x_test, 
                    idx_der = sample_t_EM_multi_lst_matern_5_2[[k]],
                    sig = EB_EM_multi_lst_matern_5_2[[k]][1], 
                    tau = EB_EM_multi_lst_matern_5_2[[k]][2],
                    l = EB_EM_multi_lst_matern_5_2[[k]][3],
                    nu = 2.5)
})


# ===========
## DGP-single
# ===========
# 
# system.time(pred_dgp_single_lst_matern_3_2 <- foreach(k = 1:no_data) %dopar% {
#     get_pi_t_matern(yJ = c(YY[[k]]$y, rep(0, 1)),
#              x = YY[[k]]$x, 
#              x_test = x_test, 
#              idx_der = sample_t_EM_single_lst_matern_3_2[[k]],
#              sig = EB_EM_single_lst_matern_3_2[[k]][1], 
#              tau = EB_EM_single_lst_matern_3_2[[k]][2],
#              l = EB_EM_single_lst_matern_3_2[[k]][3],
#              nu = 1.5)
# })
# 
# system.time(pred_dgp_single_lst_matern_5_2 <- foreach(k = 1:no_data) %dopar% {
#     get_pi_t_matern(yJ = c(YY[[k]]$y, rep(0, 1)),
#              x = YY[[k]]$x, 
#              x_test = x_test, 
#              idx_der = sample_t_EM_single_lst_matern_5_2[[k]],
#              sig = EB_EM_single_lst_matern_5_2[[k]][1], 
#              tau = EB_EM_single_lst_matern_5_2[[k]][2],
#              l = EB_EM_single_lst_matern_5_2[[k]][3],
#              nu = 2.5)
# })
stopCluster(cl)

save(EB_gp_matern_3_2, EB_gp_matern_5_2, 
     EB_dgp_oracle_3_2, EB_dgp_oracle_5_2, 
     StoEM_single_matern_3_2, StoEM_single_matern_5_2, 
     StoEM_multi_matern_3_2, StoEM_multi_matern_5_2, 
     pred_gp_lst_matern_3_2, pred_gp_lst_matern_5_2, 
     pred_dgp_oracle_lst_3_2, pred_dgp_oracle_lst_5_2,
     pred_dgp_single_lst_matern_3_2, pred_dgp_single_lst_matern_5_2, 
     pred_dgp_multi_lst_matern_3_2, pred_dgp_multi_lst_matern_5_2,
     file = "./simulation_matern.RData")



load("./simulation_matern.RData", verbose = TRUE)
# =============================================================================
## Root Mean Square Error
# =============================================================================
true_fcn_val <- regfcn(x_test)
pred_mean_gp_matern_3_2 <- lapply(pred_gp_lst_matern_3_2, 
                                  function(x) x$mu_test)
pred_mean_gp_matern_5_2 <- lapply(pred_gp_lst_matern_5_2, 
                                  function(x) x$mu_test)
pred_mean_dgp_oracle_matern_3_2 <- lapply(pred_dgp_oracle_lst_3_2, 
                                   function(x) x$mu_test)
pred_mean_dgp_oracle_matern_5_2 <- lapply(pred_dgp_oracle_lst_5_2, 
                                   function(x) x$mu_test)
# pred_mean_dgp_multi <- lapply(pred_dgp_multi_lst, function(x) x$mu_test)
pred_mean_dgp_single_matern_3_2 <- lapply(pred_dgp_single_lst_matern_3_2, 
                                          function(x) x$mu_test)
pred_mean_dgp_single_matern_5_2 <- lapply(pred_dgp_single_lst_matern_5_2, 
                                          function(x) x$mu_test)

pred_mean_dgp_multi_matern_3_2 <- lapply(pred_dgp_multi_lst_matern_3_2, 
                                          function(x) x$mu_test)
pred_mean_dgp_multi_matern_5_2 <- lapply(pred_dgp_multi_lst_matern_5_2, 
                                          function(x) x$mu_test)




rmse_gp_matern_3_2 <- sapply(pred_mean_gp_matern_3_2, 
                             function(x) {rmse_f(x, true_fcn_val)})
mean(rmse_gp_matern_3_2)

rmse_gp_matern_5_2 <- sapply(pred_mean_gp_matern_5_2, 
                             function(x) {rmse_f(x, true_fcn_val)})
mean(rmse_gp_matern_5_2)


rmse_dgp_oracle_3_2 <- sapply(pred_mean_dgp_oracle_3_2, function(x){
    rmse_f(x, true_fcn_val)})
mean(rmse_dgp_oracle_3_2)

rmse_dgp_oracle_5_2 <- sapply(pred_mean_dgp_oracle_5_2, function(x){
    rmse_f(x, true_fcn_val)})
mean(rmse_dgp_oracle_5_2)

# 
# rmse_dgp_multi <- sapply(pred_mean_dgp_multi, function(x){
#     rmse_f(x, true_fcn_val)})
# mean(rmse_dgp_multi)

rmse_dgp_single_matern_3_2 <- sapply(pred_mean_dgp_single_matern_3_2, function(x){
    rmse_f(x, true_fcn_val)})
mean(rmse_dgp_single_matern_3_2)

rmse_dgp_single_matern_5_2 <- sapply(pred_mean_dgp_single_matern_5_2, function(x){
    rmse_f(x, true_fcn_val)})
mean(rmse_dgp_single_matern_5_2)



# =============================================================================
## Root Mean Square Error a function of x and Plotting 
# =============================================================================
## GPR
diff_sq_gp_matern_3_2 <- sapply(pred_mean_gp_matern_3_2, function(x) (x - true_fcn_val) ^ 2)
rmse_gp_x_matern_3_2 <- sqrt(apply(diff_sq_gp_matern_3_2, 1, mean))

diff_sq_gp_matern_5_2 <- sapply(pred_mean_gp_matern_5_2, function(x) (x - true_fcn_val) ^ 2)
rmse_gp_x_matern_5_2 <- sqrt(apply(diff_sq_gp_matern_5_2, 1, mean))

## Oracle-DGP
diff_sq_oracle_matern_3_2 <- sapply(pred_mean_dgp_oracle_matern_3_2, function(x) (x - true_fcn_val) ^ 2)
rmse_dgp_oracle_x_matern_3_2 <- sqrt(apply(diff_sq_oracle_matern_3_2, 1, mean))

diff_sq_oracle_matern_5_2 <- sapply(pred_mean_dgp_oracle_matern_5_2, function(x) (x - true_fcn_val) ^ 2)
rmse_dgp_oracle_x_matern_5_2 <- sqrt(apply(diff_sq_oracle_matern_5_2, 1, mean))



## Multiple-DGP
diff_sq_multi_matern_3_2 <- sapply(pred_mean_dgp_multi_matern_3_2, function(x) (x - true_fcn_val) ^ 2)
rmse_dgp_multi_x_matern_3_2 <- sqrt(apply(diff_sq_multi_matern_3_2, 1, mean))
diff_sq_multi_matern_5_2 <- sapply(pred_mean_dgp_multi_matern_5_2, function(x) (x - true_fcn_val) ^ 2)
rmse_dgp_multi_x_matern_5_2 <- sqrt(apply(diff_sq_multi_matern_5_2, 1, mean))


## Single-DGP
diff_sq_single_matern_3_2 <- sapply(pred_mean_dgp_single_matern_3_2, function(x) (x - true_fcn_val) ^ 2)
rmse_dgp_single_x_matern_3_2 <- sqrt(apply(diff_sq_single_matern_3_2, 1, mean))
diff_sq_single_matern_5_2 <- sapply(pred_mean_dgp_single_matern_5_2, function(x) (x - true_fcn_val) ^ 2)
rmse_dgp_single_x_matern_5_2 <- sqrt(apply(diff_sq_single_matern_5_2, 1, mean))

## plotting
## rmse_compare.png (need simulation_compare_song.r)
par(mfrow = c(1, 1))
par(mar = c(4, 4, 2, 1))
plot(x_test, rmse_gp_x_matern_3_2, xlab = "x", ylab = "RMSE", main = "RMSE comparison Matern_3_2", 
     col = 1, lty = 1, lwd = 3, type = "l", ylim = c(0, 0.25), axes = F)
axis(1)
axis(2, las = 2)
lines(x_test, rmse_dgp_oracle_x_matern_3_2, col = "purple", lty = 2, lwd = 3)
lines(x_test, rmse_dgp_multi_x_matern_3_2, col = "red", lty = 3, lwd = 3)
lines(x_test, rmse_dgp_single_x_matern_3_2, col = "blue", lty = 4, lwd = 3)
# lines(YY[[2]]$x, rmse_song_x, col = "green3", lty = 5, lwd = 3)
legend("bottomright", 
       legend = c("GPR Matern_3_2", "Oracle-DGP Matern_3_2", "Multiple-DGP Matern_3_2", "Single-DGP Matern_3_2"),
       col = c("black", "purple", "red", "blue"),
       lty = c(1, 2, 3, 4), 
       lwd = c(2, 2, 2, 2) + 1,
       bty = "n")

par(mfrow = c(1, 1))
par(mar = c(4, 4, 2, 1))
plot(x_test, rmse_gp_x_matern_5_2, xlab = "x", ylab = "RMSE", main = "RMSE comparison Matern_5_2", 
     col = 1, lty = 1, lwd = 3, type = "l", ylim = c(0, 0.25), axes = F)
axis(1)
axis(2, las = 2)
lines(x_test, rmse_dgp_oracle_x_matern_5_2, col = "purple", lty = 2, lwd = 3)
lines(x_test, rmse_dgp_multi_x_matern_5_2, col = "red", lty = 3, lwd = 3)
lines(x_test, rmse_dgp_single_x_matern_5_2, col = "blue", lty = 4, lwd = 3)
# lines(YY[[2]]$x, rmse_song_x, col = "green3", lty = 5, lwd = 3)
legend("bottomright", 
       legend = c("GPR Matern_3_2", "Oracle-DGP Matern_5_2", "Multiple-DGP Matern_5_2", "Single-DGP Matern_5_2"),
       col = c("black", "purple", "red", "blue"),
       lty = c(1, 2, 3, 4), 
       lwd = c(2, 2, 2, 2) + 1,
       bty = "n")

par(mfrow = c(2, 2))
par(mar = c(4, 4, 2, 1))
plot(x_test, rmse_gp_x_matern_3_2, xlab = "x", ylab = "RMSE", main = "GPR RMSE Matern", 
     col = 1, lty = 1, lwd = 3, type = "l", ylim = c(0, 0.25), axes = F)
axis(1)
axis(2, las = 2)
lines(x_test, rmse_gp_x_matern_5_2, col = "purple", lty = 2, lwd = 3)
# lines(x_test, rmse_dgp_multi_x_matern_5_2, col = "red", lty = 3, lwd = 3)
# lines(x_test, rmse_dgp_single_x_matern_5_2, col = "blue", lty = 4, lwd = 3)
# lines(YY[[2]]$x, rmse_song_x, col = "green3", lty = 5, lwd = 3)
legend("bottomright", 
       legend = c("GPR Matern_3_2", "GPR Matern_5_2"),
       col = c("black", "purple"),
       lty = c(1, 2), 
       lwd = c(2, 2) + 1,
       bty = "n")

# par(mfrow = c(1, 1))
par(mar = c(4, 4, 2, 1))
plot(x_test, rmse_dgp_oracle_x_matern_3_2, xlab = "x", ylab = "RMSE", main = "Oracle-DGP RMSE Matern", 
     col = 1, lty = 1, lwd = 3, type = "l", ylim = c(0, 0.25), axes = F)
axis(1)
axis(2, las = 2)
lines(x_test,rmse_dgp_oracle_x_matern_5_2, col = "purple", lty = 2, lwd = 3)
legend("bottomright", 
       legend = c("Oracle-DGP Matern_3_2", "Oracle-DGP Matern_5_2"),
       col = c("black", "purple"),
       lty = c(1, 2), 
       lwd = c(2, 2) + 1,
       bty = "n")

# par(mfrow = c(1, 1))
par(mar = c(4, 4, 2, 1))
plot(x_test, rmse_dgp_multi_x_matern_3_2, xlab = "x", ylab = "RMSE", main = "Multiple-DGP RMSE Matern", 
     col = 1, lty = 1, lwd = 3, type = "l", ylim = c(0, 0.25), axes = F)
axis(1)
axis(2, las = 2)
lines(x_test, rmse_dgp_multi_x_matern_5_2, col = "purple", lty = 2, lwd = 3)
# lines(x_test, rmse_dgp_multi_x_matern_5_2, col = "red", lty = 3, lwd = 3)
# lines(x_test, rmse_dgp_single_x_matern_5_2, col = "blue", lty = 4, lwd = 3)
# lines(YY[[2]]$x, rmse_song_x, col = "green3", lty = 5, lwd = 3)
legend("bottomright", 
       legend = c("Multiple-DGP Matern_3_2", "Multiple-DGP Matern_5_2"),
       col = c("black", "purple"),
       lty = c(1, 2), 
       lwd = c(2, 2) + 1,
       bty = "n")

# par(mfrow = c(1, 1))
par(mar = c(4, 4, 2, 1))
plot(x_test, rmse_dgp_single_x_matern_3_2, xlab = "x", ylab = "RMSE", main = "Single-DGP RMSE Matern", 
     col = 1, lty = 1, lwd = 3, type = "l", ylim = c(0, 0.25), axes = F)
axis(1)
axis(2, las = 2)
lines(x_test, rmse_dgp_single_x_matern_5_2, col = "purple", lty = 2, lwd = 3)
# lines(x_test, rmse_dgp_multi_x_matern_5_2, col = "red", lty = 3, lwd = 3)
# lines(x_test, rmse_dgp_single_x_matern_5_2, col = "blue", lty = 4, lwd = 3)
# lines(YY[[2]]$x, rmse_song_x, col = "green3", lty = 5, lwd = 3)
legend("bottomright", 
       legend = c("Single-DGP Matern_3_2", "Single-DGP Matern_5_2"),
       col = c("black", "purple"),
       lty = c(1, 2), 
       lwd = c(2, 2) + 1,
       bty = "n")


# par(mfrow = c(3, 4))
par(mfrow = c(3, 2))
# pred_hist_matern_3_2.png
for (k in 2) {
    # idx <- seq(min(YY[[k]]$x), max(YY[[k]]$x), length.out = 500)
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_gp_lst_matern_3_2[[k]]$mu_test,
                     CI_Low_f = pred_gp_lst_matern_3_2[[k]]$ci_low,
                     CI_High_f = pred_gp_lst_matern_3_2[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = .6, 
                     plot.type = "p", col.poly = rgb(0, 0, 0, 0.2), 
                     pred_lwd = 2, der_line_col = "black",
                     reg_fcn_col = "black", pred_col = "black",
                     # title = paste("No der: data", k), 
                     title = paste("GPR Matern_3_2"), 
                     # der_line_col = "black",
                     is.legend = FALSE)
    
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_oracle_lst_3_2[[k]]$mu_test,
                     CI_Low_f = pred_dgp_oracle_lst_3_2[[k]]$ci_low,
                     CI_High_f = pred_dgp_oracle_lst_3_2[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = .6, 
                     plot.type = "p", col.poly = rgb(0, 0, 0, 0.2), 
                     pred_lwd = 2, der_line_col = "black",
                     reg_fcn_col = "black", pred_col = "black",
                     # title = paste("(a): data", k), 
                     title = paste("oracle-DGP Matern_3_2"), 
                     # der_line_col = "black",
                     is.legend = FALSE, legend.loc = "bottomleft")
    
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_multi_lst_matern_3_2[[k]]$mu_test,
                     CI_Low_f = pred_dgp_multi_lst_matern_3_2[[k]]$ci_low,
                     CI_High_f = pred_dgp_multi_lst_matern_3_2[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = .6, 
                     plot.type = "p", col.poly = rgb(0,0,0,0.2), 
                     pred_lwd = 2, der_line_col = "black",
                     reg_fcn_col = "black", pred_col = "black",
                     # title = paste("(b): data", k),
                     title = paste("multiple-DGP Matern_3_2"), 
                     # der_line_col = "black",
                     is.legend = FALSE)
    
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_single_lst_matern_3_2[[k]]$mu_test,
                     CI_Low_f = pred_dgp_single_lst_matern_3_2[[k]]$ci_low,
                     CI_High_f = pred_dgp_single_lst_matern_3_2[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = .6, 
                     plot.type = "p", col.poly = rgb(0, 0, 0, 0.2), 
                     pred_lwd = 2, der_line_col = "black",
                     reg_fcn_col = "black", pred_col = "black",
                     # title = paste("(c): data", k), 
                     title = paste("single-DGP Matern_3_2"), 
                     # der_line_col = "black",
                     is.legend = FALSE)
    hist_t(sample_t_EM_multi_lst_matern_3_2[[k]], 0, 2, ylim = c(0, 1), col = "black", 
           den.line = FALSE, breaks = 20)
    title(list(paste("Distribution of t: multiple-DGP"), cex = 1.5))
    abline(v = cri_pts, col = "black", lty = 1, lwd = 1)
    hist_t(sample_t_EM_single_lst_matern_3_2[[k]], 0, 2, ylim = c(0, 1), col = "black", 
           den.line = FALSE, breaks = 20)
    title(list(paste("Distribution of t: single-DGP"), cex = 1.5))
    abline(v = cri_pts, col = "black", lty = 1, lwd = 1)
}

# par(mfrow = c(3, 4))
par(mfrow = c(3, 2))
# pred_hist_matern_5_2.png
for (k in 2) {
    # idx <- seq(min(YY[[k]]$x), max(YY[[k]]$x), length.out = 500)
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_gp_lst_matern_5_2[[k]]$mu_test,
                     CI_Low_f = pred_gp_lst_matern_5_2[[k]]$ci_low,
                     CI_High_f = pred_gp_lst_matern_5_2[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = .6, 
                     plot.type = "p", col.poly = rgb(0, 0, 0, 0.2), 
                     pred_lwd = 2, der_line_col = "black",
                     reg_fcn_col = "black", pred_col = "black",
                     # title = paste("No der: data", k), 
                     title = paste("GPR Matern_5_2"), 
                     # der_line_col = "black",
                     is.legend = FALSE)
    
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_oracle_lst_5_2[[k]]$mu_test,
                     CI_Low_f = pred_dgp_oracle_lst_5_2[[k]]$ci_low,
                     CI_High_f = pred_dgp_oracle_lst_5_2[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = .6, 
                     plot.type = "p", col.poly = rgb(0, 0, 0, 0.2), 
                     pred_lwd = 2, der_line_col = "black",
                     reg_fcn_col = "black", pred_col = "black",
                     # title = paste("(a): data", k), 
                     title = paste("oracle-DGP Matern_5_2"), 
                     # der_line_col = "black",
                     is.legend = FALSE, legend.loc = "bottomleft")
    
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_multi_lst_matern_5_2[[k]]$mu_test,
                     CI_Low_f = pred_dgp_multi_lst_matern_5_2[[k]]$ci_low,
                     CI_High_f = pred_dgp_multi_lst_matern_5_2[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = .6, 
                     plot.type = "p", col.poly = rgb(0,0,0,0.2), 
                     pred_lwd = 2, der_line_col = "black",
                     reg_fcn_col = "black", pred_col = "black",
                     # title = paste("(b): data", k),
                     title = paste("multiple-DGP Matern_5_2"), 
                     # der_line_col = "black",
                     is.legend = FALSE)
    
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_single_lst_matern_3_2[[k]]$mu_test,
                     CI_Low_f = pred_dgp_single_lst_matern_3_2[[k]]$ci_low,
                     CI_High_f = pred_dgp_single_lst_matern_3_2[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = .6, 
                     plot.type = "p", col.poly = rgb(0, 0, 0, 0.2), 
                     pred_lwd = 2, der_line_col = "black",
                     reg_fcn_col = "black", pred_col = "black",
                     # title = paste("(c): data", k), 
                     title = paste("single-DGP Matern_5_2"), 
                     # der_line_col = "black",
                     is.legend = FALSE)
    hist_t(sample_t_EM_multi_lst_matern_5_2[[k]], 0, 2, ylim = c(0, 3), col = "black", 
           den.line = FALSE, breaks = 20)
    title(list(paste("Distribution of t: multiple-DGP"), cex = 1.5))
    abline(v = cri_pts, col = "black", lty = 1, lwd = 1)
    hist_t(sample_t_EM_single_lst_matern_5_2[[k]], 0, 2, ylim = c(0, 3), col = "black", 
           den.line = FALSE, breaks = 20)
    title(list(paste("Distribution of t: single-DGP"), cex = 1.5))
    abline(v = cri_pts, col = "black", lty = 1, lwd = 1)
}
################################################################################
### Compare distribution of t
################################################################################
# ===============
# GPR no der
# ===============
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
# (1) RMSE
cl1_single_t_lst_matern_3_2 <- lapply(sample_t_EM_single_lst_matern_3_2, function(x){ x[x < 1] })
cl2_single_t_lst_matern_3_2 <- lapply(sample_t_EM_single_lst_matern_3_2, function(x){ x[x > 1] })

# cl1_rmse_within_single_lst <- lapply(cl1_single_t_lst, function(x) {
#     sqrt(mean((x - cri_pts[1]) ^ 2))
# })
# cl2_rmse_within_single_lst <- lapply(cl2_single_t_lst, function(x) {
#     sqrt(mean((x - cri_pts[2]) ^ 2))
# })

# cl1_sd_within_single_lst <- lapply(cl1_single_t_lst, function(x) sd(x))
# cl2_sd_within_single_lst <- lapply(cl2_single_t_lst, function(x) sd(x))

cl1_mean_single_vec_matern_3_2 <- unlist(lapply(cl1_single_t_lst_matern_3_2, mean))
cl2_mean_single_vec_matern_3_2 <- unlist(lapply(cl2_single_t_lst_matern_3_2, mean))

cl1_mode_single_vec_matern_3_2 <- unlist(lapply(cl1_single_t_lst_matern_3_2, function(mm) {
    density(mm)$x[which(density(mm)$y == max(density(mm)$y))]
}))

cl2_mode_single_vec_matern_3_2 <- unlist(lapply(cl2_single_t_lst_matern_3_2, function(mm) {
    density(mm)$x[which(density(mm)$y == max(density(mm)$y))]
}))


(cl1_rmse_between_single_matern_3_2 <- sqrt(mean((cl1_mean_single_vec_matern_3_2 - cri_pts[1]) ^ 2)))
(cl2_rmse_between_single_matern_3_2 <- sqrt(mean((cl2_mean_single_vec_matern_3_2 - cri_pts[2]) ^ 2)))

(cl1_rmse_between_single_mode_matern_3_2 <- sqrt(mean((cl1_mode_single_vec_matern_3_2 - cri_pts[1]) ^ 2)))
(cl2_rmse_between_single_mode_matern_3_2 <- sqrt(mean((cl2_mode_single_vec_matern_3_2 - cri_pts[2]) ^ 2)))
# 
# cl1_sd_between_single <- sd(cl1_mean_single_vec)
# cl2_sd_between_single <- sd(cl2_mean_single_vec)
# 
# cl1_overall_mean_single <- mean(cl1_mean_single_vec)
# cl2_overall_mean_single <- mean(cl2_mean_single_vec)
# -------------------------------------------------------
cl1_single_t_lst_matern_5_2 <- lapply(sample_t_EM_single_lst_matern_5_2, function(x){ x[x < 1] })
cl2_single_t_lst_matern_5_2 <- lapply(sample_t_EM_single_lst_matern_5_2, function(x){ x[x > 1] })


cl1_mean_single_vec_matern_5_2 <- unlist(lapply(cl1_single_t_lst_matern_5_2, mean))
cl2_mean_single_vec_matern_5_2 <- unlist(lapply(cl2_single_t_lst_matern_5_2, mean))

(cl1_rmse_between_single_matern_5_2 <- sqrt(mean((cl1_mean_single_vec_matern_5_2 - cri_pts[1]) ^ 2)))
(cl2_rmse_between_single_matern_5_2 <- sqrt(mean((cl2_mean_single_vec_matern_5_2 - cri_pts[2]) ^ 2)))


cl1_mode_single_vec_matern_5_2 <- unlist(lapply(cl1_single_t_lst_matern_5_2, function(mm) {
    density(mm)$x[which(density(mm)$y == max(density(mm)$y))]
}))

cl2_mode_single_vec_matern_5_2 <- unlist(lapply(cl2_single_t_lst_matern_5_2, function(mm) {
    density(mm)$x[which(density(mm)$y == max(density(mm)$y))]
}))

(cl1_rmse_between_single_mode_matern_5_2 <- sqrt(mean((cl1_mode_single_vec_matern_5_2 - cri_pts[1]) ^ 2)))
(cl2_rmse_between_single_mode_matern_5_2 <- sqrt(mean((cl2_mode_single_vec_matern_5_2 - cri_pts[2]) ^ 2)))


# ===============
# DGP-multiple
# ===============
# (1) RMSE
cl1_multi_t_lst_matern_3_2 <- lapply(sample_t_EM_multi_lst_matern_3_2, function(x){ x[x < 1] })
cl2_multi_t_lst_matern_3_2 <- lapply(sample_t_EM_multi_lst_matern_3_2, function(x){ x[x > 1] })

cl1_mean_multi_vec_matern_3_2 <- unlist(lapply(cl1_multi_t_lst_matern_3_2, mean))
cl2_mean_multi_vec_matern_3_2 <- unlist(lapply(cl2_multi_t_lst_matern_3_2, mean))

cl1_mode_multi_vec_matern_3_2 <- unlist(lapply(cl1_multi_t_lst_matern_3_2, function(mm) {
    density(mm)$x[which(density(mm)$y == max(density(mm)$y))]
}))

cl2_mode_multi_vec_matern_3_2 <- unlist(lapply(cl2_multi_t_lst_matern_3_2, function(mm) {
    density(mm)$x[which(density(mm)$y == max(density(mm)$y))]
}))


(cl1_rmse_between_multi_matern_3_2 <- sqrt(mean((cl1_mean_multi_vec_matern_3_2 - cri_pts[1]) ^ 2)))
(cl2_rmse_between_multi_matern_3_2 <- sqrt(mean((cl2_mean_multi_vec_matern_3_2 - cri_pts[2]) ^ 2)))

(cl1_rmse_between_multi_mode_matern_3_2 <- sqrt(mean((cl1_mode_multi_vec_matern_3_2 - cri_pts[1]) ^ 2)))
(cl2_rmse_between_multi_mode_matern_3_2 <- sqrt(mean((cl2_mode_multi_vec_matern_3_2 - cri_pts[2]) ^ 2)))

# -------------------------------------------------------
cl1_multi_t_lst_matern_5_2 <- lapply(sample_t_EM_multi_lst_matern_5_2, function(x){ x[x < 1] })
cl2_multi_t_lst_matern_5_2 <- lapply(sample_t_EM_multi_lst_matern_5_2, function(x){ x[x > 1] })


cl1_mean_multi_vec_matern_5_2 <- unlist(lapply(cl1_multi_t_lst_matern_5_2, mean))
cl2_mean_multi_vec_matern_5_2 <- unlist(lapply(cl2_multi_t_lst_matern_5_2, mean))

(cl1_rmse_between_multi_matern_5_2 <- sqrt(mean((cl1_mean_multi_vec_matern_5_2 - cri_pts[1]) ^ 2)))
(cl2_rmse_between_multi_matern_5_2 <- sqrt(mean((cl2_mean_multi_vec_matern_5_2 - cri_pts[2]) ^ 2)))


cl1_mode_multi_vec_matern_5_2 <- unlist(lapply(cl1_multi_t_lst_matern_5_2, function(mm) {
    density(mm)$x[which(density(mm)$y == max(density(mm)$y))]
}))

cl2_mode_multi_vec_matern_5_2 <- unlist(lapply(cl2_multi_t_lst_matern_5_2, function(mm) {
    density(mm)$x[which(density(mm)$y == max(density(mm)$y))]
}))

(cl1_rmse_between_multi_mode_matern_5_2 <- sqrt(mean((cl1_mode_multi_vec_matern_5_2 - cri_pts[1]) ^ 2)))
(cl2_rmse_between_multi_mode_matern_5_2 <- sqrt(mean((cl2_mode_multi_vec_matern_5_2 - cri_pts[2]) ^ 2)))








# ===============
# Plotting
# ===============
par(mfrow = c(2, 3))
for (k in 1:6) {
    hist_t(sample_t_EM_single_lst_matern_3_2[[k]], 0, 2, 
           ylim = c(0, 10), col = "blue", 
           den.line = FALSE)
    title(list(paste("Dist. of t of DGP-single-matern_3_2: Data", k), cex = 2))
    # hist_t(sample_t_EM_multi_lst_new[[k]], 0, 2, ylim = c(0, 10), col = "blue", 
    #        den.line = FALSE)
    # title(list(paste("Dist. of t of DGP-multiple: Data", k), cex = 2))
}

for (k in 1:6) {
    hist_t(sample_t_EM_single_lst_matern_5_2[[k]], 0, 2, 
           ylim = c(0, 10), col = "blue", 
           den.line = FALSE)
    title(list(paste("Dist. of t of DGP-single-matern_5_2: Data", k), cex = 2))
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
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
# pred_data1.png
k <- 2
# GPR no der 
plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                 mu_test = pred_gp_lst_matern_3_2[[k]]$mu_test,
                 CI_Low_f = pred_gp_lst_matern_3_2[[k]]$ci_low,
                 CI_High_f = pred_gp_lst_matern_3_2[[k]]$ci_high,
                 xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                 ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                 plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
                 pred_lwd = 2, title = paste("GPR_matern_3_2"), 
                 is.legend = FALSE, true_fcn = regfcn)

plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                 mu_test = pred_gp_lst_matern_5_2[[k]]$mu_test,
                 CI_Low_f = pred_gp_lst_matern_5_2[[k]]$ci_low,
                 CI_High_f = pred_gp_lst_matern_5_2[[k]]$ci_high,
                 xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                 ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                 plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
                 pred_lwd = 2, title = paste("GPR_matern_5_2"), is.legend = FALSE)
# DGP-oracle
plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                 mu_test = pred_dgp_oracle_lst_3_2[[k]]$mu_test,
                 CI_Low_f = pred_dgp_oracle_lst_3_2[[k]]$ci_low,
                 CI_High_f = pred_dgp_oracle_lst_3_2[[k]]$ci_high,
                 xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                 ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1,
                 plot.type = "p", col.poly = rgb(0, 0, 1, 0.1),
                 pred_lwd = 2, title = paste("oracle-DGP_matern_3_2"),
                 is.legend = FALSE, legend.loc = "bottomleft")

plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                 mu_test = pred_dgp_oracle_lst_5_2[[k]]$mu_test,
                 CI_Low_f = pred_dgp_oracle_lst_5_2[[k]]$ci_low,
                 CI_High_f = pred_dgp_oracle_lst_5_2[[k]]$ci_high,
                 xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                 ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1,
                 plot.type = "p", col.poly = rgb(0, 0, 1, 0.1),
                 pred_lwd = 2, title = paste("oracle-DGP_matern_5_2"),
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
plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                 mu_test = pred_dgp_single_lst_matern_3_2[[k]]$mu_test,
                 CI_Low_f = pred_dgp_single_lst_matern_3_2[[k]]$ci_low,
                 CI_High_f = pred_dgp_single_lst_matern_3_2[[k]]$ci_high,
                 xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                 ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                 plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
                 pred_lwd = 2, title = paste("single-DGP_matern_3_2"), 
                 is.legend = FALSE)
plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                 mu_test = pred_dgp_single_lst_matern_5_2[[k]]$mu_test,
                 CI_Low_f = pred_dgp_single_lst_matern_5_2[[k]]$ci_low,
                 CI_High_f = pred_dgp_single_lst_matern_5_2[[k]]$ci_high,
                 xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                 ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                 plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
                 pred_lwd = 2, title = paste("single-DGP_matern_5_2"), 
                 is.legend = FALSE)

# hist_t(sample_t_EM_multi_lst[[k]], 0, 2, ylim = c(0, 5), col = "blue", 
#        den.line = FALSE)
# title(list(paste("Distribution of t: multiple-DGP"), cex = 1.5))
# abline(v = cri_pts, col = "darkgreen", lty = 1, lwd = 1)
hist_t(sample_t_EM_single_lst_matern_3_2[[k]], 0, 2, 
       ylim = c(0, 2), col = "blue", 
       den.line = FALSE)
title(list(paste("Distribution of t: single-DGP_matern_3_2"), cex = 1.5))
abline(v = cri_pts, col = "darkgreen", lty = 1, lwd = 1)

hist_t(sample_t_EM_single_lst_matern_5_2[[k]], 0, 2, 
       ylim = c(0, 5), col = "blue", 
       den.line = FALSE)
title(list(paste("Distribution of t: single-DGP_matern_5_2"), cex = 1.5))
abline(v = cri_pts, col = "darkgreen", lty = 1, lwd = 1)

###### HPD interval ############################################################
library(HDInterval)
hdi(StoEM_one_lst[[100]]$sample_t, credMass = 0.95, allowSplit = TRUE)
xxx <- density(StoEM_one_lst[[1]]$sample_t, from = 0, to = 2, bw = 0.01)
plot(xxx)

CI_der_em_low <- matrix(0, 100, 2)
CI_der_em_high <- matrix(0, 100, 2)
for (k in 1:100) {
    hdi_res <- hdi(density(StoEM_one_lst[[k]]$sample_t, from = 0, to = 2, 
                           bw = 0.02), credMass = 0.95, allowSplit = TRUE)
    CI_der_em_low[k, ] <- hdi_res[1, ]
    CI_der_em_high[k, ] <- hdi_res[2, ]
}

apply(CI_der_em_low, 2, mean)
apply(CI_der_em_high, 2, mean)

CI_der_em <- hdi(density(StoEM_one_lst[[100]]$sample_t, 
                         from = 0, to = 2, bw = 0.01), 
                 credMass = 0.95, allowSplit = TRUE)


hdi(density(StoEM_one_lst[[100]]$sample_t, from = 0, to = 2, bw = 0.01), 
    credMass = 0.95, allowSplit = TRUE)[1, ]

##########################
hpd_interval_hist_lst_matern_3_2 <- lapply(StoEM_single_matern_3_2, 
                                           function(x) {
    get_hpd_interval_from_hist(samples = x$sample_t, breaks = 30, 
                               is.plot = FALSE)
})

hpd_interval_hist_lst_matern_5_2 <- lapply(StoEM_single_matern_5_2, 
                                           function(x) {
                                               get_hpd_interval_from_hist(samples = x$sample_t, breaks = 30, 
                                                                          is.plot = FALSE)
                                           })



par(mfcol = c(2, 3))
for (k in 1:3) {
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_single_lst_matern_3_2[[k]]$mu_test,
                     CI_Low_f = pred_dgp_single_lst_matern_3_2[[k]]$ci_low,
                     CI_High_f = pred_dgp_single_lst_matern_3_2[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.5), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 0, 1, 0.3), 
                     pred_lwd = 2, 
                     title = paste("Predictive f of DGP-single_matern_3_2: data", k), 
                     is.legend = TRUE, legend.loc = "bottomleft")
    hist_t(sample_t_EM_single_lst_matern_3_2[[k]], 0, 2, ylim = c(0, 4), col = "blue",
           den.line = FALSE, cri_pt = cri_pts)
    title(list(paste("Distribution of t: data", k), cex = 1.8))
    segments(x0 = hpd_interval_hist_lst_matern_3_2[[k]]$ci_lower[[1]], 
             y0 = hpd_interval_hist_lst_matern_3_2[[k]]$den_value, 
             x1 = hpd_interval_hist_lst_matern_3_2[[k]]$ci_upper[[1]], 
             y1 = hpd_interval_hist_lst_matern_3_2[[k]]$den_value, col = "red", lwd = 4)
    segments(x0 = hpd_interval_hist_lst_matern_3_2[[k]]$ci_lower[[2]], 
             y0 = hpd_interval_hist_lst_matern_3_2[[k]]$den_value, 
             x1 = hpd_interval_hist_lst_matern_3_2[[k]]$ci_upper[[2]], 
             y1 = hpd_interval_hist_lst_matern_3_2[[k]]$den_value, col = "red", lwd = 4)
} 

for (k in 1:3) {
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_single_lst_matern_5_2[[k]]$mu_test,
                     CI_Low_f = pred_dgp_single_lst_matern_5_2[[k]]$ci_low,
                     CI_High_f = pred_dgp_single_lst_matern_5_2[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.5), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 0, 1, 0.3), 
                     pred_lwd = 2, 
                     title = paste("Predictive f of DGP-single_matern_5_2: data", k), 
                     is.legend = TRUE, legend.loc = "bottomleft")
    hist_t(sample_t_EM_single_lst_matern_5_2[[k]], 0, 2, ylim = c(0, 4), col = "blue",
           den.line = FALSE, cri_pt = cri_pts)
    title(list(paste("Distribution of t: data", k), cex = 1.8))
    segments(x0 = hpd_interval_hist_lst_matern_5_2[[k]]$ci_lower[[1]], 
             y0 = hpd_interval_hist_lst_matern_5_2[[k]]$den_value, 
             x1 = hpd_interval_hist_lst_matern_5_2[[k]]$ci_upper[[1]], 
             y1 = hpd_interval_hist_lst_matern_5_2[[k]]$den_value, col = "red", lwd = 4)
    segments(x0 = hpd_interval_hist_lst_matern_5_2[[k]]$ci_lower[[2]], 
             y0 = hpd_interval_hist_lst_matern_5_2[[k]]$den_value, 
             x1 = hpd_interval_hist_lst_matern_5_2[[k]]$ci_upper[[2]], 
             y1 = hpd_interval_hist_lst_matern_5_2[[k]]$den_value, col = "red", lwd = 4)
} 

## pred_dist_sim1.png
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
# 
# for (i in seq(0.01, 0.1, length.out = 10)) {
#     plot(density(unlist(sample_t_EM_single_lst), from = 0, to = 2, bw = i))
#     print(hdi(density(unlist(sample_t_EM_single_lst), from = 0, to = 2, bw = i), 
#               credMass = 0.95, allowSplit = TRUE))
# }

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

################################################################################
# Comparision betweeen SE kernel and Matern kernel
# Need results from simulation_main.R
################################################################################
## RMSE

rmse_se_matern <- rbind(c(mean(rmse_gp), 
                          mean(rmse_dgp_oracle), 
                          mean(rmse_dgp_single)),
                        c(mean(rmse_gp_matern_5_2), 
                          mean(rmse_dgp_oracle_5_2), 
                          mean(rmse_dgp_single_matern_5_2)),
                        c(mean(rmse_gp_matern_3_2), 
                          mean(rmse_dgp_oracle_3_2), 
                          mean(rmse_dgp_single_matern_3_2)))

rownames(rmse_se_matern) <- c("Gaussian(SE)", "Matern_nu_5_2", "Matern_nu_3_2")
colnames(rmse_se_matern) <- c("GPR", "Oracle", "DGP-single")


rmse_se_matern
# ===============================================================
## mean predictive f
# ===============================================================
############
# GPR
############
# the true regression function is infinitely differentiable! 
# So SE kernel should be the best
par(mfrow = c(2, 3))
for (k in 1:6) {
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_gp_lst[[k]]$mu_test,
                     CI_Low_f = pred_gp_lst[[k]]$ci_low,
                     CI_High_f = pred_gp_lst[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
                     pred_lwd = 3, title = paste("GPR", k), is.legend = FALSE,
                     is.pred.f = FALSE, pred_lty = 1)
    lines(x_test, pred_gp_lst_matern_3_2[[k]]$mu_test, 
          col = "green3", lwd = 2.5, lty = 1)
    lines(x_test, pred_gp_lst_matern_5_2[[k]]$mu_test, 
          col = "orange", lwd = 2.5, lty = 1)
}
legend("bottomleft", c("true", "SE", "Matern_3_2", "Matern_5_2"), 
       lwd = rep(2.5, 4), bty = "n",
       col = c("red", "blue", "green3", "orange"))

par(mfrow = c(2, 3))
for (k in 1:6) {
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_gp_lst[[k]]$mu_test,
                     CI_Low_f = pred_gp_lst[[k]]$ci_low,
                     CI_High_f = pred_gp_lst[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 0, 1, 0.3), 
                     pred_lwd = 3, title = paste("GPR", k), is.legend = FALSE,
                     is.pred.f = TRUE, pred_lty = 1)
}

for (k in 1:6) {
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_gp_lst_matern_3_2[[k]]$mu_test,
                     CI_Low_f = pred_gp_lst_matern_3_2[[k]]$ci_low,
                     CI_High_f = pred_gp_lst_matern_3_2[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 1, 0, 0.3), 
                     pred_lwd = 3, title = paste("GPR_matern_3_2", k), 
                     is.legend = FALSE, pred_col = "green4",
                     is.pred.f = TRUE, pred_lty = 1)
}

for (k in 1:6) {
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_gp_lst_matern_5_2[[k]]$mu_test,
                     CI_Low_f = pred_gp_lst_matern_5_2[[k]]$ci_low,
                     CI_High_f = pred_gp_lst_matern_5_2[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(1, 0, 0, 0.3), 
                     pred_lwd = 3, title = paste("GPR_matern_5_2", k), 
                     is.legend = FALSE, pred_col = "orange",
                     is.pred.f = TRUE, pred_lty = 1)
}


############
# DGP-oracle
############
par(mfrow = c(2, 3))
for (k in 1:6) {
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_oracle_lst[[k]]$mu_test,
                     CI_Low_f = pred_dgp_oracle_lst[[k]]$ci_low,
                     CI_High_f = pred_dgp_oracle_lst[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
                     pred_lwd = 3, title = paste("Oracle-DGP", k), 
                     is.legend = FALSE, legend.loc = "bottomleft",
                     is.pred.f = FALSE, pred_lty = 1)
    
    lines(x_test, pred_dgp_oracle_lst_3_2[[k]]$mu_test, 
          col = "green3", lwd = 2.5, lty = 1)
    lines(x_test, pred_dgp_oracle_lst_5_2[[k]]$mu_test, 
          col = "orange", lwd = 2.5, lty = 1)
}

legend("bottomleft", c("true", "SE", "Matern_3_2", "Matern_5_2"), 
       lwd = rep(2.5, 4),
       col = c("red", "blue", "green3", "orange"))

par(mfrow = c(2, 3))
for (k in 1:6) {
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_oracle_lst[[k]]$mu_test,
                     CI_Low_f = pred_dgp_oracle_lst[[k]]$ci_low,
                     CI_High_f = pred_dgp_oracle_lst[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 0, 1, 0.3), 
                     pred_lwd = 3, title = paste("Oracle-DGP-SE", k), 
                     is.legend = FALSE,
                     is.pred.f = TRUE, pred_lty = 1)
}

for (k in 1:6) {
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_oracle_lst_3_2[[k]]$mu_test,
                     CI_Low_f = pred_dgp_oracle_lst_3_2[[k]]$ci_low,
                     CI_High_f = pred_dgp_oracle_lst_3_2[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 1, 0, 0.3), 
                     pred_lwd = 3, title = paste("Oracle-DGP_matern_3_2", k), 
                     is.legend = FALSE, pred_col = "green4",
                     is.pred.f = TRUE, pred_lty = 1)
}

for (k in 1:6) {
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_oracle_lst_5_2[[k]]$mu_test,
                     CI_Low_f = pred_dgp_oracle_lst_5_2[[k]]$ci_low,
                     CI_High_f = pred_dgp_oracle_lst_5_2[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(1, 0, 0, 0.3), 
                     pred_lwd = 3, title = paste("Oracle-DGP_matern_5_2", k), 
                     is.legend = FALSE, pred_col = "orange",
                     is.pred.f = TRUE, pred_lty = 1)
}

par(mfrow = c(1, 3))
k = 3
plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                 mu_test = pred_dgp_oracle_lst[[k]]$mu_test,
                 CI_Low_f = pred_dgp_oracle_lst[[k]]$ci_low,
                 CI_High_f = pred_dgp_oracle_lst[[k]]$ci_high,
                 xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                 ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                 plot.type = "p", col.poly = rgb(0, 0, 1, 0.3), 
                 pred_lwd = 3, title = paste("Oracle-DGP-SE dataset", k), 
                 is.legend = FALSE,
                 is.pred.f = TRUE, pred_lty = 1)
plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                 mu_test = pred_dgp_oracle_lst_5_2[[k]]$mu_test,
                 CI_Low_f = pred_dgp_oracle_lst_5_2[[k]]$ci_low,
                 CI_High_f = pred_dgp_oracle_lst_5_2[[k]]$ci_high,
                 xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                 ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                 plot.type = "p", col.poly = rgb(1, 0, 0, 0.3), 
                 pred_lwd = 3, title = paste("Oracle-DGP-matern_5_2 dataset", k), 
                 is.legend = FALSE, pred_col = "orange",
                 is.pred.f = TRUE, pred_lty = 1)
plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                 mu_test = pred_dgp_oracle_lst_3_2[[k]]$mu_test,
                 CI_Low_f = pred_dgp_oracle_lst_3_2[[k]]$ci_low,
                 CI_High_f = pred_dgp_oracle_lst_3_2[[k]]$ci_high,
                 xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                 ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                 plot.type = "p", col.poly = rgb(0, 1, 0, 0.3), 
                 pred_lwd = 3, title = paste("Oracle-DGP-matern_3_2 dataset", k), 
                 is.legend = FALSE, pred_col = "green4",
                 is.pred.f = TRUE, pred_lty = 1)


############
# DGP-single
############
par(mfrow = c(2, 3))
for (k in 1:6) {
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_single_lst[[k]]$mu_test,
                     CI_Low_f = pred_dgp_single_lst[[k]]$ci_low,
                     CI_High_f = pred_dgp_single_lst[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
                     pred_lwd = 3, title = paste("Single-DGP", k), 
                     is.legend = FALSE, legend.loc = "bottomleft",
                     is.pred.f = FALSE, pred_lty = 1)
    
    lines(x_test, pred_dgp_single_lst_matern_3_2[[k]]$mu_test, 
          col = "green3", lwd = 2.5, lty = 1)
    lines(x_test, pred_dgp_single_lst_matern_5_2[[k]]$mu_test, 
          col = "orange", lwd = 2.5, lty = 1)
}

legend("bottomleft", c("true", "SE", "Matern_3_2", "Matern_5_2"), 
       lwd = rep(2.5, 4),
       col = c("red", "blue", "green3", "orange"))

par(mfrow = c(2, 3))
for (k in 1:6) {
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_single_lst[[k]]$mu_test,
                     CI_Low_f = pred_dgp_single_lst[[k]]$ci_low,
                     CI_High_f = pred_dgp_single_lst[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 0, 1, 0.3), 
                     pred_lwd = 3, title = paste("Single-DGP-SE", k), 
                     is.legend = FALSE,
                     is.pred.f = TRUE, pred_lty = 1)
}

for (k in 1:6) {
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_single_lst_matern_3_2[[k]]$mu_test,
                     CI_Low_f = pred_dgp_single_lst_matern_3_2[[k]]$ci_low,
                     CI_High_f = pred_dgp_single_lst_matern_3_2[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 1, 0, 0.3), 
                     pred_lwd = 3, title = paste("Single-DGP_matern_3_2", k), 
                     is.legend = FALSE, pred_col = "green4",
                     is.pred.f = TRUE, pred_lty = 1)
}

for (k in 1:6) {
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_single_lst_matern_5_2[[k]]$mu_test,
                     CI_Low_f = pred_dgp_single_lst_matern_5_2[[k]]$ci_low,
                     CI_High_f = pred_dgp_single_lst_matern_5_2[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(1, 0, 0, 0.3), 
                     pred_lwd = 3, title = paste("Single-DGP_matern_5_2", k), 
                     is.legend = FALSE, pred_col = "orange",
                     is.pred.f = TRUE, pred_lty = 1)
}





















# ===============================================================
## 10 predictive f realizations
# ===============================================================
par(mfrow = c(2, 3))
for(k in 1:6) {
    matplot(x_test, t(pred_dgp_oracle_lst[[k]]$pred_f)[, 1:10], type = "l", lwd = 2,
            main = paste("SE 10 f realizations", k), xlab = "x", ylab = "y",
            ylim = c(0.3, 2.3))
    lines(x_test, pred_dgp_oracle_lst[[k]]$mu_test, 
          col = rgb(0, 0, 0, alpha = 0.6), lwd = 5)
    points(x = YY[[k]]$x, y = YY[[k]]$y, pch = 1, 
           col = rgb(0, 0, 0, alpha = 0.6), lwd = 2)
    abline(v = cri_pts, col = "purple")
}

for(k in 1:6) {
    matplot(x_test, t(pred_dgp_oracle_lst_3_2[[k]]$pred_f)[, 1:10], 
            type = "l", lwd = 2, ylim = c(0.3, 2.3),
            main = paste("Ma_3_2 10 f realizations", k), xlab = "x", ylab = "y")
    lines(x_test, pred_dgp_oracle_lst_3_2[[k]]$mu_test, 
          col = rgb(0, 0, 0, alpha = 0.6), lwd = 5)
    points(x = YY[[k]]$x, y = YY[[k]]$y, pch = 1, 
           col = rgb(0, 0, 0, alpha = 0.6), lwd = 2)
    abline(v = cri_pts, col = "purple")
}

for(k in 1:6) {
    matplot(x_test, t(pred_dgp_oracle_lst_5_2[[k]]$pred_f)[, 1:10], 
            type = "l", lwd = 2, ylim = c(0.3, 2.3),
            main = paste("Ma_5_2 10 f realizations", k), xlab = "x", ylab = "y")
    lines(x_test, pred_dgp_oracle_lst_5_2[[k]]$mu_test, 
          col = rgb(0, 0, 0, alpha = 0.6), lwd = 5)
    points(x = YY[[k]]$x, y = YY[[k]]$y, pch = 1, 
           col = rgb(0, 0, 0, alpha = 0.6), lwd = 2)
    abline(v = cri_pts, col = "purple")
}

par(mfrow = c(1, 3))
k = 3
matplot(x_test, t(pred_dgp_oracle_lst[[k]]$pred_f)[, 1:10], type = "l", lwd = 2,
        main = paste("SE 10 f realizations data", k), xlab = "x", ylab = "y",
        ylim = c(0.3, 2.3))
# lines(x_test, pred_dgp_oracle_lst[[k]]$mu_test, 
#       col = rgb(0, 0, 0, alpha = 0.6), lwd = 5)
points(x = YY[[k]]$x, y = YY[[k]]$y, pch = 1, 
       col = rgb(0, 0, 0, alpha = 0.6), lwd = 2)
abline(v = cri_pts, col = "purple")

matplot(x_test, t(pred_dgp_oracle_lst_5_2[[k]]$pred_f)[, 1:10], 
        type = "l", lwd = 2, ylim = c(0.3, 2.3),
        main = paste("Matern_5_2 10 f realizations data", k), xlab = "x", ylab = "y")
# lines(x_test, pred_dgp_oracle_lst_5_2[[k]]$mu_test, 
#       col = rgb(0, 0, 0, alpha = 0.6), lwd = 5)
points(x = YY[[k]]$x, y = YY[[k]]$y, pch = 1, 
       col = rgb(0, 0, 0, alpha = 0.6), lwd = 2)
abline(v = cri_pts, col = "purple")

matplot(x_test, t(pred_dgp_oracle_lst_3_2[[k]]$pred_f)[, 1:10], 
        type = "l", lwd = 2, ylim = c(0.3, 2.3),
        main = paste("Matern_3_2 10 f realizations data", k), xlab = "x", ylab = "y")
# lines(x_test, pred_dgp_oracle_lst_3_2[[k]]$mu_test, 
#       col = rgb(0, 0, 0, alpha = 0.6), lwd = 5)
points(x = YY[[k]]$x, y = YY[[k]]$y, pch = 1, 
       col = rgb(0, 0, 0, alpha = 0.6), lwd = 2)
abline(v = cri_pts, col = "purple")



par(mfrow = c(2, 3))
for(k in 1:6) {
    matplot(x_test, t(pred_gp_lst[[k]]$pred_f)[, 1:100], type = "l", lwd = 2,
            main = paste("SE 10 f realizations", k), xlab = "x", ylab = "y",
            ylim = c(0.3, 2.3))
    lines(x_test, pred_gp_lst[[k]]$mu_test, 
          col = rgb(0, 0, 0, alpha = 0.6), lwd = 5)
    points(x = YY[[k]]$x, y = YY[[k]]$y, pch = 1, 
           col = rgb(0, 0, 0, alpha = 0.6), lwd = 2)
}

for(k in 1:6) {
    matplot(x_test, t(pred_gp_lst_matern_3_2[[k]]$pred_f)[, 1:100], 
            type = "l", lwd = 2, ylim = c(0.3, 2.3),
            main = paste("Ma_3_2 10 f realizations", k), xlab = "x", ylab = "y")
    lines(x_test, pred_gp_lst_matern_3_2[[k]]$mu_test, 
          col = rgb(0, 0, 0, alpha = 0.6), lwd = 5)
    points(x = YY[[k]]$x, y = YY[[k]]$y, pch = 1, 
           col = rgb(0, 0, 0, alpha = 0.6), lwd = 2)
}

for(k in 1:6) {
    matplot(x_test, t(pred_gp_lst_matern_5_2[[k]]$pred_f)[, 1:100], 
            type = "l", lwd = 2, ylim = c(0.3, 2.3),
            main = paste("Ma_5_2 10 f realizations", k), xlab = "x", ylab = "y")
    lines(x_test, pred_gp_lst_matern_5_2[[k]]$mu_test, 
          col = rgb(0, 0, 0, alpha = 0.6), lwd = 5)
    points(x = YY[[k]]$x, y = YY[[k]]$y, pch = 1, 
           col = rgb(0, 0, 0, alpha = 0.6), lwd = 2)
}


par(mfrow = c(2, 3))
for(k in 1:6) {
    matplot(x_test, pred_dgp_single_lst[[k]]$pred_f[, 1:100], type = "l", lwd = 2,
            main = paste("SE 10 f realizations", k), xlab = "x", ylab = "y",
            ylim = c(0.3, 2.3))
    lines(x_test, pred_dgp_single_lst[[k]]$mu_test, 
          col = rgb(0, 0, 0, alpha = 0.6), lwd = 5)
    points(x = YY[[k]]$x, y = YY[[k]]$y, pch = 1, 
           col = rgb(0, 0, 0, alpha = 0.6), lwd = 2)
}

for(k in 1:6) {
    matplot(x_test, pred_dgp_single_lst_matern_3_2[[k]]$pred_f[, 1:100], 
            type = "l", lwd = 2, ylim = c(0.3, 2.3),
            main = paste("Ma_3_2 10 f realizations", k), xlab = "x", ylab = "y")
    lines(x_test, pred_dgp_single_lst_matern_3_2[[k]]$mu_test, 
          col = rgb(0, 0, 0, alpha = 0.6), lwd = 5)
    points(x = YY[[k]]$x, y = YY[[k]]$y, pch = 1, 
           col = rgb(0, 0, 0, alpha = 0.6), lwd = 2)
}

for(k in 1:6) {
    matplot(x_test, pred_dgp_single_lst_matern_5_2[[k]]$pred_f[, 1:100], 
            type = "l", lwd = 2, ylim = c(0.3, 2.3),
            main = paste("Ma_5_2 10 f realizations", k), xlab = "x", ylab = "y")
    lines(x_test, pred_dgp_single_lst_matern_5_2[[k]]$mu_test, 
          col = rgb(0, 0, 0, alpha = 0.6), lwd = 5)
    points(x = YY[[k]]$x, y = YY[[k]]$y, pch = 1, 
           col = rgb(0, 0, 0, alpha = 0.6), lwd = 2)
}




rmse_dgp_single_matern_3_2_each_f <- sapply(pred_dgp_single_lst_matern_3_2, 
                                            function(x) {
                                                apply(x$pred_f, 2, function(s) sqrt(mean((s - true_fcn_val) ^ 2)))
                                            })


mean(apply(rmse_dgp_single_matern_3_2_each_f, 2, mean))


rmse_dgp_single_matern_5_2_each_f <- sapply(pred_dgp_single_lst_matern_5_2, 
       function(x) {
           apply(x$pred_f, 2, function(s) sqrt(mean((s - true_fcn_val) ^ 2)))
        })


mean(apply(rmse_dgp_single_matern_5_2_each_f, 2, mean))


rmse_dgp_single_each_f <- sapply(pred_dgp_single_lst, function(x) {
                                                apply(x$pred_f, 2, function(s) sqrt(mean((s - true_fcn_val) ^ 2)))
                                            })


mean(apply(rmse_dgp_single_each_f, 2, mean))

rmse_table <- rbind(rmse_se_matern[, 3], c(mean(apply(rmse_dgp_single_each_f, 2, mean)),
                             mean(apply(rmse_dgp_single_matern_5_2_each_f, 2, mean)),
                             mean(apply(rmse_dgp_single_matern_3_2_each_f, 2, mean))))
rownames(rmse_table) <- c("RMSE_1", "RMSE_2")

library(xtable)

xtable(rmse_table, digits = 3, caption = "afw")
# ===============================================================
## Distribution of t
# ===============================================================
par(mfrow = c(2, 3))
for (k in 1:6) {
    hist_t(sample_t_EM_single_lst_matern_3_2[[k]], 0, 2, ylim = c(0, 4.5), 
           col = rgb(0, 1, 0, 1), 
           den.line = FALSE)
    title(list(paste("Dist. of t: single-DGP", k), cex = 1))
    abline(v = cri_pts, col = "purple", lty = 1, lwd = 2)
    Sys.sleep(1)
    hist(sample_t_EM_single_lst_matern_5_2[[k]], breaks = 30,
         freq = FALSE, add = TRUE, col = rgb(1, 0, 0, 0.4), border = "white")
    Sys.sleep(1)
    hist(sample_t_EM_single_lst[[k]], breaks = 30,
         freq = FALSE, add = TRUE, col = rgb(0, 0, 1, 0.4), border = "white")
    Sys.sleep(1)
}

legend(0.5, 4, c("Matern_3_2", "Matern_5_2", "SE"),
       col = c(rgb(0, 1, 0, 1), rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.6)),
       lwd = c(2,2,2), bty = "n", cex = 1)


par(mfrow = c(1, 3))
k = 3
hist_t(sample_t_EM_single_lst_matern_3_2[[k]], 0, 2, ylim = c(0, 6.5), 
       col = rgb(0, 1, 0, 1), breaks = 60,
       den.line = FALSE)
title(list(paste("Dist. of t of Matern_3_2 dataset", k), cex = 1))
abline(v = cri_pts, col = "purple", lty = 1, lwd = 2)
# Sys.sleep(1)
hist(sample_t_EM_single_lst_matern_5_2[[k]], breaks = 60, ylim = c(0, 6.5),
     freq = FALSE, add = FALSE, col = rgb(1, 0, 0, 0.8), border = "white",
     main = "", xlab = "t", xlim = c(0, 2))
abline(v = cri_pts, col = "purple", lty = 1, lwd = 2)
title(list(paste("Dist. of t of Matern_5_2 dataset", k), cex = 1))
# Sys.sleep(1)
hist(sample_t_EM_single_lst[[k]], breaks = 60, ylim = c(0, 6.5),
     freq = FALSE, add = FALSE, col = rgb(0, 0, 1, 0.8), border = "white",
     main = "", xlab = "t", xlim = c(0, 2))
abline(v = cri_pts, col = "purple", lty = 1, lwd = 2)
title(list(paste("Dist. of t of SE dataset", k), cex = 1))
# Sys.sleep(1)









