### original get_pred_ci_t()
# get_pred_ci_t(yJ = c(y_voiced[, data_idx[k]], 0, 0),
#               x = x_s, 
#               x_test = x_test_s, 
#               idx_der = multi_t_EM_voiced_lst[[k]],
#               sig = multi_EB_EM_voiced_lst[[k]][1], 
#               tau = multi_EB_EM_voiced_lst[[k]][2],
#               h = multi_EB_EM_voiced_lst[[k]][3])


# yJ = c(y_voiced[, 11], 0, 0)
# 
# x = x_s
# x_test = x_test_s
# idx_der = multi_t_EM_voiced_lst[[1]]
# sig = multi_EB_EM_voiced_lst[[1]][1] 
# tau = multi_EB_EM_voiced_lst[[1]][2]
# h = multi_EB_EM_voiced_lst[[1]][3]

get_pred_ci_t <- function(yJ, x, x_test, idx_der, sig, tau, h) {
    # if (is.null(phi)) {
    #     phi <- 1 / h ^ 2
    # } 
    len_test <- length(x_test)
    mm <- dim(idx_der)[1]
    predictive_f_der <- matrix(0, mm, len_test)
    
    Kff <- compute_cov_1d(idx1 = x, idx2 = x, tau = tau, h = h)
    Kffnew <- compute_cov_1d(idx1 = x, idx2 = x_test, tau = tau, h = h)
    Kfnewfnew <- compute_cov_1d(idx1 = x_test, idx2 = x_test, tau = tau, h = h)
    
    for (i in 1:mm) {
        # if (i %% 100 == 0) print(i)
        dd <- idx_der[i, ]
        K_joint <- compute_joint_cov(idx_obs = x, idx_der = dd,
                                     tau = tau, h = h, sig = sig)
        Kdfnew <- computeCovDer1_h(idx1 = dd, idx2 = x_test, 
                                   tau = tau, h = h) 
        Kbar <- rbind(Kffnew, Kdfnew)
        # if(!is.positive.definite(K_joint)){
        #     K_joint <- as.matrix(nearPD(K_joint)$mat)
        # }
        # KJ_Kbar <- quad.form.inv(K_joint, Kbar)
        # Sig_test_der <- Kfnewfnew - KJ_Kbar
        Sig_test_der <- Kfnewfnew - quad.form.inv(K_joint, Kbar)
        # if (!matrixcalc::is.positive.definite(Sig_test_der)) {
        #     Sig_test_der <- as.matrix(nearPD(Sig_test_der)$mat)
        # }
        Sig_test_der <- as.matrix(nearPD(Sig_test_der)$mat)
        R <- chol(K_joint)
        # B <- forwardsolve(t(R), yJ)
        # aa <- backsolve(R, B)
        # mu_test_der <- crossprod(Kbar, aa)
        mu_test_der <- crossprod(Kbar, backsolve(R, forwardsolve(t(R), yJ)))
        # mean_var <- cond_norm_mean_var(mean_vec_1 = 0, mean_vec_2 = 0, 
        #                                obs_2 = yJ, cov_mat_1 = Kfnewfnew, 
        #                                cov_mat_2 = K_joint,
        #                                cov_mat_12 = t(Kbar))
        # predictive_f_der[i, ] <- mvnfast::rmvn(1, mu = mean_var$condMean, 
        #                             sigma = mean_var$condVar)
        predictive_f_der[i, ] <- rmvn(1, mu = mu_test_der, sigma = Sig_test_der)

    }

    PI_f_der <- apply(predictive_f_der, 2, quantile, prob = c(0.025, 0.975))
    mu_test_der_mean <- apply(predictive_f_der, 2, function(x) {sum(x) / mm})

    return(list(pred_f = predictive_f_der, ci_low = PI_f_der[1, ],
                ci_high = PI_f_der[2, ],
                mu_test = mu_test_der_mean))
}


get_pi_t <- function(yJ, x, x_test, idx_der, sig, tau, h) {

    len_test <- length(x_test)
    mm <- dim(idx_der)[1]
    
    Kff <- compute_cov_1d(idx1 = x, idx2 = x, tau = tau, h = h)
    Kffnew <- compute_cov_1d(idx1 = x, idx2 = x_test, tau = tau, h = h)
    Kfnewfnew <- compute_cov_1d(idx1 = x_test, idx2 = x_test, tau = tau, h = h)
    
    # K_joint_lst <- lapply(1:mm, function(i) {
    #     K_joint <- compute_joint_cov(idx_obs = x, idx_der = idx_der[i, ],
    #                                  sig = sig, tau = tau, h = h)
    #     # K_joint <- as.matrix(nearPD(K_joint)$mat)
    #     return(K_joint)
    # })
    # 
    # Kbar_lst <- lapply(1:mm, function(i) {
    #     Kdfnew <- computeCovDer1(idx1 = idx_der[i, ], idx2 = x_test, 
    #                              tau = tau, h = h)
    #     return(rbind(Kffnew, Kdfnew))
    # })
    # 
    # mean_var_lst <- lapply(1:mm, function(j) {
    #     return(cond_norm_mean_var(mean_vec_1 = 0, mean_vec_2 = 0, obs_2 = yJ,
    #                        cov_mat_1 = Kfnewfnew, cov_mat_2 = K_joint_lst[[j]],
    #                        cov_mat_12 = t(Kbar_lst[[j]])))
    # })
    
    
    mean_var_lst <- lapply(1:mm, function(i) {
        K_joint <- compute_joint_cov(idx_obs = x, idx_der = idx_der[i, ],
                                     sig = sig, tau = tau, h = h)
        Kdfnew <- computeCovDer1(idx1 = idx_der[i, ], idx2 = x_test, 
                                 tau = tau, h = h)
        # Kbar <- rbind(Kffnew, Kdfnew)
        return(cond_norm_mean_var(mean_vec_1 = 0, mean_vec_2 = 0, obs_2 = yJ,
                                  cov_mat_1 = Kfnewfnew, cov_mat_2 = K_joint,
                                  cov_mat_12 = t(rbind(Kffnew, Kdfnew))))
    })
    
    # mean_vec_1 = 0
    # mean_vec_2 = 0
    # obs_2 = yJ
    # cov_mat_1 = Kfnewfnew
    # cov_mat_2 = K_joint_lst[[j]]
    # cov_mat_12 = t(Kbar_lst[[j]])
    # 
    # B <- cov_mat_12 %*% chol2inv(chol(cov_mat_2))
    # cMu <- mean_vec_1 + B %*% (obs_2 - mean_vec_2)
    # # cVar <- cov_mat_1 - B %*% t(cov_mat_12)
    # cVar <- cov_mat_1 - tcrossprod(B, cov_mat_12)
    # if (!matrixcalc::is.positive.definite(cVar)) {
    #     cVar <- as.matrix(Matrix::nearPD(cVar)$mat)
    # }
    pred_f_der_mat <- sapply(mean_var_lst, function(x) {
        mvnfast::rmvn(1, mu = x$condMean,
                      sigma = x$condVar)})
    
    # pred_f_der_mat <- sapply(1:mm, function(i) {
    #     K_joint <- compute_joint_cov(idx_obs = x, idx_der = idx_der[i, ], 
    #                                  sig = sig, tau = tau, h = h)
    #     Kdfnew <- computeCovDer1_h(idx_der[i, ], idx2 = x_test, tau = tau, 
    #                                h = h) 
    #     Kbar <- rbind(Kffnew, Kdfnew)
    #     
    #     mean_var <- cond_norm_mean_var(mean_vec_1 = 0, mean_vec_2 = 0, 
    #                                    obs_2 = yJ, cov_mat_1 = Kfnewfnew, 
    #                                    cov_mat_2 = K_joint,
    #                                    cov_mat_12 = t(Kbar))
    #     
    #     pred_f_der <- mvnfast::rmvn(1, mu = mean_var$condMean, 
    #                                 sigma = mean_var$condVar)
    # })
    
    
    # Kbar_lst <- apply(idx_der, 1, function(i) {
    #     Kdfnew <- computeCovDer1_h(idx1 = i, idx2 = x_test, tau, h) 
    #     return(rbind(Kffnew, Kdfnew))
    # })
    # 
    # mean_var_lst <- lapply(1:dim(idx_der)[1], function(j) {
    #     return(cond_norm_mean_var(mean_vec_1 = 0, mean_vec_2 = 0, obs_2 = yJ, 
    #                        cov_mat_1 = Kfnewfnew, cov_mat_2 = K_joint_lst[[j]],
    #                        cov_mat_12 = t(Kbar_lst[[j]])))
    # })
    # 
    # pred_f_der_lst <- lapply(1:dim(idx_der)[1], function(j) {
    #     mvnfast::rmvn(1, mu = mean_var_lst[[j]]$condMean, 
    #                   sigma = mean_var_lst[[j]]$condVar)})
    
    # for (i in 1:dim(idx_der)[1]) {
    #     if (i %% 50 == 0) print(i)
    #     
    # 
    #     # Kdf <- computeCovDer1_h(idx1 = idx_der[i, ], idx2 = x, tau, h) 
    #     # Kdd <- computeCovDer2_h(idx1 = idx_der[i, ], tau, h)
    #     # K_joint <- rbind(cbind(Kff + diag(sig ^ 2, length(x)), t(Kdf)), 
    #     #                  cbind(Kdf, Kdd))
    #     
    # 
    #     # if(!matrixcalc::is.positive.definite(K_joint)){
    #     #     K_joint <- as.matrix(Matrix::nearPD(K_joint)$mat)
    #     # }
    #     
    #     # KJ_Kbar <- quad.form.inv(K_joint, Kbar) 
    #     # Sig_test_der <- Kfnewfnew - KJ_Kbar
    #     # 
    #     # Sig_test_der <- as.matrix(Matrix::nearPD(Sig_test_der)$mat)
    #     # 
    #     # R <- chol(K_joint)
    #     # B <- forwardsolve(t(R), yJ)
    #     # aa <- backsolve(R, B)
    #     # mu_test_der <- crossprod(Kbar, aa)
    #     
    #     
    #     
    #     # predict_mean_var_der <- cond_norm_mean_var(mean_vec_1 = 0, 
    #     #                                            mean_vec_2 = 0,
    #     #                                            obs_2 = yJ, 
    #     #                                            cov_mat_1 = Kfnewfnew,
    #     #                                            cov_mat_2 = K_joint,
    #     #                                            cov_mat_12 = t(Kbar))
    #     
    #     
    #     # predictive_f_der[i, ] <- mvnfast::rmvn(1, mu = mu_test_der, 
    #     #                                        sigma = Sig_test_der)
    # 
    # }
    
    PI_f_der <- apply(pred_f_der_mat, 1, quantile, prob = c(0.025, 0.975))
    mu_test_der_mean <- apply(pred_f_der_mat, 1, mean)
    
    return(list(pred_f = pred_f_der_mat, ci_low = PI_f_der[1, ],
                ci_high = PI_f_der[2, ],
                mu_test = mu_test_der_mean))
}


get_pi_t_matern <- function(yJ, x, x_test, idx_der, sig, tau, l, nu = 1.5) {
    
    len_test <- length(x_test)
    mm <- dim(idx_der)[1]
    
    Kff <- compute_cov_1d(idx1 = x, idx2 = x, ker_fcn = matern_ker, 
                          tau = tau, l = l, nu = nu)
    Kffnew <- compute_cov_1d(idx1 = x, idx2 = x_test, ker_fcn = matern_ker,
                             tau = tau, l = l, nu = nu)
    Kfnewfnew <- compute_cov_1d(idx1 = x_test, idx2 = x_test, ker_fcn = matern_ker,
                                tau = tau, l = l, nu = nu)

    mean_var_lst <- lapply(1:mm, function(i) {
        K_joint <- compute_joint_cov_matern(idx_obs = x, idx_der = idx_der[i, ],
                                     sig = sig, tau = tau, l = l, nu = nu)
        Kdfnew <- computeCovDer1(idx1 = idx_der[i, ], idx2 = x_test, 
                                 ker_der1_fcn = matern_ker_der1,
                                 tau = tau, l = l, nu = nu)
        # Kbar <- rbind(Kffnew, Kdfnew)
        return(cond_norm_mean_var(mean_vec_1 = 0, mean_vec_2 = 0, obs_2 = yJ,
                                  cov_mat_1 = Kfnewfnew, cov_mat_2 = K_joint,
                                  cov_mat_12 = t(rbind(Kffnew, Kdfnew))))
    })

    pred_f_der_mat <- sapply(mean_var_lst, function(x) {
        mvnfast::rmvn(1, mu = x$condMean,
                      sigma = x$condVar)})

    PI_f_der <- apply(pred_f_der_mat, 1, quantile, prob = c(0.025, 0.975))
    mu_test_der_mean <- apply(pred_f_der_mat, 1, mean)
    
    return(list(pred_f = pred_f_der_mat, ci_low = PI_f_der[1, ],
                ci_high = PI_f_der[2, ],
                mu_test = mu_test_der_mean))
}


# get_pred_ci_t(yJ, x, x_test, idx_der, sig, tau, phi = NULL, h)
# 
# j = 1
# 
# system.time(pred_old <- get_pred_ci_t(yJ = c(YY[[j]]$y, rep(0, 2)),
#               x = YY[[j]]$x, 
#               x_test = x_test, 
#               idx_der = sample_t_EM_multi_lst[[j]],
#               sig = EB_EM_multi_lst[[j]][1], 
#               tau = EB_EM_multi_lst[[j]][2],
#               h = EB_EM_multi_lst[[j]][3]))
# 
# system.time(pred_new <- get_pi_t(yJ = c(YY[[j]]$y, rep(0, 2)),
#                                       x = YY[[j]]$x, 
#                                       x_test = x_test, 
#                                       idx_der = sample_t_EM_multi_lst[[j]],
#                                       sig = EB_EM_multi_lst[[j]][1], 
#                                       tau = EB_EM_multi_lst[[j]][2],
#                                       h = EB_EM_multi_lst[[j]][3]))



# microbenchmark(pred_old <- get_pred_ci_t(yJ = c(YY[[j]]$y, rep(0, 2)),
#                                          x = YY[[j]]$x, 
#                                          x_test = x_test, 
#                                          idx_der = sample_t_EM_multi_lst[[j]],
#                                          sig = EB_EM_multi_lst[[j]][1], 
#                                          tau = EB_EM_multi_lst[[j]][2],
#                                          h = EB_EM_multi_lst[[j]][3]),
#                pred_new <- get_pi_t(yJ = c(YY[[j]]$y, rep(0, 2)),
#                                          x = YY[[j]]$x, 
#                                          x_test = x_test, 
#                                          idx_der = sample_t_EM_multi_lst[[j]],
#                                          sig = EB_EM_multi_lst[[j]][1], 
#                                          tau = EB_EM_multi_lst[[j]][2],
#                                          h = EB_EM_multi_lst[[j]][3]), 
#                times = 5)


# yJ = c(YY[[j]]$y, rep(0, 2))
# x = YY[[j]]$x
# x_test = x_test
# idx_der = sample_t_EM_multi_lst[[j]]
# sig = EB_EM_multi_lst[[j]][1]
# tau = EB_EM_multi_lst[[j]][2]
# h = EB_EM_multi_lst[[j]][3]

# system.time(pred_new <- get_pi_t(yJ = c(YY[[k]]$y, rep(0, 2)),
#          x = YY[[k]]$x, 
#          x_test = x_test, 
#          idx_der = sample_t_EM_multi_lst[[k]],
#          sig = EB_EM_multi_lst[[k]][1], 
#          tau = EB_EM_multi_lst[[k]][2],
#          h = EB_EM_multi_lst[[k]][3]))
# 
# system.time(pred_old <- get_pred_ci_t(yJ = c(YY[[k]]$y, rep(0, 2)),
#                                  x = YY[[k]]$x, 
#                                  x_test = x_test, 
#                                  idx_der = sample_t_EM_multi_lst[[k]],
#                                  sig = EB_EM_multi_lst[[k]][1], 
#                                  tau = EB_EM_multi_lst[[k]][2],
#                                  h = EB_EM_multi_lst[[k]][3]))


