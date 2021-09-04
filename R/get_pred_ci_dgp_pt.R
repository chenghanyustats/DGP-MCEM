### original get_pred_ci_pt()


# get_pred_ci_pt <- function(yJ, x, x_test, idx_der, sig, tau, h) {
#     if (is.null(phi)) {
#         phi <- 1 / h ^ 2
#     } 
#     ## treat yJ as c(y, f') = c(y , 0)
#     Kff <- computeCov_new(idx1 = x, idx2 = x, tau, h)
#     Kdf <- computeCovDer1_h(idx1 = idx_der, idx2 = x, tau, h) 
#     # Kfd <- t(Kdf)
#     Kdd <- computeCovDer2_h(idx1 = idx_der, idx2 = idx_der, tau, h)
#     K_joint <- rbind(cbind(Kff + sig ^ 2 * diag(length(x)), t(Kdf)), cbind(Kdf, Kdd))
#     Kffnew <- computeCov_new(idx1 = x, idx2 = x_test, tau, phi)
#     Kdfnew <- computeCovDer1_h(idx1 = idx_der, idx2 = x_test, tau, h) 
#     Kbar <- rbind(Kffnew, Kdfnew) 
#     Kfnewfnew <- computeCov_new(idx1 = x_test, idx2 = x_test, tau, phi)
#     if(!is.positive.definite(K_joint)){
#         K_joint <- as.matrix(nearPD(K_joint)$mat)
#     }
#     KJ_Kbar <- quad.form.inv(K_joint, Kbar) 
#     Sig_test_der <- Kfnewfnew - KJ_Kbar
#     Sig_test_der <- as.matrix(nearPD(Sig_test_der)$mat)
#     
#     Ky_star <- Kfnewfnew + diag(sig ^ 2, length(x_test)) 
#     Sig_test_y_der <- Ky_star - KJ_Kbar
#     Sig_test_y_der <- as.matrix(nearPD(Sig_test_y_der)$mat)
#     
#     R <- chol(K_joint) ## t(R) %*% R = K_joint
#     B <- forwardsolve(t(R), yJ) ## B = t(R)^(-1) * yJ
#     aa <- backsolve(R, B) ## aa = (t(R) %*% R) ^ (-1) yJ = K_joint ^ (-1) * yJ
#     mu_test_der <- crossprod(Kbar, aa)
#     
#     predictive_f_der <- rmvn(1000, mu = mu_test_der, sigma = Sig_test_der)
#     CI_Low_f_der <- apply(predictive_f_der, 2, quantile, prob = 0.025)
#     CI_High_f_der <- apply(predictive_f_der, 2, quantile, prob = 0.975)
#     
#     return(list(pred_f = predictive_f_der, ci_low = CI_Low_f_der,
#                 ci_high = CI_High_f_der,
#                 mu_test = mu_test_der))
# }


get_pred_ci_dgp_pt <- function(yJ, x, x_test, idx_der, sig, tau, h, 
                               n_draw = 1000) {
    ## treat yJ as c(y, f') = c(y , 0)
    # Kff <- compute_cov_1d(idx1 = x, idx2 = x, tau = tau, h = h)
    # Kdf <- computeCovDer1_h(idx1 = idx_der, idx2 = x, tau = tau, h = h) 
    # Kdd <- computeCovDer2_h(idx1 = idx_der, idx2 = idx_der, tau = tau, h = h)
    # K_joint <- rbind(cbind(Kff + sig ^ 2 * diag(length(x)), t(Kdf)), 
    #                  cbind(Kdf, Kdd))
    
    K_joint <- compute_joint_cov(idx_obs = x, idx_der = idx_der, sig = sig, 
                                 tau = tau, h = h)
    Kffnew <- compute_cov_1d(idx1 = x, idx2 = x_test, tau = tau, h = h)
    Kdfnew <- computeCovDer1_h(idx1 = idx_der, idx2 = x_test, tau = tau, h = h) 
    Kbar <- rbind(Kffnew, Kdfnew) 
    Kfnewfnew <- compute_cov_1d(idx1 = x_test, idx2 = x_test, tau = tau, h = h)
    # if(!is.positive.definite(K_joint)){
    #     K_joint <- as.matrix(nearPD(K_joint)$mat)
    # }
    KJ_Kbar <- emulator::quad.form.inv(K_joint, Kbar) 
    Sig_test_der <- Kfnewfnew - KJ_Kbar
    Sig_test_der <- (Sig_test_der + t(Sig_test_der)) / 2
    if (!is.positive.definite(Sig_test_der)) {
        Sig_test_der <- as.matrix(Matrix::nearPD(Sig_test_der)$mat)
    }
    
    # Ky_star <- Kfnewfnew + diag(sig ^ 2, length(x_test)) 
    # Sig_test_y_der <- Ky_star - KJ_Kbar
    # Sig_test_y_der <- as.matrix(Matrix::nearPD(Sig_test_y_der)$mat)
    
    R <- chol(K_joint) ## t(R) %*% R = K_joint
    # B <- forwardsolve(t(R), yJ) ## B = t(R)^(-1) * yJ
    # aa <- backsolve(R, B) ## aa = (t(R) %*% R) ^ (-1) yJ = K_joint ^ (-1) * yJ

    mu_test_der <- crossprod(Kbar, backsolve(R, forwardsolve(t(R), yJ)))
    
    predictive_f_der <- mvnfast::rmvn(n_draw, mu = mu_test_der, 
                                      sigma = Sig_test_der)
    PI_f_der <- apply(predictive_f_der, 2, stats::quantile, 
                      prob = c(0.025, 0.975))
    # CI_Low_f_der <- apply(predictive_f_der, 2, quantile, prob = 0.025)
    # CI_High_f_der <- apply(predictive_f_der, 2, quantile, prob = 0.975)
    
    return(list(pred_f = predictive_f_der, ci_low = PI_f_der[1, ],
                ci_high = PI_f_der[2, ],
                mu_test = mu_test_der))
}


get_pred_ci_dgp_pt_matern <- function(yJ, x, x_test, idx_der, sig, tau, l, 
                               n_draw = 1000, nu = 1.5) {
    ## treat yJ as c(y, f') = c(y , 0)
    # Kff <- compute_cov_1d(idx1 = x, idx2 = x, tau = tau, h = h)
    # Kdf <- computeCovDer1_h(idx1 = idx_der, idx2 = x, tau = tau, h = h) 
    # Kdd <- computeCovDer2_h(idx1 = idx_der, idx2 = idx_der, tau = tau, h = h)
    # K_joint <- rbind(cbind(Kff + sig ^ 2 * diag(length(x)), t(Kdf)), 
    #                  cbind(Kdf, Kdd))
    
    K_joint <- compute_joint_cov_matern(idx_obs = x, idx_der = idx_der, sig = sig, 
                                 tau = tau, l = l, nu = nu)
    Kffnew <- compute_cov_1d(idx1 = x, idx2 = x_test, ker_fcn = matern_ker,
                             tau = tau, l = l, nu = nu)
    Kdfnew <- computeCovDer1(idx1 = idx_der, idx2 = x_test, 
                             ker_der1_fcn = matern_ker_der1,
                             tau = tau, l = l, nu = nu) 
    Kbar <- rbind(Kffnew, Kdfnew) 
    Kfnewfnew <- compute_cov_1d(idx1 = x_test, idx2 = x_test, 
                                ker_fcn = matern_ker,
                                tau = tau, l = l, nu = nu)
    # if(!is.positive.definite(K_joint)){
    #     K_joint <- as.matrix(nearPD(K_joint)$mat)
    # }
    KJ_Kbar <- emulator::quad.form.inv(K_joint, Kbar) 
    Sig_test_der <- Kfnewfnew - KJ_Kbar
    Sig_test_der <- (Sig_test_der + t(Sig_test_der)) / 2
    if (!is.positive.definite(Sig_test_der)) {
        Sig_test_der <- as.matrix(Matrix::nearPD(Sig_test_der)$mat)
    }

    
    # Ky_star <- Kfnewfnew + diag(sig ^ 2, length(x_test)) 
    # Sig_test_y_der <- Ky_star - KJ_Kbar
    # Sig_test_y_der <- as.matrix(Matrix::nearPD(Sig_test_y_der)$mat)
    
    R <- chol(K_joint) ## t(R) %*% R = K_joint
    # B <- forwardsolve(t(R), yJ) ## B = t(R)^(-1) * yJ
    # aa <- backsolve(R, B) ## aa = (t(R) %*% R) ^ (-1) yJ = K_joint ^ (-1) * yJ
    
    mu_test_der <- crossprod(Kbar, backsolve(R, forwardsolve(t(R), yJ)))
    
    predictive_f_der <- mvnfast::rmvn(n_draw, mu = mu_test_der, 
                                      sigma = Sig_test_der)
    PI_f_der <- apply(predictive_f_der, 2, stats::quantile, 
                      prob = c(0.025, 0.975))
    # CI_Low_f_der <- apply(predictive_f_der, 2, quantile, prob = 0.025)
    # CI_High_f_der <- apply(predictive_f_der, 2, quantile, prob = 0.975)
    
    return(list(pred_f = predictive_f_der, ci_low = PI_f_der[1, ],
                ci_high = PI_f_der[2, ],
                mu_test = mu_test_der))
}
