## original marg_lik_gp_der_mc_new_h_1()


# log_mar_lik_gp_der_mc <- function(theta, y, x_vec, der_mc_mat, H0,
#                                        ga_shape = 6, ga_rate = 5,
#                                        is.sig.par = TRUE) {
#     nn <- ncol(der_mc_mat)
#     n <- length(y)
#     mu_vec <- rep(0L, n)
# 
#     Kff <- power_expo_h(H0, theta[2], theta[3])
#     
#     Kdf_mc <- apply(der_mc_mat, 1, computeCovDer1_h, idx2 = x_vec, 
#                      tau = theta[2], h = theta[3])
#     
#     den_value_mc <- apply(Kdf_mc, 2, function(x) {
#         Sigma <- Kff - tcrossprod(x) / ((theta[2] / theta[3]) ^ 2)
#         mvnfast::dmvn(y, mu = mu_vec, sigma = Sigma + diag(theta[1] ^ 2, n))
#     })
#     A <- sum(den_value_mc) / length(den_value_mc)
# 
#     if (is.sig.par) {
#         B <- invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape, 
#                                  rate = ga_rate, log = TRUE)
#         return(-log(A) - B)
#     } else {
#         return(-log(A))
#     }
# }

marg_lik_gp_der_mc_new_h_1 <- function(theta, y, x_vec, der_mc_mat, H0,
                                       ga_shape = 6, ga_rate = 5,
                                       is.sig.par = TRUE) {
    # dim_der <- dim(der_mc_mat)
    nn <- ncol(der_mc_mat)
    n <- length(y)
    # mu_vec <- rep(0, n)
    Kff <- se_ker(H0 = H0, tau = theta[2], h = theta[3])
    
    Kdf_lst <- as.list(data.frame(apply(der_mc_mat, 2, computeCovDer1_h, 
                                        idx2 = x_vec, 
                                        tau = theta[2], h = theta[3])))
    Kdf_lst_old <- lapply(Kdf_lst, function(x){
        matrix(x, nrow = nn, ncol = n)
    })
    
    # den_value_lst <- lapply(1:mm, function(m) {
    #     Sigma <- Kff - crossprod(Kdf_lst[[m]]) / ((theta[2] / theta[3]) ^ 2)
    #     diag(Sigma) <- diag(Sigma) + theta[1] ^ 2
    #     dmvn(y, mu = mu_vec, sigma = Sigma)
    # })
    den_value_lst_old <- lapply(Kdf_lst_old, function(m) {
        Sigma <- Kff_old - crossprod(m) / ((theta[2] / theta[3]) ^ 2)
        diag(Sigma) <- diag(Sigma) + theta[1] ^ 2
        dmvn(y, mu = rep(0, n), sigma = Sigma)
        # Sigma
    })
    # , log = TRUE
    A <- sum(unlist(den_value_lst)) / nrow(der_mc_mat)
    if (is.sig.par) {
        B <- invgamma::dinvgamma(theta[1]^2, shape = ga_shape, rate = ga_rate, log = TRUE)
        return(-log(A)-B)
    } else {
        return(-log(A))
    }
    # B <- dinvgamma(theta[1]^2, shape = 3, scale = 1, log = TRUE)
    # return(-log(A))
    # return(-A)
}


log_mar_lik_gp_der_mc <- function(theta, y, x_vec, der_mc_mat, H0,
                                  ga_shape = 6, ga_rate = 5,
                                  a_h = 1, b_h = 1,
                                  is.sig.par = TRUE,
                                  is.h.par = TRUE) {
    nn <- nrow(der_mc_mat)
    n <- length(y)
    # mu_vec <- rep(0L, n)
    
    Kff <- se_ker(H0 = H0, tau = theta[2], h = theta[3])
    # Kff <- compute_cov_1d(idx1 = x_vec, tau = theta[2], h = theta[3])
    
    # Kdf_mc <- apply(der_mc_mat, 2, computeCovDer1, idx2 = x_vec, 
    #                 tau = theta[2], h = theta[3])
    
    
    Kdf_lst <- as.list(data.frame(apply(der_mc_mat, 2, computeCovDer1, 
                                        idx2 = x_vec, 
                                        tau = theta[2], h = theta[3])))
    
    Kdf_lst_new <- lapply(Kdf_lst, function(x){
        matrix(x, nrow = nn, ncol = n)
    })
    
    den_value_lst_new <- lapply(Kdf_lst_new, function(m) {
        # m <- matrix(m, nrow = nn, ncol = n)
        Sigma <- Kff - crossprod(m) / ((theta[2] / theta[3]) ^ 2)
        Psi <- Sigma + diag(theta[1] ^ 2, n)
        
        Psi <- (Psi + t(Psi)) / 2
        if (!is.positive.definite(Psi)) {
            Psi <- as.matrix(nearPD(Psi)$mat)
        }
        mvnfast::dmvn(y, mu = rep(0L, n), sigma = Psi)
        # Sigma <- Sigma + diag(theta[1] ^ 2, n)
        # Sigma
    })
    
    
    # den_value_mc <- apply(Kdf_mc, 2, function(x) {
    #     x <- matrix(x, nrow = nn, ncol = n)
    #     Sigma <- Kff - crossprod(x) / ((theta[2] / theta[3]) ^ 2)
    #     mvnfast::dmvn(y, mu = mu_vec, sigma = Sigma + diag(theta[1] ^ 2, n))
    # })
    
    # A <- sum(den_value_mc) / length(den_value_mc)
    
    # A <- sum(unlist(den_value_lst)) / ncol(der_mc_mat)
    if (is.sig.par) {
        B <- invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape, 
                                 rate = ga_rate, log = TRUE)
        return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - B)
    } else {
        return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)))
    }
    
    
    # if (is.sig.par) {
    #     # B <- invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape, 
    #     #                          rate = ga_rate, log = TRUE)
    #     return(-log(sum(den_value_mc) / length(den_value_mc)) - 
    #                invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape, 
    #                                    rate = ga_rate, log = TRUE))
    # } else {
    #     return(-log(sum(den_value_mc) / length(den_value_mc)))
    # }
}


# res <- Rsolnp::solnp(pars = theta_k, fun = log_mar_lik_gp_der_mc,
#                      LB = lower, control = ctrl, y = y, x_vec = x, 
#                      H0 = H0, der_mc_mat = sample_t,
#                      ga_shape = ga_shape, ga_rate = ga_rate,
#                      is.sig.par = is.sig.par,
#                      is.h.par = is.h.par,
#                      a_h = a_h, b_h = b_h)
# 
# der_mc_mat = sample_t
# theta <- theta_k

log_mar_lik_gp_der_mc <- function(theta, y, x_vec, der_mc_mat, H0,
                                  ga_shape = 6, ga_rate = 5,
                                  a_h = 1, b_h = 1,
                                  is.sig.par = TRUE,
                                  is.h.par = TRUE) {
    nn <- nrow(der_mc_mat)
    n <- length(y)
    # mu_vec <- rep(0L, n)
    
    Kff <- se_ker(H0 = H0, tau = theta[2], h = theta[3])
    # Kff <- compute_cov_1d(idx1 = x_vec, tau = theta[2], h = theta[3])
    
    # Kdf_mc <- apply(der_mc_mat, 2, computeCovDer1, idx2 = x_vec, 
    #                 tau = theta[2], h = theta[3])
    
    
    Kdf_lst <- as.list(data.frame(apply(der_mc_mat, 2, computeCovDer1, 
                                        idx2 = x_vec, 
                                        tau = theta[2], h = theta[3])))
    
    Kdf_lst_new <- lapply(Kdf_lst, function(x){
        matrix(x, nrow = nn, ncol = n)
    })
    
    den_value_lst_new <- lapply(Kdf_lst_new, function(m) {
        # m <- matrix(m, nrow = nn, ncol = n)
        Sigma <- Kff - crossprod(m) / ((theta[2] / theta[3]) ^ 2)
        Psi <- Sigma + diag(theta[1] ^ 2, n)
        # Psi <- Sigma + diag(theta[1] ^ 2, n)
        if (!is.positive.definite(Psi)) {
            Psi <- as.matrix(Matrix::nearPD(Psi)$mat)
        }
        mvnfast::dmvn(y, mu = rep(0L, n), sigma = Psi)
        # Sigma <- Sigma + diag(theta[1] ^ 2, n)
        # Sigma
    })
    
    
    # den_value_mc <- apply(Kdf_mc, 2, function(x) {
    #     x <- matrix(x, nrow = nn, ncol = n)
    #     Sigma <- Kff - crossprod(x) / ((theta[2] / theta[3]) ^ 2)
    #     mvnfast::dmvn(y, mu = mu_vec, sigma = Sigma + diag(theta[1] ^ 2, n))
    # })
    
    # A <- sum(den_value_mc) / length(den_value_mc)
    
    # A <- sum(unlist(den_value_lst)) / ncol(der_mc_mat)
    if (is.sig.par) {
        B <- invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape, 
                                 rate = ga_rate, log = TRUE) + log(2 * theta[1])
        if (is.h.par) {
            C <- dgamma(theta[3], shape = a_h, rate = b_h, log = TRUE)
            return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - B - C)
        } else {
            return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - B)
        }
    } else {
        if (is.h.par) {
            C <- dgamma(theta[3], shape = a_h, rate = b_h, log = TRUE)
            return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - C)
        } else {
            return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)))
        }
        # return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)))
    }
}




log_mar_lik_gp_der_mc_multi <- function(theta, y, x_vec, der_mc_mat, H0,
                                        ga_shape = 6, ga_rate = 5,
                                        a_h = 1, b_h = 1,
                                        is.sig.par = TRUE,
                                        is.h.par = TRUE) {
    nn <- nrow(der_mc_mat)
    n <- length(y)
    
    Kff <- se_ker(H0 = H0, tau = theta[2], h = theta[3])
    Kdf_lst <- as.list(data.frame(apply(der_mc_mat, 2, computeCovDer1, 
                                        idx2 = x_vec, 
                                        tau = theta[2], h = theta[3])))
    
    Kdf_lst_new <- lapply(Kdf_lst, function(x){
        matrix(x, nrow = nn, ncol = n)
    })
    
    
    Kdd_lst <- as.list(data.frame(apply(der_mc_mat, 2, computeCovDer2,
                                        tau = theta[2], h = theta[3])))
    Kdd_lst_new <- lapply(Kdd_lst, function(x){
        matrix(x, nrow = nn, ncol = nn)
    })
    
    den_value_lst_new <- lapply(1:ncol(der_mc_mat), function(m) {
        Sigma <- Kff - quad.form.inv(Kdd_lst_new[[m]], Kdf_lst_new[[m]])
        Psi <- Sigma + diag(theta[1] ^ 2, n)
        Psi <- (Psi + t(Psi))/2
        if (!is.positive.definite(Psi)) {
            Psi <- as.matrix(Matrix::nearPD(Psi)$mat)
        }
        mvnfast::dmvn(y, mu = rep(0L, n), sigma = Psi)
    })
    
    # den_value_lst_new <- lapply(Kdf_lst_new, function(m) {
    #     # m <- matrix(m, nrow = nn, ncol = n)
    #     Sigma <- Kff - crossprod(m) / ((theta[2] / theta[3]) ^ 2)
    #     Psi <- Sigma + diag(theta[1] ^ 2, n)
    #     if (!is.positive.definite(Psi)) {
    #         Psi <- as.matrix(Matrix::nearPD(Psi)$mat)
    #     }
    #     mvnfast::dmvn(y, mu = rep(0L, n), sigma = Psi)
    # })
    
    if (is.sig.par) {
        B <- invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape, 
                                 rate = ga_rate, log = TRUE) + log(2 * theta[1])
        if (is.h.par) {
            C <- dgamma(theta[3], shape = a_h, rate = b_h, log = TRUE)
            return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - B - C)
        } else {
            return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - B)
        }
    } else {
        if (is.h.par) {
            C <- dgamma(theta[3], shape = a_h, rate = b_h, log = TRUE)
            return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - C)
        } else {
            return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)))
        }
    }
}




log_mar_lik_gp_der_mc_matern <- function(theta, y, x_vec, der_mc_mat, H0,
                                  ga_shape = 6, ga_rate = 5,
                                  a_h = 1, b_h = 1, nu,
                                  is.sig.par = TRUE,
                                  is.h.par = TRUE) {
    nn <- nrow(der_mc_mat)
    n <- length(y)
    Kff <- matern_ker(H0 = H0, nu = nu, tau = theta[2], l = theta[3])
    Kdf_lst <- as.list(data.frame(apply(der_mc_mat, 2, computeCovDer1, 
                                        idx2 = x_vec, 
                                        ker_der1_fcn = matern_ker_der1,
                                        tau = theta[2], l = theta[3],
                                        nu = nu)))
    
    Kdf_lst_new <- lapply(Kdf_lst, function(x){
        matrix(x, nrow = nn, ncol = n)
    })
    
    den_value_lst_new <- lapply(Kdf_lst_new, function(m) {
        if (nu == 1.5) {
            Sigma <- Kff - crossprod(m) / (3 * theta[2] ^ 2 / (theta[3] ^ 2))
        } 
        if (nu == 2.5) {
            Sigma <- Kff - crossprod(m) / (5 * theta[2] ^ 2 / (3 * theta[3] ^ 2))
        } 
        if (nu == 0.5) {
            Sigma <- Kff + crossprod(m) / (theta[2] ^ 2 / theta[3] ^ 2)
        }
        
        Psi <- Sigma + diag(theta[1] ^ 2, n)
        if (!is.positive.definite(Psi)) {
            Psi <- as.matrix(Matrix::nearPD(Psi)$mat)
        }
        mvnfast::dmvn(y, mu = rep(0L, n), sigma = Psi)

    })

    if (is.sig.par) {
        B <- invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape, 
                                 rate = ga_rate, log = TRUE) + log(2 * theta[1])
        if (is.h.par) {
            C <- dgamma(theta[3], shape = a_h, rate = b_h, log = TRUE)
            return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - B - C)
        } else {
            return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - B)
        }
    } else {
        if (is.h.par) {
            C <- dgamma(theta[3], shape = a_h, rate = b_h, log = TRUE)
            return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - C)
        } else {
            return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)))
        }
    }
}





log_mar_lik_gp_der_mc_multi_matern <- function(theta, y, x_vec, der_mc_mat, H0,
                                        ga_shape = 6, ga_rate = 5,
                                        a_h = 1, b_h = 1, nu = 1.5,
                                        is.sig.par = TRUE,
                                        is.h.par = TRUE) {
    nn <- nrow(der_mc_mat)
    n <- length(y)
    
    Kff <- matern_ker(H0 = H0, tau = theta[2], l = theta[3], nu = nu)
    Kdf_lst <- as.list(data.frame(apply(der_mc_mat, 2, computeCovDer1, 
                                        idx2 = x_vec, 
                                        ker_der1_fcn = matern_ker_der1,
                                        tau = theta[2], l = theta[3], nu = nu)))
    
    Kdf_lst_new <- lapply(Kdf_lst, function(x){
        matrix(x, nrow = nn, ncol = n)
    })
    
    
    Kdd_lst <- as.list(data.frame(apply(der_mc_mat, 2, computeCovDer2,
                                        ker_der2_fcn = matern_ker_der2,
                                        tau = theta[2], l = theta[3],
                                        nu = nu)))
    Kdd_lst_new <- lapply(Kdd_lst, function(x){
        matrix(x, nrow = nn, ncol = nn)
    })
    
    den_value_lst_new <- lapply(1:ncol(der_mc_mat), function(m) {
        Sigma <- Kff - quad.form.inv(Kdd_lst_new[[m]], Kdf_lst_new[[m]])
        Psi <- Sigma + diag(theta[1] ^ 2, n)
        Psi <- (Psi + t(Psi))/2
        if (!is.positive.definite(Psi)) {
            Psi <- as.matrix(Matrix::nearPD(Psi)$mat)
        }
        mvnfast::dmvn(y, mu = rep(0L, n), sigma = Psi)
    })
    
    # den_value_lst_new <- lapply(Kdf_lst_new, function(m) {
    #     # m <- matrix(m, nrow = nn, ncol = n)
    #     Sigma <- Kff - crossprod(m) / ((theta[2] / theta[3]) ^ 2)
    #     Psi <- Sigma + diag(theta[1] ^ 2, n)
    #     if (!is.positive.definite(Psi)) {
    #         Psi <- as.matrix(Matrix::nearPD(Psi)$mat)
    #     }
    #     mvnfast::dmvn(y, mu = rep(0L, n), sigma = Psi)
    # })
    
    if (is.sig.par) {
        B <- invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape, 
                                 rate = ga_rate, log = TRUE) + log(2 * theta[1])
        if (is.h.par) {
            C <- dgamma(theta[3], shape = a_h, rate = b_h, log = TRUE)
            return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - B - C)
        } else {
            return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - B)
        }
    } else {
        if (is.h.par) {
            C <- dgamma(theta[3], shape = a_h, rate = b_h, log = TRUE)
            return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - C)
        } else {
            return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)))
        }
    }
}






log_mar_lik_gp_der_mc_new <- function(theta, y, x_vec, der_mc_mat, H0,
                                      ga_shape = 6, ga_rate = 5,
                                      a_h = 1, b_h = 1,
                                      is.sig.par = TRUE,
                                      is.h.par = TRUE) {
    nn <- nrow(der_mc_mat)
    n <- length(y)
    # mu_vec <- rep(0L, n)
    
    Kff <- se_ker(H0 = H0, tau = theta[2], h = theta[3])
    # Kff <- compute_cov_1d(idx1 = x_vec, tau = theta[2], h = theta[3])
    
    # Kdf_mc <- apply(der_mc_mat, 2, computeCovDer1, idx2 = x_vec, 
    #                 tau = theta[2], h = theta[3])
    
    
    Kdf_lst <- as.list(data.frame(apply(der_mc_mat, 2, computeCovDer1,
                                        idx2 = x_vec,
                                        tau = theta[2], h = theta[3])))

    Kdf_lst_new <- lapply(Kdf_lst, function(x){
        matrix(x, nrow = nn, ncol = n)
    })
    
    # Kdf_lst_new <- lapply(der_mc_mat, 2, computeCovDer1, 
    #                      idx2 = x_vec, 
    #                      tau = theta[2], h = theta[3])
    den_value_lst_new <- lapply(Kdf_lst_new, function(m) {
        # m <- matrix(m, nrow = nn, ncol = n)
        Sigma <- Kff - crossprod(m) / ((theta[2] / theta[3]) ^ 2)
        Psi <- Sigma + diag(theta[1] ^ 2, n)
        mvnfast::dmvn(y, mu = rep(0L, n), sigma = Psi, log = TRUE)
        # Sigma <- Sigma + diag(theta[1] ^ 2, n)
        # Sigma
    })

    if (is.sig.par) {
        B <- invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape, 
                                 rate = ga_rate, log = TRUE) + log(2 * theta[1])
        if (is.h.par) {
            C <- dgamma(theta[3], shape = a_h, rate = b_h, log = TRUE)
            # return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - B - C)'
            return(-sum(unlist(den_value_lst_new)) / ncol(der_mc_mat) - B - C)
        } else {
            # return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - B)
            return(-sum(unlist(den_value_lst_new)) / ncol(der_mc_mat) - B)
        }
    } else {
        if (is.h.par) {
            C <- dgamma(theta[3], shape = a_h, rate = b_h, log = TRUE)
            # return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - C)
            return(-sum(unlist(den_value_lst_new)) / ncol(der_mc_mat) - C)
        } else {
            # return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)))
            return(-sum(unlist(den_value_lst_new)) / ncol(der_mc_mat))
        }
        # return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)))
    }
}


log_mar_lik_gp_der_mc_new2 <- function(theta, y, x_vec, der_mc_mat, H0,
                                      ga_shape = 6, ga_rate = 5,
                                      a_h = 1, b_h = 1,
                                      is.sig.par = TRUE,
                                      is.h.par = TRUE) {
    nn <- nrow(der_mc_mat)
    n <- length(y)
    # mu_vec <- rep(0L, n)
    
    Kff <- se_ker(H0 = H0, tau = theta[2], h = theta[3])
    # Kff <- compute_cov_1d(idx1 = x_vec, tau = theta[2], h = theta[3])
    
    # Kdf_mc <- apply(der_mc_mat, 2, computeCovDer1, idx2 = x_vec, 
    #                 tau = theta[2], h = theta[3])
    
    # 
    # Kdf_lst <- as.list(data.frame(apply(der_mc_mat, 2, computeCovDer1,
    #                                     idx2 = x_vec,
    #                                     tau = theta[2], h = theta[3])))
    # 
    # Kdf_lst_new <- lapply(Kdf_lst, function(x){
    #     matrix(x, nrow = nn, ncol = n)
    # })
    # 
    Kdf_lst_new <- lapply(der_mc_mat, computeCovDer1,
                         idx2 = x_vec,
                         tau = theta[2], h = theta[3])
    
    den_value_lst_new <- lapply(Kdf_lst_new, function(m) {
        # m <- matrix(m, nrow = nn, ncol = n)
        Sigma <- Kff - crossprod(m) / ((theta[2] / theta[3]) ^ 2)
        Psi <- Sigma + diag(theta[1] ^ 2, n)
        mvnfast::dmvn(y, mu = rep(0L, n), sigma = Psi, log = TRUE)
        # Sigma <- Sigma + diag(theta[1] ^ 2, n)
        # Sigma
    })
    
    if (is.sig.par) {
        B <- invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape, 
                                 rate = ga_rate, log = TRUE) + log(2 * theta[1])
        if (is.h.par) {
            C <- dgamma(theta[3], shape = a_h, rate = b_h, log = TRUE)
            # return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - B - C)'
            return(-sum(unlist(den_value_lst_new)) / ncol(der_mc_mat) - B - C)
        } else {
            # return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - B)
            return(-sum(unlist(den_value_lst_new)) / ncol(der_mc_mat) - B)
        }
    } else {
        if (is.h.par) {
            C <- dgamma(theta[3], shape = a_h, rate = b_h, log = TRUE)
            # return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)) - C)
            return(-sum(unlist(den_value_lst_new)) / ncol(der_mc_mat) - C)
        } else {
            # return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)))
            return(-sum(unlist(den_value_lst_new)) / ncol(der_mc_mat))
        }
        # return(-log(sum(unlist(den_value_lst_new)) / ncol(der_mc_mat)))
    }
}

log_mar_lik_gp_der_mc_2_par <- function(theta, y, x_vec, der_mc_mat, H0,
                                        a_h = 1, b_h = 1, sig2, 
                                        is.h.par = TRUE) {
    nn <- nrow(der_mc_mat)
    n <- length(y)
    mu_y <- rep(0L, n)
    Kff_no_tau <- se_ker_no_tau(H0 = H0, h = theta[2])
    A <- Kff_no_tau + diag(1 / theta[1], n)
    
    # Kdf_lst_no_tau <- as.list(data.frame(apply(der_mc_mat, 2, computeCovDer1, 
    #                                            idx2 = x,
    #                                            ker_der1_fcn = se_ker_der1_no_tau,
    #                                            h = theta[2])))
    # 
    # Kdf_lst_no_tau <- lapply(Kdf_lst_no_tau, function(x){
    #     matrix(x, nrow = nn, ncol = n)
    # })
    
    Kdf_lst_no_tau <- lapply(der_mc_mat, computeCovDer1, 
                             idx2 = x,
                             ker_der1_fcn = se_ker_der1_no_tau,
                             h = theta[2])
    
    B_lst_no_tau <- lapply(Kdf_lst_no_tau, create_B, h = theta[2])
    
    den_value <- sapply(B_lst_no_tau, function(B) {
        # m <- matrix(m, nrow = nn, ncol = n)
        # Sigma <- Kff - crossprod(m) / ((theta[2] / theta[3]) ^ 2)
        Psi <- sig2 * theta[1] * (A + B)
        mvnfast::dmvn(y, mu = mu_y, sigma = Psi, log = TRUE)
        # Sigma <- Sigma + diag(theta[1] ^ 2, n)
        # Sigma
    })
    return(sum(den_value))
}


# if (is.sig.par) {
#     B <- invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape, 
#                              rate = ga_rate, log = TRUE) + log(2 * theta[1])
#     if (is.h.par) {
#         C <- dgamma(theta[3], shape = a_h, rate = b_h, log = TRUE)
#         return(-log_summ - B - C)
#     } else {
#         return(-log_summ - B)
#     }
# } else {
#     if (is.h.par) {
#         C <- dgamma(theta[3], shape = a_h, rate = b_h, log = TRUE)
#         return(-log_summ - C)
#     } else {
#         return(-log_summ)
#     }
#     # return(-log_summ)
# }

# Reduce("+", your.list) 
# test_val1 <- 1.4 + rnorm(100, 0, 0.1)
# test_val2 <- 1.4 + rnorm(200, 0, 0.1)
# 
# test_mat1 <- matrix(test_val1, 1, 100)
# test_mat2 <- matrix(test_val2, 2, 100)
# 
# 
# 
# test <- apply(test_mat2, 2, computeCovDer1_h,
#       idx2 = YY[[1]]$x, tau = 1, h = 1)
# 
# # test_lst0 <- as.list(test)
# 
# 
# test_lst <- (as.list(data.frame(apply(test_mat2, 2, computeCovDer1_h,
#                                       idx2 = YY[[1]]$x, tau = 1, h = 1))))
# 
# # test_lst1 <- lapply(test_lst, function(x) matrix(x, nrow = 1, ncol = 50))
# 
# 
# microbenchmark(test <- apply(test_mat2, 2, computeCovDer1_h,
#                              idx2 = YY[[1]]$x, tau = 1, h = 1),
#                test_lst <- (as.list(data.frame(apply(t(test_mat2), 1, computeCovDer1_h,
#                                                      idx2 = YY[[1]]$x, tau = 1, h = 1)))), times = 100)
# # 
# # test_apply <- apply(test, 2, function(x) {
#     x <- matrix(x, nrow = 2, ncol = 50)
#     Sigma <- Kff - crossprod(x)
#     mu_vec <- rep(0L, 50)
#     Psi <- Sigma + diag(1, 50)
#     mvnfast::dmvn(y, mu = mu_vec, sigma = Psi)
# })
# # 
# microbenchmark(test_apply <- apply(test, 2, function(x) {
#     x <- matrix(x, nrow = 2, ncol = 50)
#     Sigma <- Kff - crossprod(x)
#     mu_vec <- rep(0L, 50)
#     Psi <- Sigma + diag(1, 50)
#     # mvnfast::dmvn(y, mu = mu_vec, sigma = Psi)
# }),
# den_value_lst <- lapply(test_lst, function(m) {
#     m <- matrix(m, nrow = 2, ncol = 50)
#     Sigma <- Kff - crossprod(m)
#     mu_vec <- rep(0L, 50)
#     Psi <- Sigma + diag(1, 50)
#     # mvnfast::dmvn(y, mu = mu_vec, sigma = Psi)
# }), times = 10)


# microbenchmark({Kdf_lst <- as.list(data.frame(apply(der_mc_mat, 1, computeCovDer1_h, idx2 = x_vec, 
#                                                    tau = theta[2], h = theta[3])))
#                Kdf_lst <- lapply(Kdf_lst, function(x){
#                    matrix(x, nrow = nn, ncol = n)
#                })
#                den_value_lst <- lapply(Kdf_lst, function(m) {
#                    Sigma <- Kff - crossprod(m) / ((theta[2] / theta[3]) ^ 2)
#                    diag(Sigma) <- diag(Sigma) + theta[1] ^ 2
#                    # dmvn(y, mu = mu_vec, sigma = Sigma)
#                })},
#                {Kdf_mc <- apply(der_mc_mat, 2, computeCovDer1, idx2 = x_vec, 
#                                 tau = theta[2], h = theta[3])
#                den_value_mc <- apply(Kdf_mc, 2, function(x) {
#                    x <- matrix(x, nrow = nn, ncol = n)
#                    Sigma <- Kff - crossprod(x) / ((theta[2] / theta[3]) ^ 2)
#                    Sigma <- Sigma + diag(theta[1] ^ 2, n)
#                    # mvnfast::dmvn(y, mu = mu_vec, sigma = Sigma + diag(theta[1] ^ 2, n))
#                })},
#                times = 100)



# Kdf_lst <- as.list(data.frame(apply(der_mc_mat, 1, computeCovDer1_h, idx2 = x_vec, 
#                                     tau = theta[2], h = theta[3])))
# Kdf_lst <- lapply(Kdf_lst, function(x){
#     matrix(x, nrow = nn, ncol = n)
# })
# den_value_lst <- lapply(Kdf_lst, function(m) {
#     Sigma <- Kff - crossprod(m) / ((theta[2] / theta[3]) ^ 2)
#     diag(Sigma) <- diag(Sigma) + theta[1] ^ 2
#     # dmvn(y, mu = mu_vec, sigma = Sigma)
# })




# microbenchmark(sum(unlist(den_value_lst)) / 100, 
#                sum(test_apply) / 100)
# 
# 
# 
# den_value_lst <- lapply(test_lst, function(m) {
#     m <- matrix(m, nrow = 1, ncol = 50)
#     Sigma <- Kff - crossprod(m)
#     mu_vec <- rep(0L, 50)
#     Psi <- Sigma + diag(1, 50)
#     mvnfast::dmvn(y, mu = mu_vec, sigma = Psi)
# })
# 
# 
# microbenchmark(test_apply <- apply(test, 2, function(x) {
#     x <- matrix(x, nrow = 1, ncol = 50)
#     Sigma <- Kff - crossprod(x)
#     mu_vec <- rep(0L, 50)
#     Psi <- Sigma + diag(1, 50)
#     mvnfast::dmvn(y, mu = mu_vec, sigma = Psi)
# }),
# test_apply <- apply(test, 2, function(x) {
#     # x <- matrix(x, nrow = 1, ncol = 50)
#     Sigma <- Kff - tcrossprod(x)
#     mu_vec <- rep(0L, 50)
#     Psi <- Sigma + diag(1, 50)
#     mvnfast::dmvn(y, mu = mu_vec, sigma = Psi)
# }), times = 1000)
# 
# 
# microbenchmark(mean(test_apply), sum(test_apply) / length(100))
# system.time(for (i in 1:100000) {
#     X <- mvnfast::dmvn(YY[[k]]$y, mu = rep(0L, n), sigma = diag(n), log = TRUE)
# })
# 
# system.time(for (i in 1:100000) {
#     X <- mvnfast::dmvn(YY[[k]]$y, mu = rep(0L, n), sigma = diag(n))
# })
# 
# 
# system.time(
#     for(i in 1:1000000) Reduce(`+`, list(runif(100))))
# 
# system.time(
#     for(i in 1:1000000) sum(unlist(list(runif(100)))))
