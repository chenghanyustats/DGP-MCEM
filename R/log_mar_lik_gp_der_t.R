### original log_marginal_lik_gp_der_t_new_h()



# log_mar_lik_gp_der_t <- function(theta, y, x_vec, H0 = NULL,
#                                  ga_shape = 6, ga_rate = 5, is.sig.par = TRUE) {
#     # theta = c(sig, tau, h, t) (no mu)
#     mm <- length(theta)
#     if(is.null(H0)) {
#         H0 <- outer(as.vector(x_vec), as.vector(x_vec),
#                     FUN = function(x1, x2) (x1 - x2))
#     }
#     Kff <- power_expo_h(H0, theta[2], theta[3])
#     Kdf <- computeCovDer1_h(idx1 = theta[4:mm], idx2 = x_vec,
#                             tau = theta[2], h = theta[3])
#     Kdd <- computeCovDer2_h(idx1 = theta[4:mm], tau = theta[2], h = theta[3])
#     n <- length(x_vec)
#     Sigma <- Kff - emulator::quad.form.inv(Kdd, Kdf)
#     # Psi <- diag(theta[1] ^ 2, n) + Sigma
#     if(is.sig.par) {
#         return(-mvnfast::dmvn(y, mu = rep(0L, n),
#                               sigma = diag(theta[1] ^ 2, n) + Sigma,
#                               log = TRUE) -
#             invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape,
#                                 rate = ga_rate, log = TRUE))
#     } else {
#         return(-mvnfast::dmvn(y, mu = rep(0L, n),
#                               sigma = diag(theta[1] ^ 2, n) + Sigma,
#                               log = TRUE))
#     }
# }


log_mar_lik_gp_der_t <- function(theta, y, x_vec, H0, ga_shape = 6, ga_rate = 5, 
                                 is.sig.par = TRUE) {
    # theta = c(sig, tau, h, t) (no mu)
    mm <- length(theta)
    # if(is.null(H0)) {
    #     H0 <- outer(as.vector(x_vec), as.vector(x_vec), 
    #                 FUN = function(x1, x2) (x1 - x2))
    # }
    Kff <- se_ker(H0 = H0, tau = theta[2], h = theta[3])
    # Kff <- compute_cov_1d(idx1 = x_vec, tau = theta[2], h = theta[3])
    # Kdf <- computeCovDer1_test(idx1 = theta[4:mm], idx2 = x_vec,
    #                         tau = theta[2], h = theta[3])
    # Kdd <- computeCovDer2_test(idx1 = theta[4:mm], tau = theta[2], h = theta[3])
    Kdf <- computeCovDer1(idx1 = theta[4:mm], idx2 = x_vec,
                               tau = theta[2], h = theta[3])
    Kdd <- computeCovDer2(idx1 = theta[4:mm], tau = theta[2], h = theta[3])
    n <- length(x_vec)
    Sigma <- Kff - emulator::quad.form.inv(Kdd, Kdf)
    
    Psi <- diag(theta[1] ^ 2, n) + Sigma
    Psi <- (Psi + t(Psi)) / 2
    if (!is.positive.definite(Psi)) {
        Psi <- as.matrix(nearPD(Psi)$mat)
    }
    Psi <- as.matrix(Matrix::nearPD(Psi)$mat)
    
    if(is.sig.par) {
        return(-mvnfast::dmvn(y, mu = rep(0L, n),
                              sigma = Psi,
                              log = TRUE) -
                   invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape,
                                       rate = ga_rate, log = TRUE))
    } else {
        return(-mvnfast::dmvn(y, mu = rep(0L, n),
                              sigma = Psi,
                              log = TRUE))
    }
    
    
    # if(is.sig.par) {
    #     return(-mvnfast::dmvn(y, mu = rep(0L, n),
    #                           sigma = diag(theta[1] ^ 2, n) + Sigma,
    #                           log = TRUE) -
    #                invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape,
    #                                    rate = ga_rate, log = TRUE))
    # } else {
    #     return(-mvnfast::dmvn(y, mu = rep(0L, n),
    #                           sigma = diag(theta[1] ^ 2, n) + Sigma,
    #                           log = TRUE))
    # }
    # if(is.sig.par) {
    #     ll <- mvnfast::dmvn(y, mu = rep(0L, n), sigma = Psi, log = TRUE) +
    #         invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape, rate = ga_rate, log = TRUE)
    # } else {
    #     ll <- mvnfast::dmvn(y, mu = rep(0L, n), sigma = Psi, log = TRUE)
    # }
    
    # return(-dmvn(y, mu = rep(0, n), sigma = Psi, log = TRUE))
    # return(-ll)
}


log_mar_lik_gp_der_t_matern <- function(theta, y, x_vec, H0, nu,
                                        ga_shape = 6, ga_rate = 5, 
                                        is.sig.par = TRUE) {
    # theta = c(sig, tau, h, t) (no mu)
    mm <- length(theta)
    # if(is.null(H0)) {
    #     H0 <- outer(as.vector(x_vec), as.vector(x_vec), 
    #                 FUN = function(x1, x2) (x1 - x2))
    # }
    Kff <- matern_ker(H0 = H0, nu = nu, tau = theta[2], l = theta[3])
    # Kff <- compute_cov_1d(idx1 = x_vec, tau = theta[2], h = theta[3])
    # Kdf <- computeCovDer1_test(idx1 = theta[4:mm], idx2 = x_vec,
    #                         tau = theta[2], h = theta[3])
    # Kdd <- computeCovDer2_test(idx1 = theta[4:mm], tau = theta[2], h = theta[3])
    Kdf <- computeCovDer1(idx1 = theta[4:mm], idx2 = x_vec,
                          ker_der1_fcn = matern_ker_der1,
                          tau = theta[2], l = theta[3], nu = nu)
    Kdd <- computeCovDer2(idx1 = theta[4:mm], ker_der2_fcn = matern_ker_der2,
                          tau = theta[2], l = theta[3], nu = nu)
    n <- length(x_vec)
    Sigma <- Kff - emulator::quad.form.inv(Kdd, Kdf)
    Psi <- diag(theta[1] ^ 2, n) + Sigma
    # Psi <- as.matrix(Matrix::nearPD(Psi)$mat)
    if(is.sig.par) {
        return(-mvnfast::dmvn(y, mu = rep(0L, n),
                              sigma = Psi,
                              log = TRUE) -
                   invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape,
                                       rate = ga_rate, log = TRUE))
    } else {
        return(-mvnfast::dmvn(y, mu = rep(0L, n),
                              sigma = Psi,
                              log = TRUE))
    }
    # if(is.sig.par) {
    #     ll <- mvnfast::dmvn(y, mu = rep(0L, n), sigma = Psi, log = TRUE) +
    #         invgamma::dinvgamma(theta[1] ^ 2, shape = ga_shape, rate = ga_rate, log = TRUE)
    # } else {
    #     ll <- mvnfast::dmvn(y, mu = rep(0L, n), sigma = Psi, log = TRUE)
    # }
    
    # return(-dmvn(y, mu = rep(0, n), sigma = Psi, log = TRUE))
    # return(-ll)
}


# Kff <- power_expo_h(H0_diff[[1]], 1, 1)
# Kdf <- computeCovDer1_h(idx1 = cri_pts, idx2 = YY[[1]]$x, 
#                         tau = 1, h = 1) 
# Kdd <- computeCovDer2_h(idx1 = cri_pts, tau = 1, h = 1)
# 
# system.time(
#     for (i in 1:10000) {
#         emulator::quad.form.inv(Kdd, Kdf)
#     }
# )
# system.time(
#     for (i in 1:10000) {
#         t(Kdf) %*% chol2inv(chol(Kdd)) %*% Kdf
#     }
# )
# library(microbenchmark)
# microbenchmark(emulator::quad.form.inv(Kdd, Kdf), t(Kdf) %*% chol2inv(chol(Kdd)) %*% Kdf)
# microbenchmark(A <- diag(2, n) + Sigma, diag(Sigma) <- diag(Sigma) + 2)
# microbenchmark(power_expo_derivative2_h(H0_diff[[1]], 1, 1), power_expo_derivative2_h1(H0_diff[[1]], 1, 1),times = 1000)
# 
# 
# power_expo_derivative2_h <- function(idx12_diff, tau, h) {
#     # k11 <- power_expo_h(idx12_diff, tau, h)
#     phi <- 1 / h ^ 2
#     # k22 <- k11 * phi * (1 - phi * (idx12[1] - idx12[2]) ^ 2)
#     # k22 <- power_expo_h(idx12_diff, tau, h) * phi * (1 - phi * idx12_diff ^ 2)
#     k22 <- power_expo_h(idx12_diff, tau, h) * (phi - (phi * idx12_diff) ^ 2)
#     return(k22)
# }
# 
# power_expo_derivative2_h <- function(idx12_diff, tau, h) {
#     # k11 <- power_expo_h(idx12_diff, tau, h)
#     phi <- 1 / h ^ 2
#     # k22 <- k11 * phi * (1 - phi * (idx12[1] - idx12[2]) ^ 2)
#     # k22 <- power_expo_h(idx12_diff, tau, h) * phi * (1 - phi * idx12_diff ^ 2)
#     return(power_expo_h(idx12_diff, tau, h) * (phi - (phi * idx12_diff) ^ 2))
#     # return(k22)
# }


