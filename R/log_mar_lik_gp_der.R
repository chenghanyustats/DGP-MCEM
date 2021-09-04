### original log_marginal_lik_gp_der_new_h()


# log_mar_lik_gp_der <- function(theta, y, x_vec, der_vec,
#                                           H0 = NULL,
#                                           ga_shape = 6, ga_rate = 5,
#                                           is.sig.par = TRUE) {
#     # theta = c(sig, tau, h) (no mu)
#     if (is.null(H0)) {
#         H0 <- outer(as.vector(x_vec), as.vector(x_vec), 
#                     FUN = function(x1, x2) (x1 - x2))
#     }
#     Kff <- power_expo_h(H0, theta[2], theta[3])
#     Kdf <- computeCovDer1_h(idx1 = der_vec, idx2 = x_vec, 
#                             tau = theta[2], h = theta[3]) 
#     Kdd <- computeCovDer2_h(idx1 = der_vec, tau = theta[2], h = theta[3])
#     n <- length(x_vec)
#     Sigma <- Kff - emulator::quad.form.inv(Kdd, Kdf)
#     Psi <- diag(theta[1] ^ 2, n) + Sigma
#     if (is.sig.par) {
#         return(-mvnfast::dmvn(y, mu = rep(0L, n), sigma = Psi, log = TRUE) -
#             dinvgamma(theta[1] ^ 2, shape = ga_shape, rate = ga_rate, 
#                       log = TRUE))
#     } else {
#         return(-mvnfast::dmvn(y, mu = rep(0L, n), sigma = Psi, log = TRUE))
#     }
# }


log_mar_lik_gp_der <- function(theta, y, x_vec, der_vec, H0,
                               ga_shape = 6, ga_rate = 5,
                               is.sig.par = TRUE) {
    # theta = c(sig, tau, h) (no mu)
    # if (is.null(H0)) {
    #     H0 <- outer(as.vector(x_vec), as.vector(x_vec), 
    #                 FUN = function(x1, x2) (x1 - x2))
    # }
    Kff <- se_ker(H0 = H0, tau = theta[2], h = theta[3])
    # Kff <- compute_cov_1d(idx1 = x_vec, tau = theta[2], h = theta[3])
    Kdf <- computeCovDer1(idx1 = der_vec, idx2 = x_vec, 
                            tau = theta[2], h = theta[3]) 
    Kdd <- computeCovDer2(idx1 = der_vec, tau = theta[2], h = theta[3])
    n <- length(x_vec)
    Sigma <- Kff - emulator::quad.form.inv(Kdd, Kdf)
    # Psi <- diag(theta[1] ^ 2, n) + Sigma
    Psi <- diag(theta[1] ^ 2, n) + Sigma
    Psi <- (Psi + t(Psi)) / 2
    if (!is.positive.definite(Psi)) {
        Psi <- as.matrix(nearPD(Psi)$mat)
    }
    # if (is.sig.par) {
    #     return(-mvnfast::dmvn(y, mu = rep(0L, n),
    #                           sigma = diag(theta[1] ^ 2, n) + Sigma,
    #                           log = TRUE) -
    #                dinvgamma(theta[1] ^ 2, shape = ga_shape, rate = ga_rate,
    #                          log = TRUE))
    # } else {
    #     return(-mvnfast::dmvn(y, mu = rep(0L, n),
    #                           sigma = diag(theta[1] ^ 2, n) + Sigma,
    #                           log = TRUE))
    # }
    
    if (is.sig.par) {
        return(-mvnfast::dmvn(y, mu = rep(0L, n),
                              sigma = Psi,
                              log = TRUE) -
                   dinvgamma(theta[1] ^ 2, shape = ga_shape, rate = ga_rate,
                             log = TRUE))
    } else {
        return(-mvnfast::dmvn(y, mu = rep(0L, n),
                              sigma = Psi,
                              log = TRUE))
    }
}

log_mar_lik_gp_der_matern <- function(theta, y, x_vec, der_vec, H0, nu,
                               ga_shape = 6, ga_rate = 5,
                               is.sig.par = TRUE) {
    # theta = c(sig, tau, h) (no mu)
    # if (is.null(H0)) {
    #     H0 <- outer(as.vector(x_vec), as.vector(x_vec), 
    #                 FUN = function(x1, x2) (x1 - x2))
    # }
    Kff <- matern_ker(H0 = H0, nu = nu, tau = theta[2], l = theta[3])
    # Kff <- compute_cov_1d(idx1 = x_vec, tau = theta[2], h = theta[3])
    Kdf <- computeCovDer1(idx1 = der_vec, idx2 = x_vec, 
                          ker_der1_fcn = matern_ker_der1,
                          tau = theta[2], l = theta[3], nu = nu) 
    Kdd <- computeCovDer2(idx1 = der_vec, ker_der2_fcn = matern_ker_der2, 
                          tau = theta[2], l = theta[3], nu = nu)
    n <- length(x_vec)
    Sigma <- Kff - emulator::quad.form.inv(Kdd, Kdf)
    # Psi <- diag(theta[1] ^ 2, n) + Sigma
    if (is.sig.par) {
        return(-mvnfast::dmvn(y, mu = rep(0L, n), 
                              sigma = diag(theta[1] ^ 2, n) + Sigma, 
                              log = TRUE) -
                   dinvgamma(theta[1] ^ 2, shape = ga_shape, rate = ga_rate, 
                             log = TRUE))
    } else {
        return(-mvnfast::dmvn(y, mu = rep(0L, n), 
                              sigma = diag(theta[1] ^ 2, n) + Sigma, 
                              log = TRUE))
    }
}




