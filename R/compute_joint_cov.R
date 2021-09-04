# original compute_Kjoint()
# depends on compute_cov_1d(), computeCovDer1_h(), computeCovDer2_h()


# compute_joint_cov <- function(idx_obs, idx_der, sig, tau, phi, w = 0) {
#     Kff <- compute_cov_1d(idx1 = idx_obs, tau = tau, phi = phi)
#     Kdf <- computeCovDer1_h(idx1 = idx_der, idx2 = idx_obs, tau,
#                             h = sqrt(1 / phi))
#     Kdd <- computeCovDer2_h(idx1 = idx_der, idx2 = idx_der, tau,
#                             h = sqrt(1 / phi))
#     # Noises
#     if (w != 0) {
#         K_joint <- rbind(cbind(Kff + diag(sig ^ 2, length(idx_obs)), t(Kdf)),
#                          cbind(Kdf, Kdd + diag(w ^ 2, length(idx_der))))
#     } else {
#         K_joint <- rbind(cbind(Kff + diag(sig ^ 2, length(idx_obs)), t(Kdf)),
#                          cbind(Kdf, Kdd))
#     }
#     return(K_joint)
# }



compute_joint_cov <- function(idx_obs, idx_der, sig, tau, phi = NULL, h, 
                              w = 0) {
    
    if (missing(h)) {
        if (is.null(phi)) {
            stop("must provide either phi or h")
        } else {
            h <- sqrt(1 / phi)
        }
    }
    if (!(is.null(phi) || missing(h))) {
        if(h != sqrt(1 / phi)) {
            stop("phi and h are not consistent")
        }
    }
    
    Kff <- compute_cov_1d(idx1 = idx_obs, tau = tau, h = h)
    Kdf <- computeCovDer1(idx1 = idx_der, idx2 = idx_obs, tau = tau, h = h) 
    Kdd <- computeCovDer2(idx1 = idx_der, idx2 = idx_der, tau = tau, h = h)
    
    # Noises 
    if (w != 0) {
        return(rbind(cbind(Kff + diag(sig ^ 2, length(idx_obs)), t(Kdf)), 
                         cbind(Kdf, Kdd + diag(w ^ 2, length(idx_der)))))
    } else {
        return(rbind(cbind(Kff + diag(sig ^ 2, length(idx_obs)), t(Kdf)), 
                         cbind(Kdf, Kdd)))
    }
    # return(K_joint)
}

compute_joint_cov_matern <- function(idx_obs, idx_der, sig, tau, l, nu = 1.5,
                              w = 0) {
    
    Kff <- compute_cov_1d(idx1 = idx_obs, ker_fcn = matern_ker, 
                          tau = tau, l = l, nu = nu)
    Kdf <- computeCovDer1(idx1 = idx_der, idx2 = idx_obs, 
                          ker_der1_fcn = matern_ker_der1, 
                          tau = tau, l = l, nu = nu) 
    Kdd <- computeCovDer2(idx1 = idx_der, idx2 = idx_der, 
                          ker_der2_fcn = matern_ker_der2, 
                          tau = tau, l = l, nu = nu)
    
    # Noises 
    if (w != 0) {
        return(rbind(cbind(Kff + diag(sig ^ 2, length(idx_obs)), t(Kdf)), 
                     cbind(Kdf, Kdd + diag(w ^ 2, length(idx_der)))))
    } else {
        return(rbind(cbind(Kff + diag(sig ^ 2, length(idx_obs)), t(Kdf)), 
                     cbind(Kdf, Kdd)))
    }
    # return(K_joint)
}

# system.time(
#     for (i in 1:1000) compute_Kjoint(idx_obs, idx_der, sig, tau, phi, w = 0)
# )

# function(idx1, idx2 = NULL, tau, phi = NULL, h = NULL) {
#     # idx1: 1d input vector
#     # idx2: 1d input vector
#     if (is.null(phi)) phi <- 1 / h ^ 2
#     if (is.null(idx2)) {
#         grid <- expand.grid(idx1, idx1)
#     } else {
#         grid <- expand.grid(idx1, idx2)
#     }
#     matrix(power_expo(grid[, 1] - grid[, 2], tau = tau, phi = phi),
#            nrow = length(idx1))
# }



# compute_joint_cov(idx_obs = 1:5, idx_der = 1:4, sig = 1, tau = 1, h = 1)







