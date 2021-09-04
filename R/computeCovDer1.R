# depends on se_ker_der1()

computeCovDer1_h <- function(idx1, idx2, tau, h) {
    # take derivative w.r.t. the 1st argument
    grid <- expand.grid(idx1, idx2)
    K12 <- matrix(power_expo_derivative1_h(idx12_diff = grid[[1]] - grid[[2]],
                                           tau = tau, h = h),
                  nrow = length(idx1),
                  ncol = length(idx2))
    return(K12)
}


computeCovDer1 <- function(idx1, idx2, ker_der1_fcn = se_ker_der1, ...) {
    # take derivative w.r.t. the 1st argument
    return(outer(idx1, idx2, function(x1, x2) ker_der1_fcn(x1, x2, ...)))
}


# computeCovDer1_test <- function(idx1, idx2, tau, h) {
#     # take derivative w.r.t. the 1st argument
#     K12 <- outer(idx1, idx2, function(x1, x2) se_ker_der1(x1 = idx1, 
#                                                           x2 = idx2, 
#                                                           tau = tau,
#                                                           h = h))
#     return(K12)
# }
# 
# computeCovDer1_h(idx1 = 1:5, idx2 = 1:4, tau = 1, h = 1)
# 
# computeCovDer1(1:5, idx2 = 1:4, tau = 1, h = 1)
# 
# computeCovDer1_test(1:5, idx2 = 1:4, tau = 1, h = 1)