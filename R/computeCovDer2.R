# depends on se_ker_der2()


computeCovDer2_h <- function(idx1, idx2 = NULL, tau, h) {
    if(is.null(idx2)) {
        grid <- expand.grid(idx1, idx1)
    } else {
        grid <- expand.grid(idx1, idx2)
    }
    K22 <- matrix(power_expo_derivative2_h(idx12_diff = grid[[1]] - grid[[2]],
                                           tau = tau, h = h),
                  nrow = length(idx1))
    return(K22)
}

computeCovDer2 <- function(idx1, idx2, ker_der2_fcn = se_ker_der2, ...) {
    if(missing(idx2)) {
        return(outer(idx1, idx1, function(x1, x2) ker_der2_fcn(x1, x2, ...)))
    } else {
        return(outer(idx1, idx2, function(x1, x2) ker_der2_fcn(x1, x2, ...)))
    }
}

# computeCovDer2_test <- function(idx1, idx2, tau, h) {
#     if(missing(idx2)) {
#         K22 <- outer(idx1, idx1, function(x1, x2) se_ker_der2(x1 = idx1, 
#                                                             x2 = idx1, 
#                                                             tau = tau, 
#                                                             h = h))
#     } else {
#         K22 <- outer(idx1, idx2, function(x1, x2) se_ker_der2(x1 = idx1, 
#                                                             x2 = idx2, 
#                                                             tau = tau, 
#                                                             h = h))
#     }
#     return(K22)
# }

# microbenchmark(computeCovDer2_h(idx1 = 1:5, idx2 = 1:4, tau = 1, h = 1),
#                computeCovDer2(1:5, 1:4, tau = 1, h = 1), times = 1000)




# computeCovDer2_h(idx1 = 1:5, 4:7, tau = 1, h = 1)
# computeCovDer2(idx1 = 1:5, 4:7,tau = 1, h = 1)
# computeCovDer2_test(idx1 = 1:5, 4:7, tau = 1, h = 1)






