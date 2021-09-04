### original computeCov_new() function 
### depends on se_ker()

# compute_cov_1d <- function(idx1, idx2 = NULL, tau, phi = NULL, h = NULL) {
#     # idx1: 1d input vector
#     # idx2: 1d input vector
#     if (is.null(phi)) phi <- 1 / h ^ 2
#     if (is.null(idx2)) {
#         grid <- expand.grid(idx1, idx1)
#     } else {
#         grid <- expand.grid(idx1, idx2)
#     }
#     return(matrix(power_expo(grid[[1]]-grid[[2]], tau = tau, phi = phi),
#                   nrow = length(idx1)))
# }



compute_cov_1d <- function(idx1, idx2, H0, ker_fcn = se_ker, ...) {
    # if (missing(idx1) || missing(idx2)) {
    #     if(!missing(H0)) {
    #         return(tau ^ 2 * exp(- (H0 / h) ^ 2 / 2))
    #     } else {
    #         stop("must provide either (x1, x2) or H0 = x1 - x2")
    #     }
    # } else {
    #     return(tau ^ 2 * exp(- ((x1 - x2) / h) ^ 2 / 2))
    # }
    # 
    
    if(missing(idx2)) {
        return(outer(idx1, idx1, function(x1, x2) ker_fcn(x1, x2, ...)))
    } else {
        return(outer(idx1, idx2, function(x1, x2) ker_fcn(x1, x2, ...)))
    }
}


# compute_cov_1d(1:5, 1:5)
# compute_cov_1d(idx1 = 1:5, idx2 = 1:3, tau = 1, h = 1)
