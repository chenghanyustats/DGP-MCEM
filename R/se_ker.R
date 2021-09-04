### squared exponential kernel

se_ker <- function(x1, x2, H0, tau = 1, h = 1) {
    if (missing(x1) || missing(x2)) {
        if(!missing(H0)) {
            ## H0 can be a matrix from outer operation
            return(tau ^ 2 * exp(- (H0 / h) ^ 2 / 2))
        } else {
            stop("must provide either (x1, x2) or H0")
        }
    } else {
        return(tau ^ 2 * exp(- ((x1 - x2) / h) ^ 2 / 2))
    }
}


se_ker_no_tau <- function(x1, x2, H0, h) {
    if (missing(x1) || missing(x2)) {
        if(!missing(H0)) {
            ## H0 can be a matrix from outer operation
            return(exp(- (H0 / h) ^ 2 / 2))
        } else {
            stop("must provide either (x1, x2) or H0")
        }
    } else {
        return(exp(- ((x1 - x2) / h) ^ 2 / 2))
    }
}


# se_ker <- function(x1, x2, H0, tau = 1, h = 1) {
#     if(!missing(H0)) {
#             ## H0 can be a matrix from outer operation
#             return(tau ^ 2 * exp(- (H0 / h) ^ 2 / 2))
#     } else {
#         return(tau ^ 2 * exp(- ((x1 - x2) / h) ^ 2 / 2))
#     }
# }
# 
# se_ker <- function(x1, x2, H0 = NULL, tau = 1, h = 1) {
#     if(!is.null(H0)) {
#         ## H0 can be a matrix from outer operation
#         A = (tau ^ 2 * exp(- (H0 / h) ^ 2 / 2))
#     } else {
#         A = (tau ^ 2 * exp(- ((x1 - x2) / h) ^ 2 / 2))
#     }
#     return(A)
# }
# 
# se_ker <- function(x1, x2, H0 = NULL, tau = 1, h = 1) {
#     return(tau ^ 2 * exp(- ((x1 - x2) / h) ^ 2 / 2))
# }
# 
# se_ker <- function(x1, x2, tau = 1, h = 1) {
#     return(tau ^ 2 * exp(- ((x1 - x2) / h) ^ 2 / 2))
# }



