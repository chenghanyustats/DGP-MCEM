se_ker_der1 <- function(x1, x2, H0, tau = 1, h = 1) {
    if (missing(x1) || missing(x2)) {
        if(!missing(H0)) {
            ### H0 can be a matrix
            return(se_ker(H0, tau, h) * (H0) / (h ^ 2))
        } else {
            stop("must provide either (x1, x2) or H0")
        }
    } else {
        return(se_ker(x1, x2, tau = tau, h = h) * (x1 - x2) / (h ^ 2))
    }
}

se_ker_der1_no_tau <- function(x1, x2, H0, h = 1) {
    if (missing(x1) || missing(x2)) {
        if(!missing(H0)) {
            ### H0 can be a matrix
            return(se_ker_no_tau(H0, h) * (H0) / (h ^ 2))
        } else {
            stop("must provide either (x1, x2) or H0")
        }
    } else {
        return(se_ker_no_tau(x1, x2, h = h) * (x1 - x2) / (h ^ 2))
    }
}

# se_ker_der1 <- function(x1, x2, H0, tau = 1, h = 1) {
#     if(!missing(H0)) {
#             ### H0 can be a matrix
#             return(se_ker(H0, tau, h) * (H0) / (h ^ 2))
#     } else {
#         return(se_ker(x1, x2, tau, h) * (x1 - x2) / (h ^ 2))
#     }
# }
# 
# se_ker_der1 <- function(x1, x2, H0 = NULL, tau = 1, h = 1) {
#     if(!is.null(H0)) {
#         ### H0 can be a matrix
#         A = (se_ker(H0, tau, h) * (H0) / (h ^ 2))
#     } else {
#         A = (se_ker(x1, x2, tau = tau, h = h) * (x1 - x2) / (h ^ 2))
#     }
#     return(A)
# }
# 
# se_ker_der1 <- function(x1, x2, H0 = NULL, tau = 1, h = 1) {
#     return(se_ker(x1, x2, tau = tau, h = h) * (x1 - x2) / (h ^ 2))
# }
# 
# se_ker_der1 <- function(x1, x2, tau = 1, h = 1) {
#     return(se_ker(x1, x2, tau, h) * (x1 - x2) / (h ^ 2))
# }


