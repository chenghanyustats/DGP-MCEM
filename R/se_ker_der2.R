se_ker_der2 <- function(x1, x2, H0, tau = 1, h = 1) {
    phi <- 1 / h ^ 2
    if (missing(x1) || missing(x2)) {
        if(!missing(H0)) {
            return(se_ker(H0, tau, h) * (phi - (phi * H0) ^ 2))
        } else {
            stop("must provide either (x1, x2) or H0 = x1 - x2")
        }
    } else {
        return(se_ker(x1, x2, tau = tau, h = h) * (phi - (phi * (x1 - x2)) ^ 2))
    }
}


se_ker_der2_no_tau <- function(x1, x2, H0, h = 1) {
    phi <- 1 / h ^ 2
    if (missing(x1) || missing(x2)) {
        if(!missing(H0)) {
            return(se_ker_no_tau(H0, h) * (phi - (phi * H0) ^ 2))
        } else {
            stop("must provide either (x1, x2) or H0 = x1 - x2")
        }
    } else {
        return(se_ker_no_tau(x1, x2, h = h) * (phi - (phi * (x1 - x2)) ^ 2))
    }
}
# 
# 
# se_ker_der2 <- function(x1, x2, H0, tau = 1, h = 1) {
#     phi <- 1 / h ^ 2
#     if(!missing(H0)) {
#         return(se_ker(H0, tau, h) * (phi - (phi * H0) ^ 2))
#     } else {
#         return(se_ker(x1, x2, tau, h) * (phi - (phi * (x1 - x2)) ^ 2))
#     }
# }
# 
# se_ker_der2 <- function(x1, x2, H0 = NULL, tau = 1, h = 1) {
#     phi <- 1 / h ^ 2
#     if(!is.null(H0)) {
#         A = (se_ker(H0, tau, h) * (phi - (phi * H0) ^ 2))
#     } else {
#         A = (se_ker(x1, x2, tau = tau, h = h) * (phi - (phi * (x1 - x2)) ^ 2))
#     }
#     return(A)
# }
# 
# 
# se_ker_der2 <- function(x1, x2, tau = 1, h = 1) {
#     phi <- 1 / h ^ 2
#     return(se_ker(x1, x2, tau, h) * (phi - (phi * (x1 - x2)) ^ 2))
# }


