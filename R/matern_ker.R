#####
# matern kernel
#####
# se_ker <- function(x1, x2, H0, tau = 1, h = 1) {
#     if (missing(x1) || missing(x2)) {
#         if(!missing(H0)) {
#             ## H0 can be a matrix from outer operation
#             return(tau ^ 2 * exp(- (H0 / h) ^ 2 / 2))
#         } else {
#             stop("must provide either (x1, x2) or H0")
#         }
#     } else {
#         return(tau ^ 2 * exp(- ((x1 - x2) / h) ^ 2 / 2))
#     }
# }


matern_ker <- function(x1, x2, H0, nu = 1.5, tau = 1, l = 1) {
    ## H0 is difference matrix [xi - xj]
    if (!(nu %in% c(0.5, 1.5, 2.5))) {
        stop("nu must be equal to 0.5, 1.5 or 2.5")
    }
    p <- nu - 0.5
    
    if (missing(x1) || missing(x2)) {
        if(!missing(H0)) {
            ## H0 can be a matrix from outer operation
            if (p == 0) {
                return(tau ^ 2 * exp(- abs(H0) / l))
            } else if (p == 1) {
                b <- sqrt(3) * abs(H0) / l
                return(tau ^ 2 * (1 + b) * exp(-b))
            } else {
                b <- sqrt(5) * abs(H0) / l
                return(tau ^ 2 * (1 + b + b ^ 2 / 3) * exp(-b))
            }
        } else {
            stop("must provide either (x1, x2) or H0")
        }
    } else {
        d <- abs(x1 - x2)
        if (p == 0) {
            return(tau ^ 2 * exp(- d / l))
        } else if (p == 1) {
            b <- sqrt(3) * d / l
            return(tau ^ 2 * (1 + b) * exp(-b))
        } else {
            b <- sqrt(5) * d / l
            return(tau ^ 2 * (1 + b + b ^ 2 / 3) * exp(-b))
        } 
    }
}




