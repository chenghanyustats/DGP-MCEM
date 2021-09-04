
matern_ker_der1 <- function(x1, x2, H0, nu = 1.5, tau = 1, l = 1) {
    if (!(nu %in% c(0.5, 1.5, 2.5))) {
        stop("nu must be equal to 0.5, 1.5 or 2.5")
    }
    p <- nu - 0.5
    
    if (missing(x1) || missing(x2)) {
        if(!missing(H0)) {
            ## H0 can be a matrix from outer operation
            if (p == 0) {
                return(tau ^ 2 (-H0 / (abs(H0) * l)) * exp(- abs(H0) / l))
            } else if (p == 1) {
                b <- sqrt(3) * abs(H0) / l
                return(tau ^ 2 * (-3 / l ^ 2 * H0) * exp(-b))
            } else {
                b <- sqrt(5) * abs(H0) / l
                c <- -5 * H0 / (3 * l ^ 3 * abs(H0))
                g <- l * abs(H0) + sqrt(5) * H0 ^ 2
                return(tau ^ 2 * (c * g) * exp(-b))
            }
        } else {
            stop("must provide either (x1, x2) or H0")
        }
    } else {
        d <- abs(x1 - x2)
        if (p == 0) {
            return(tau ^ 2 * (-(x1 - x2) / (d * l)) * exp(- d / l))
        } else if (p == 1) {
            b <- sqrt(3) * d / l
            return(tau ^ 2 * (-3 / l ^ 2 * (x1 - x2)) * exp(-b))
        } else {
            b <- sqrt(5) * d / l
            c <- -5 * (x1 - x2) / (3 * l ^ 3 * d)
            g <- l * d + sqrt(5) * d ^ 2
            return(tau ^ 2 * (c * g) * exp(-b))
        } 
    }
}

