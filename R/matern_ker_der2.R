
matern_ker_der2 <- function(x1, x2, H0, nu = 1.5, tau = 1, l = 1) {
    if (!(nu %in% c(0.5, 1.5, 2.5))) {
        stop("nu must be equal to 0.5, 1.5 or 2.5")
    }
    p <- nu - 0.5
    
    if (missing(x1) || missing(x2)) {
        if(!missing(H0)) {
            ## H0 can be a matrix from outer operation
            if (p == 0) {
                return(tau ^ 2 (-1 / l ^ 2) * exp(- abs(H0) / l))
            } else if (p == 1) {
                b <- sqrt(3) * abs(H0) / l
                c <-  3 * sqrt(3) * (- abs(H0)) / (l ^ 3) 
                return(tau ^ 2 * (3 / l ^ 2 + c) * exp(-b))
            } else {
                b <- sqrt(5) * abs(H0) / l
                c <- 5 * sqrt(5) * l * abs(H0) - 25 * H0 ^ 2 + 5 * l ^ 2
                return(tau ^ 2 * (1 / (3 * l ^ 4)) * c * exp(-b))
            }
        } else {
            stop("must provide either (x1, x2) or H0")
        }
    } else {
        d <- abs(x1 - x2)
        if (p == 0) {
            return(tau ^ 2 * (-1 / l ^ 2) * exp(- d / l))
        } else if (p == 1) {
            b <- sqrt(3) * d / l
            c <- 3 * sqrt(3) * (-d) / (l ^ 3)
            return(tau ^ 2 * (3 / l ^ 2 + c) * exp(-b))
        } else {
            b <- sqrt(5) * d / l
            c <- 5 * sqrt(5) * l * d - 25 * d ^ 2 + 5 * l ^ 2
            return(tau ^ 2 * (1 / (3 * l ^ 4)) * c * exp(-b))
        } 
    }
}




