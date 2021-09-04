## original log_general_dbeta() 
## used in mcem_dgp()

log_dgbeta <- function(x, shape1, shape2, a, b) {
    if (shape1 == 1 && shape2 == 1) {
        if ((b - a) == 1) {
            return(0)
        } else {
            return(log(1 / (b - a)))
        }
    } else {
        return((shape1 - 1) * log(x - a) + (shape2 - 1) * log(b - x) - 
                   lbeta(shape1, shape2) - (shape1 + shape2 - 1) * log(b - a))
    }
}



