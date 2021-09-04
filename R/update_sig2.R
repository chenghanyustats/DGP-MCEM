update_sig2 <- function(y, lambda, A, avgB, is.sig.par,
                        ga_shape = 5, ga_rate = 5) {
    # A: function of h
    # B: function of t and h
    n <- length(y)
    y_A_avgB_y <- emulator::quad.form(A + avgB, y)
    if (is.sig.par) {
        denom <- n / 2 + ga_shape + 1
        numer <- (1 / (2 * lambda)) * y_A_avgB_y + ga_rate
        return(numer / denom)
    } else {
        return(y_A_avgB_y / (lambda * (n + 2)))
    }

}









