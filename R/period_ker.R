#####
# PERIODIC KERNEL
#####
period_ker <- function(x, y, p = 1, tau = 1, h = 1) {
    return(tau ^ 2 * exp(-2 * (sin(pi * abs(x - y) / p) / h) ^ 2))
}
