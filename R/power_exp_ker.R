### power-exponential kernel

power_exp_ker <- function(x, y, tau = 1, h = 1, alpha = 1) {
    return(tau ^ 2 * exp(- (abs(x - y) / h) ^ alpha))
}