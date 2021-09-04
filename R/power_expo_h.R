power_expo_h <- function(idx12_diff, tau, h) {
    return(tau ^ 2 * exp(-(idx12_diff / h) ^ 2 / 2))
}