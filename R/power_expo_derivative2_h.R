# depends on power_expo_h()


power_expo_derivative2_h <- function(idx12_diff, tau, h) {
    phi <- 1 / h ^ 2
    k22 <- power_expo_h(idx12_diff, tau, h) * (phi - (phi * idx12_diff) ^ 2)
    return(k22)
}