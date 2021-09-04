# depends on power_expo_h()

power_expo_derivative1_h <- function(idx12_diff, tau, h) {
    return(power_expo_h(idx12_diff, tau, h) * (idx12_diff) / (h ^ 2))
}