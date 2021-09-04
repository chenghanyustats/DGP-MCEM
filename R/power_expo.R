## original power_expo()


power_expo <- function(idx12_diff, tau, phi) {
    return(tau ^ 2 * exp(-phi / 2 * idx12_diff ^ 2))
}
