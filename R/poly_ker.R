#####
# POLYNOMIAL KERNEL
#####
poly_ker <- function(x, y, tau = 1, d = 1) {
    return((tau ^ 2 + x * y) ^ d)
}