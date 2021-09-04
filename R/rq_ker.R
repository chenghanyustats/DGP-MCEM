# Rational quadratic kernel

rq_ker <- function(x, y, alpha = 1, tau = 1, h = 1) {
    tau ^ 2 * (1 + (x - y) ^ 2 / (2 * alpha * h ^ 2)) ^ (-alpha)
}