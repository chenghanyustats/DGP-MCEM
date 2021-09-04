## code to prepare `sim_compare` dataset goes here

library(Deriv)
## Seeds to generate simulated data sets
seeds <- c(8710, 6422, 5403, 1074, 2114, 5065, 8086, 4862, 2898, 8941, 
           7328, 5974, 8922, 7589, 1699, 7910, 7761, 8942, 4659, 8816,
           9300, 7719, 8928, 6094, 6857, 3086, 1394, 3497, 1613, 1777,
           8270, 8660, 7928, 2784, 1661, 8146, 7374, 3052, 1891, 1560, 
           2800, 5382, 6262, 5093, 1653, 9808, 5867, 3508, 9801, 2534,
           4817, 7236, 9046, 5950, 6724, 2074, 5409, 5741, 3067, 2866, 
           1618, 5228, 3307, 7845, 4130, 4770, 5720, 1721, 2387, 5442, 
           8119, 2649, 9277, 1114, 8876, 3791, 5149, 4653, 8712, 9848,
           8788, 5027, 2761, 4887, 7317, 6457, 3353, 8482, 9962, 7971,
           1912, 1971, 9426, 1163, 3035, 7683, 6387, 7087, 4704, 2485)



## parameters
no_data <- 100
sig <- 1 / 4
n <- 50
x_a <- 0
x_b <- 2


## regression function, derivative function and stationary points
## Original amplitude is 0.5
## amplitude is 0.1
regfcn_tiny_amp <- function(x) {
    0.3 + 0.4 * x + 0.15 * sin(3.2 * x) + (1.1 / (1 + x ^ 2))
}


## amplitude is 0.2
regfcn_small_amp <- function(x) {
    0.3 + 0.4 * x + 0.2 * sin(3.2 * x) + (1.1 / (1 + x ^ 2))
}

## amplitude is 0.8
regfcn_large_amp <- function(x) {
    0.3 + 0.4 * x + 0.8 * sin(3.2 * x) + (1.1 / (1 + x ^ 2))
}

regfcn_der_tiny_amp <- Deriv::Deriv(regfcn_tiny_amp)
regfcn_der_small_amp <- Deriv::Deriv(regfcn_small_amp)
regfcn_der_large_amp <- Deriv::Deriv(regfcn_large_amp)

cri_pts_tiny_amp <- c(stats::uniroot(regfcn_der_tiny_amp, c(x_a, 1))$root, 
                       stats::uniroot(regfcn_der_tiny_amp, c(1, x_b))$root)
cri_pts_small_amp <- c(stats::uniroot(regfcn_der_small_amp, c(x_a, 1))$root, 
             stats::uniroot(regfcn_der_small_amp, c(1, x_b))$root)
cri_pts_large_amp <- c(stats::uniroot(regfcn_der_large_amp, c(x_a, 1))$root, 
             stats::uniroot(regfcn_der_large_amp, c(1, x_b))$root)

## create 100 simulated data sets
YY_tiny_amp <- list()
YY_small_amp <- list()
YY_large_amp <- list()

for (k in 1:no_data) {
    x <- sort(runif(n, x_a, x_b))
    x <- sort(x)
    set.seed(seeds[k])
    y_tiny_amp <- regfcn_tiny_amp(x) + rnorm(n, 0, sig)
    y_small_amp <- regfcn_small_amp(x) + rnorm(n, 0, sig)
    y_large_amp <- regfcn_large_amp(x) + rnorm(n, 0, sig)
    YY_tiny_amp[[k]] <- list(x = x, y = y_small_amp)
    YY_small_amp[[k]] <- list(x = x, y = y_small_amp)
    YY_large_amp[[k]] <- list(x = x, y = y_large_amp)
}


## store data sets
sim_compare_data_amp <- list(YY_tiny_amp = YY_tiny_amp,
                             YY_small_amp = YY_small_amp, 
                             YY_large_amp = YY_large_amp, 
                             no_data = no_data, 
                             sig = sig, 
                             n = n, 
                             x_a = x_a, 
                             x_b = x_b, 
                             regfcn_tiny_amp = regfcn_tiny_amp, 
                             regfcn_der_tiny_amp = regfcn_der_tiny_amp,
                             cri_pts_tiny_amp = cri_pts_tiny_amp,
                             regfcn_small_amp = regfcn_small_amp, 
                             regfcn_der_small_amp = regfcn_der_small_amp,
                             cri_pts_small_amp = cri_pts_small_amp,
                             regfcn_large_amp = regfcn_large_amp, 
                             regfcn_der_large_amp = regfcn_der_large_amp,
                             cri_pts_large_amp = cri_pts_large_amp)


usethis::use_data(sim_compare_data_amp, compress = "xz")

idx <- seq(x_a, x_b, length.out = 200)
truefcn_tiny_amp <- regfcn_tiny_amp(idx)
par(mfrow = c(3, 3))
par(mar = c(4, 4, 2, 1))
for(k in 1:9) {
    plot(idx, truefcn_tiny_amp, col = "red", type = "l", xlim = c(x_a, x_b),
         ylim = c(0.3, 2.3), lwd = 2, main = k)
    points(YY_tiny_amp[[k]]$x, YY_tiny_amp[[k]]$y, pch = 3)
    points(cri_pts_tiny_amp, regfcn_tiny_amp(cri_pts_tiny_amp), 
           pch = "_", lwd = .1,
           col = "green", cex = 3)
}


idx <- seq(x_a, x_b, length.out = 200)
truefcn_small_amp <- regfcn_small_amp(idx)
par(mfrow = c(3, 3))
par(mar = c(4, 4, 2, 1))
for(k in 1:9) {
    plot(idx, truefcn_small_amp, col = "red", type = "l", xlim = c(x_a, x_b),
         ylim = c(0.3, 2.3), lwd = 2, main = k)
    points(YY_small_amp[[k]]$x, YY_small_amp[[k]]$y, pch = 3)
    points(cri_pts_small_amp, regfcn_small_amp(cri_pts_small_amp), 
           pch = "_", lwd = .1,
           col = "green", cex = 3)
}


idx <- seq(x_a, x_b, length.out = 200)
truefcn_large_amp <- regfcn_large_amp(idx)
par(mfrow = c(3, 3))
par(mar = c(4, 4, 2, 1))
for(k in 1:9) {
    plot(idx, truefcn_large_amp, col = "red", type = "l", xlim = c(x_a, x_b),
         ylim = c(0, 3), lwd = 2, main = k)
    points(YY_large_amp[[k]]$x, YY_large_amp[[k]]$y, pch = 3)
    points(cri_pts_large_amp, regfcn_large_amp(cri_pts_large_amp), 
           pch = "_", lwd = .1,
           col = "green", cex = 3)
}


# idx <- seq(0, 2, length.out = 200)
# truefcn_small_amp <- regfcn_small_amp(idx)
# truefcn_large_amp <- regfcn_large_amp(idx)
# 
# plot(idx, truefcn_large_amp, col = "red", type = "l", xlim = c(x_a, x_b),
#      ylim = c(0.2, 2.3), lwd = 2)
# lines(idx, truefcn_small_amp, col = "blue")


