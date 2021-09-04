###############################################################################
# functions for computing amplitude
###############################################################################
find_closest <- function(vec, number) {
    which(abs(vec - number) == min(abs(vec - number)))
}
get_amp_at_mean_lat <- function(pred_f, x_test, r1_sample, r2_sample) {
    ## fix r at posterior mean
    mm <- length(r1_sample)
    # amplitude_sample <- matrix(0, mm, 2)
    amplitude_sample_mean <- matrix(0, mm, 2)
    r1_post_mean <- mean(r1_sample)
    r2_post_mean <- mean(r2_sample)
    r1_mean_idx <- find_closest(x_test, r1_post_mean)
    r2_mean_idx <- find_closest(x_test, r2_post_mean)
    for (m in 1:mm) {
        # r1_idx <- find_closest(x_test, r1_sample[m])
        # r2_idx <- find_closest(x_test, r2_sample[m])
        # amplitude_sample[m, 1] <- abs(pred_f[r1_idx, m])
        # amplitude_sample[m, 2] <- abs(pred_f[r2_idx, m])
        amplitude_sample_mean[m, 1] <- abs(pred_f[r1_mean_idx, m])
        amplitude_sample_mean[m, 2] <- abs(pred_f[r2_mean_idx, m])
    }
    return(amplitude_sample_mean)
}

## max peak over range of latency r (Max Peak)
get_amp_max <- function(pred_f, x_test, lat_sample_1, lat_sample_2,
                        is.abs = FALSE, is.old = FALSE) {
    mm <- length(lat_sample_1)
    amplitude_sample <- matrix(0, mm, 2)
    
    a1 <- find_closest(x_test, range(lat_sample_1)[1])
    b1 <- find_closest(x_test, range(lat_sample_1)[2])
    for (m in 1:mm) {
        if(is.abs) {
            pred_r1 <- abs(pred_f[a1:b1, m])
            amplitude_sample[m, 1] <- max(pred_r1)
        } else {
            pred_r1 <- pred_f[a1:b1, m]
            amplitude_sample[m, 1] <- pred_r1[which.max(abs(pred_r1))]
        }
    }
    
    a2 <- find_closest(x_test, range(lat_sample_2)[1])
    b2 <- find_closest(x_test, range(lat_sample_2)[2])
    for (m in 1:mm) {
        if(is.abs) {
            pred_r2 <- abs(pred_f[a2:b2, m])
            amplitude_sample[m, 2] <- max(pred_r2)
        } else {
            pred_r2 <- pred_f[a2:b2, m]
            if(is.old) {
                amplitude_sample[m, 2] <- pred_r2[which.min(abs(pred_r2))]
            } else {
                amplitude_sample[m, 2] <- pred_r2[which.max(abs(pred_r2))]
            }
        }
    }
    return(amplitude_sample)
}



## peak at each draw of r (Predicted Peak)
get_amplitude <- function(pred_f, x_test, lat_sample_1, lat_sample_2,
                          is.abs = FALSE) {
    mm <- length(lat_sample_1)
    amplitude_sample <- matrix(0, mm, 2)
    # amplitude_sample_mean <- matrix(0, mm, 2)
    # r1_post_mean <- mean(r1_sample)
    # r2_post_mean <- mean(r2_sample)
    # r1_mean_idx <- find_closest(x_test, r1_post_mean)
    # r2_mean_idx <- find_closest(x_test, r2_post_mean)
    for (m in 1:mm) {
        r1_idx <- find_closest(x_test, lat_sample_1[m])
        r2_idx <- find_closest(x_test, lat_sample_2[m])
        if (is.abs) {
            amplitude_sample[m, 1] <- abs(pred_f[r1_idx, m])
            amplitude_sample[m, 2] <- abs(pred_f[r2_idx, m])
        } else {
            amplitude_sample[m, 1] <- pred_f[r1_idx, m]
            amplitude_sample[m, 2] <- pred_f[r2_idx, m]
        }
        
        # amplitude_sample_mean[m, 1] <- abs(pred_f[r1_mean_idx, m])
        # amplitude_sample_mean[m, 2] <- abs(pred_f[r2_mean_idx, m])
    }
    return(amplitude_sample)
}


## mean peak over range of latency r (Mean Peak)
get_amp_mean <- function(pred_f, x_test, lat_sample_1, lat_sample_2,
                         is.abs = FALSE) {
    mm <- length(lat_sample_1)
    amplitude_sample <- matrix(0, mm, 2)
    
    a1 <- find_closest(x_test, range(lat_sample_1)[1])
    b1 <- find_closest(x_test, range(lat_sample_1)[2])
    for (m in 1:mm) {
        if(is.abs) {
            pred_r1 <- abs(pred_f[a1:b1, m])
        } else {
            pred_r1 <- pred_f[a1:b1, m]
        }
        amplitude_sample[m, 1] <- mean(pred_r1)
    }
    
    a2 <- find_closest(x_test, range(lat_sample_2)[1])
    b2 <- find_closest(x_test, range(lat_sample_2)[2])
    for (m in 1:mm) {
        if(is.abs) {
            pred_r2 <- abs(pred_f[a2:b2, m])
        } else {
            pred_r2 <- pred_f[a2:b2, m]
        }
        amplitude_sample[m, 2] <- mean(pred_r2)
    }
    return(amplitude_sample)
}

## peak at mid latency r (Integral Peak)
get_amp_integrate <- function(pred_f, x_test, lat_sample_1, lat_sample_2,
                              is.abs = FALSE) {
    mm <- length(lat_sample_1)
    amplitude_sample <- matrix(0, mm, 2)
    
    a1 <- find_closest(x_test, range(lat_sample_1)[1])
    b1 <- find_closest(x_test, range(lat_sample_1)[2])
    for (m in 1:mm) {
        if(is.abs) {
            pred_r1 <- abs(pred_f[a1:b1, m])
        } else {
            pred_r1 <- pred_f[a1:b1, m]
        }
        
        all_sum_pred_1 <- sum(pred_r1)
        sum_pred_1 <- pred_r1[1]
        for (i in 2:length(pred_r1)) {
            sum_pred_1 <- sum_pred_1 + pred_r1[i]
            if (abs(sum_pred_1) > abs(all_sum_pred_1 / 2)) break
        }
        amplitude_sample[m, 1] <- (pred_r1[i-1] + pred_r1[i]) / 2
    }
    
    a2 <- find_closest(x_test, range(lat_sample_2)[1])
    b2 <- find_closest(x_test, range(lat_sample_2)[2])
    for (m in 1:mm) {
        if(is.abs) {
            pred_r2 <- abs(pred_f[a2:b2, m])
        } else {
            pred_r2 <- pred_f[a2:b2, m]
        }
        all_sum_pred_2 <- sum(pred_r2)
        sum_pred_2 <- pred_r2[1]
        for (j in 2:length(pred_r2)) {
            sum_pred_2 <- sum_pred_2 + pred_r2[j]
            if (abs(sum_pred_2) > abs(all_sum_pred_2 / 2)) break
        }
        amplitude_sample[m, 2] <- (pred_r2[j-1] + pred_r2[j]) / 2
    }
    return(amplitude_sample)
}



get_amp_int_t <- function(pred_f, x_test, t_sample,
                          is.abs = FALSE, is.detrend = FALSE) {
    mm <- length(t_sample)
    amplitude_sample <- matrix(0, mm, 2)
    if (is.detrend) {
        pred_f <- apply(pred_f, 2, function(x) x - mean(x))
    }
    
    a <- find_closest(x_test, range(t_sample)[1])
    b <- find_closest(x_test, range(t_sample)[2])
    for (m in 1:mm) {
        pred <- pred_f[a:b, m]
        
        pred_1 <- pred[pred < 0]
        pred_2 <- pred[pred > 0]
        
        all_sum_pred_1 <- sum(pred_1)
        sum_pred_1 <- pred_1[1]
        for (i in 2:length(pred_1)) {
            sum_pred_1 <- sum_pred_1 + pred_1[i]
            if (abs(sum_pred_1) > abs(all_sum_pred_1 / 2)) break
        }
        amplitude_sample[m, 1] <- (pred_1[i - 1] + pred_1[i]) / 2
        
        all_sum_pred_2 <- sum(pred_2)
        sum_pred_2 <- pred_2[1]
        for (i in 2:length(pred_2)) {
            sum_pred_2 <- sum_pred_2 + pred_2[i]
            if (abs(sum_pred_2) > abs(all_sum_pred_2 / 2)) break
        }
        amplitude_sample[m, 2] <- (pred_2[i - 1] + pred_2[i]) / 2
    }
    if (is.abs) amplitude_sample <- abs(amplitude_sample)
    return(amplitude_sample)
}