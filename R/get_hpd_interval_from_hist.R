get_hpd_interval_from_hist <- function(samples, breaks = 50, max_iter = 1000, 
                                       tol = 0.01,
                                       step_size = 0.001,
                                       target_prob = 0.95,
                                       is.plot = FALSE) {
    ## get hist 
    density_x <- hist(samples, breaks = breaks, freq = FALSE, plot = is.plot)
    ## get max density
    max_density <- max(density_x$density)
    ## decrease the density until 95% of sample are covered
    current_density <- max_density
    prob_value <- 0
    count <- 1
    sample_size <- length(samples)
    while(abs(prob_value - target_prob) > tol) {
        if (prob_value > target_prob) {
            current_density <- current_density + step_size
        } else {
            current_density <- current_density / sqrt(1 + count)
        }
        high_density_idx <- which(density_x$density > current_density)
        clu_idx <- which(diff(high_density_idx) > 2)
        (prob_value <- sum(density_x$counts[high_density_idx]) / sample_size)
        ci_lower_b <- c(density_x$breaks[high_density_idx][1], 
                        density_x$breaks[high_density_idx][clu_idx + 1])
        ci_upper_b <- c(density_x$breaks[high_density_idx][clu_idx],
                        density_x$breaks[high_density_idx][length(high_density_idx)])
        no_clu <- length(ci_lower_b)
        # samples[samples > ci_lower & samples < ci_upper]
        count <- count + 1
        # print(prob_value)
        if (count == max_iter) break
    }
    sample_clu_lst <- list()
    for (k in 1:no_clu) {
        sample_clu_lst[[k]] <- samples[samples > ci_lower_b[k] & 
                                           samples < ci_upper_b[k]]
    }
    ci_lower <- lapply(sample_clu_lst, min)
    ci_upper <- lapply(sample_clu_lst, max)
    return(list(no_cluster = no_clu, sample_cluster_lst = sample_clu_lst,
                ci_lower = ci_lower, ci_upper = ci_upper, 
                prob_value = prob_value,
                den_value = current_density))
}