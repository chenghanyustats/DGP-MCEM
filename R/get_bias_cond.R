get_bias_cond <- function(data, bias_cond = c(0, 1)) {
    data_lst <- list()
    for (i in 1:length(bias_cond)) {
        idx_bias <- which(data[1, 6, ] == bias_cond[i])
        data_i <- data[, , idx_bias]
        data_lst[[i]] <- data_i
    }
    return(data_lst)
}


