detect_zero_der <- function(der_vec, x_test) {
    loc <- which(diff(sign(der_vec)) != 0)
    return(list(x_test[loc + 1]))
}



detect_zero_der <- function(der_vec, x_test) {
    loc <- which(diff(sign(der_vec)) != 0)
    loc_final <- rep(0, length(loc))
    for (i in 1:length(loc)) {
        loc_i <- loc[i]
        log_i_vec <- loc_i + c(-1, 0, 1)
        der_loc <- abs(der_vec[log_i_vec])
        loc_final[i] <- log_i_vec[which.min(der_loc)]
    }
    return(list(value = x_test[loc_final],
                loc = loc_final))
}

