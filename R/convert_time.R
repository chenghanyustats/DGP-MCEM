convert_time <- function(t_val, x_all) {
    t_val <- t_val * (max(x_all) - min(x_all)) + min(x_all)
    sec_to_msec(t_val, 200)
}