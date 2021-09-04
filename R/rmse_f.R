rmse_f <- function(pred_f, true_f) {
    return(sqrt(mean((pred_f - true_f) ^ 2)))
}