hist_t <- function(t_sample, a, b, cri_pt = NULL, ylim = NULL, main = "", 
                   is.freq = FALSE, den.line = TRUE, breaks = 30, 
                   is.ci = FALSE, col = "lightblue", xlab = "t",
                   cri_pt_col = NULL, cri_pt_lty = 1) {
    if(is.null(ylim)) {
        hist(t_sample, breaks = breaks, freq = is.freq, col = col, 
             border = FALSE, xlab = xlab, xlim = c(a, b), main = main)
    } else {
        hist(t_sample, breaks = breaks, freq = is.freq, col = col, 
             border = FALSE, xlab = xlab, ylim = ylim, xlim = c(a, b),
             main = main)
    }
    if (den.line) {
        lines(density(t_sample), col = "red", lwd = 4) 
    }
    if(is.ci) {
        abline(v = quantile(t_sample, probs = c(0.025, 0.975)), 
               col = "green", lwd = 3)
    }
    if(!is.null(cri_pt)) {
        abline(v = cri_pt, col = cri_pt_col, lwd = 1, lty = 2)
    }
}