plot_subj_hpd_bar <- function(x_msec, hpdi_lst, r1_r2_mean, r1_r2_sd, 
                              title) {
    plot(rep(10, 11), seq(11), type = "l", xlab = "Time (msec)", yaxt='n',
         ylab = "Subject", cex.axis = 1, tcl = -0.4, 
         xlim = c(min(x_msec), max(x_msec)), 
         main = title)
    axis(2, at = 1:11, las = 2, tcl = -0.3, cex.axis = 0.9)
    for (k in 1:11) {
        n_clu <- length(hpdi_lst[[k]][[1]])
        for (j in 1:n_clu) {
            segments(x0 = hpdi_lst[[k]][[1]][j], 
                     y0 = k, 
                     x1 = hpdi_lst[[k]][[2]][j], 
                     y1 = k, col = "black", lwd = 4)
        }
    }
    abline(v = r1_r2_mean, col = "black", lwd = 2, lty = 2)
    
    r1_intvl <- r1_r2_mean[1] + c(-1, 1) * qnorm(0.975) * r1_r2_sd[1]
    r2_intvl <- r1_r2_mean[2] + c(-1, 1) * qnorm(0.975) * r1_r2_sd[2]
    
    polygon(c(seq(r1_intvl[1], r1_intvl[2], length = 100), 
              rev(seq(r1_intvl[1], r1_intvl[2], length = 100))), 
            c(rep(21, 100), rep(0, 100)), 
            col = rgb(0, 0, 0, 0.1), border = NA)
    polygon(c(seq(r2_intvl[1], r2_intvl[2], length = 100), 
              rev(seq(r2_intvl[1], r2_intvl[2], length = 100))), 
            c(rep(21, 100), rep(0, 100)), 
            col = rgb(0, 0, 0, 0.1), border = NA)
}