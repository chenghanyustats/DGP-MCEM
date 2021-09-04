## plotting sample paths of GP

gp_plot <- function(idx, obs, idx_test, predict_val,
                    ci_lower = NULL, ci_upper = NULL, idx_der, 
                    ylim = c(-8, 8), tau, h, sig, obs_col = "red4", 
                    main = NULL, is.multiplepath = FALSE, n_path = NULL,
                    is.ci = TRUE, is.title = TRUE, is.derline = FALSE) {
    
    plot(idx_test, predict_val[1, ], type = "l", col = "black", ylim = ylim, 
         lwd = 2, xlab = "x", ylab = "f(x)", cex.lab = 1.3)
    
    if(is.ci) {
        if (is.null(ci_lower) || is.null(ci_upper)) {
            stop("must provide ci_lower and ci_upper")
        }
        graphics::polygon(c(idx_test, rev(idx_test)),
                          c(ci_lower, rev(ci_upper)),
                          col = rgb(0, 0, 1, 0.2), border = NA)
    }
    
    if(is.multiplepath) {
        if (is.null(n_path)) {
            n_path <- nrow(predict_val)
        }
        for (i in 2:n_path) {
            lines(idx_test, predict_val[i, ], col = i, lwd = 2)
        }
    }
    
    points(idx, obs, pch = 3, lwd = 2, col = obs_col, cex = 1.2)
    
    if (is.title) {
        if(is.null(main)) {
            title(substitute(paste(tau, " = ", t1, "  ", h, " = ", p1, " ", 
                                   sigma, " = ", s1),
                             list(t1 = round(tau, 2), p1 = round(h, 2), 
                                  s1 = round(sig, 2))), cex.main = 1.8)   
        } else {
            title(main = main)
        }
    }
    
    if(is.derline) {
        abline(v = idx_der, lwd = 2, col = "purple", lty = 1)
    }
}






