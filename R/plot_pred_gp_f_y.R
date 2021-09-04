plot_pred_gp_f_y <- function(x, y, idx = NULL, x_test, mu_test, 
                             CI_Low_f, CI_High_f, 
                             CI_Low_y = NULL, CI_High_y = NULL, 
                             alpha = 2, is.pred.f = TRUE, 
                             is.pred.y = FALSE, xlab = "x", ylab = "y",
                             is.der.line = FALSE, cri_pts = NULL, 
                             col.poly = rgb(0, 0, 1, 0.2), pch = 19,
                             xlim = c(-4, 4), ylim = c(-2, 4), cex = 0.5,
                             title = NULL, is.legend = TRUE, pred_lwd = 1, 
                             y_lwd = 1, plot.type = "p", pred_lty = 2,
                             t_est = NULL, x_a = 0, x_b = 2, 
                             legend.loc = NULL, is.true.fcn = TRUE,
                             pred_col = "blue", true_fcn = regfcn,
                             reg_fcn_col = "red",
                             der_line_col = "darkgreen") {
    plot(x, y, pch = pch, xlim = xlim, ylim = ylim, cex = cex, xlab = xlab, 
         ylab = ylab, type = plot.type, lwd = y_lwd)
    if (is.true.fcn) {
        # lines(idx, regfcn(scale_back(idx, x_b - x_a, x_a)),
        #       col = "red", lwd = 2)
        lines(idx, true_fcn(idx),
              col = reg_fcn_col, lwd = 2.5)
    }
    # lines(idx, regfcn(idx), col = "red", lwd = 2)
    lines(x_test, mu_test, col = pred_col, lwd = pred_lwd, lty = pred_lty)
    if(is.pred.f) {
        polygon(c(x_test, rev(x_test)),
                c(CI_Low_f, rev(CI_High_f)),
                col = col.poly, border = NA)
    }
    if(is.pred.y) {
        polygon(c(x_test, rev(x_test)),
                c(CI_Low_y, rev(CI_High_y)),
                col = rgb(1, 0, 0, 0.4), border = NA)
    }
    if(is.der.line) {
        abline(v = cri_pts, lwd = 1.5, col = der_line_col, lty = 1)
    }
    if (is.legend) {
        if (is.pred.y) {
            if(is.der.line) {
                if(is.true.fcn) {
                    legend_text <- c("observed", "true f", 
                                     "mean predictive f* (y*)",
                                     "CI for f*", "CI for y*", "deriv. 0")
                    col <- c("black", "red", "blue", "purple", "pink", 
                             "dark green")
                    lwd = c(0, 2, 2, 0, 0, 2)
                    lty = c(0, 2, 2, 0, 0, 2)
                    pch = c(19, NA, NA, 22, 22, NA)
                    pt.bg = c(NA, NA, NA, "purple", "pink", NA)
                } else {
                    legend_text <- c("observed", "mean predictive f* (y*)",
                                     "CI for f*", "CI for y*", "deriv. 0")
                    col <- c("black", "blue", "purple", "pink", "dark green")
                    lwd = c(0, 2, 0, 0, 2)
                    lty = c(0, 1, 0, 0, 2)
                    pch = c(19, NA, 22, 22, NA)
                    pt.bg = c(NA, NA, "purple", "pink", NA)
                }
            } else {
                if(is.true.fcn) {
                    legend_text <- c("observed", "true f", 
                                     "mean predictive f* (y*)",
                                     "CI for f*", "CI for y*")
                    col <- c("black", "red", "blue", "purple", "pink")
                    lwd = c(0, 2, 2, 0, 0)
                    lty = c(0, 1, 1, 0, 0)
                    pch = c(19, NA, NA, 22, 22)
                    pt.bg = c(NA, NA, NA, "purple", "pink")
                } else {
                    legend_text <- c("observed", "mean predictive f* (y*)",
                                     "CI for f*", "CI for y*")
                    col <- c("black", "blue", "purple", "pink")
                    lwd = c(0, 2, 0, 0)
                    lty = c(0, 1, 0, 0)
                    pch = c(19, NA, 22, 22)
                    pt.bg = c(NA, NA, "purple", "pink")
                }
            }
        } else {
            if(is.der.line) {
                if(is.true.fcn) {
                    legend_text <- c("observed", "true f", "mean predictive f*",
                                     "CI for f*", "deriv. 0")
                    col <- c("black", "red", "blue", rgb(0, 0, 1, 0.2), 
                             "dark green")
                    lwd = c(0, 2, 2, 0, 2)
                    lty = c(0, 1, 1, 0, 2)
                    pch = c(19, NA, NA, 22, NA)
                    pt.bg = c(NA, NA, NA, rgb(0, 0, 1, 0.2), NA)
                } else {
                    legend_text <- c("observed", "mean predictive f*",
                                     "CI for f*", "deriv. 0")
                    col <- c("black", "blue", rgb(0, 0, 1, 0.2), "dark green")
                    lwd = c(0, 2, 0, 1)
                    lty = c(0, 1, 0, 2)
                    pch = c(19, NA, 22, NA)
                    pt.bg = c(NA, NA, rgb(0, 0, 1, 0.2), NA)
                }
            } else {
                if(is.true.fcn) {
                    legend_text <- c("observed", "true f", "mean predictive f*",
                                     "CI for f*")
                    col <- c("black", "red", "blue", rgb(0, 0, 1, 0.2))
                    lwd = c(0, 2, 2, 0)
                    lty = c(0, 1, 1, 0)
                    pch = c(19, NA, NA, 22)
                    pt.bg = c(NA, NA, NA, rgb(0, 0, 1, 0.2))
                } else {
                    legend_text <- c("observed", "mean predictive f*",
                                     "CI for f*")
                    col <- c("black", "blue", rgb(0, 0, 1, 0.2))
                    lwd = c(0, 2, 0)
                    lty = c(0, 1, 0)
                    pch = c(19, NA, 22)
                    pt.bg = c(NA, NA, rgb(0, 0, 1, 0.2))
                }
            }
        }
        if(!is.null(legend.loc)) {
            loc = legend.loc
        } else {
            loc = "bottomright"
        }
        legend(loc, legend_text, col = col, lwd = lwd, lty = lty, bty = "n",
               pch = pch, pt.bg = pt.bg, pt.cex = 1, text.font = 1, cex = 1)
    }
    if(!is.null(title)) { title(main = list(title, cex = 1.5)) }
    if(!is.null(t_est))  abline(v = t_est, col = "purple", lwd = 2)
}



plot_pred_gp_f_y <- function(x, y, idx = NULL, x_test, mu_test, CI_Low_f, CI_High_f, 
                             CI_Low_y = NULL, CI_High_y = NULL, 
                             alpha = 2, is.pred.f = TRUE, 
                             is.pred.y = FALSE, xlab = "x", ylab = "y",
                             is.der.line = FALSE, cri_pts = NULL, 
                             col.poly = rgb(0, 0, 1, 0.2), pch = 19,
                             xlim = c(-4, 4), ylim = c(-2, 4), cex = 0.5,
                             title = NULL, is.legend = TRUE, pred_lwd = 1, 
                             y_lwd = 1, plot.type = "p", pred_lty = 2,
                             t_est = NULL, x_a = 0, x_b = 2,
                             legend.loc = NULL, is.true.fcn = TRUE,
                             reg_fcn_col = "red", pred_col = "blue",
                             der_line_col = "darkgreen") {
    plot(x, y, pch = pch, xlim = xlim, ylim = ylim, cex = cex, xlab = xlab, 
         ylab = ylab, type = plot.type, lwd = y_lwd)
    if (is.true.fcn) {
        # lines(idx, regfcn(scale_back(idx, x_b - x_a, x_a)),
        #       col = "red", lwd = 2)
        lines(idx, regfcn(idx),
              col = reg_fcn_col, lwd = 2)
    }
    # lines(idx, regfcn(idx), col = "red", lwd = 2)
    lines(x_test, mu_test, col = pred_col, lwd = pred_lwd, lty = pred_lty)
    if(is.pred.f) {
        polygon(c(x_test, rev(x_test)),
                c(CI_Low_f, rev(CI_High_f)),
                col = col.poly, border = NA)
    }
    if(is.pred.y) {
        polygon(c(x_test, rev(x_test)),
                c(CI_Low_y, rev(CI_High_y)),
                col = rgb(1, 0, 0, 0.4), border = NA)
    }
    if(is.der.line) {
        abline(v = cri_pts, lwd = 1, col = der_line_col, lty = 1)
    }
    if (is.legend) {
        if (is.pred.y) {
            if(is.der.line) {
                if(is.true.fcn) {
                    legend_text <- c("observed", "true f", "mean predictive f* (y*)",
                                     "CI for f*", "CI for y*", "deriv. 0")
                    col <- c("black", "red", "blue", "purple", "pink", "dark green")
                    lwd = c(0, 2, 2, 0, 0, 2)
                    lty = c(0, 2, 2, 0, 0, 2)
                    pch = c(19, NA, NA, 22, 22, NA)
                    pt.bg = c(NA, NA, NA, "purple", "pink", NA)
                } else {
                    legend_text <- c("observed", "mean predictive f* (y*)",
                                     "CI for f*", "CI for y*", "deriv. 0")
                    col <- c("black", "blue", "purple", "pink", "dark green")
                    lwd = c(0, 2, 0, 0, 2)
                    lty = c(0, 1, 0, 0, 2)
                    pch = c(19, NA, 22, 22, NA)
                    pt.bg = c(NA, NA, "purple", "pink", NA)
                }
            } else {
                if(is.true.fcn) {
                    legend_text <- c("observed", "true f", "mean predictive f* (y*)",
                                     "CI for f*", "CI for y*")
                    col <- c("black", "red", "blue", "purple", "pink")
                    lwd = c(0, 2, 2, 0, 0)
                    lty = c(0, 1, 1, 0, 0)
                    pch = c(19, NA, NA, 22, 22)
                    pt.bg = c(NA, NA, NA, "purple", "pink")
                } else {
                    legend_text <- c("observed", "mean predictive f* (y*)",
                                     "CI for f*", "CI for y*")
                    col <- c("black", "blue", "purple", "pink")
                    lwd = c(0, 2, 0, 0)
                    lty = c(0, 1, 0, 0)
                    pch = c(19, NA, 22, 22)
                    pt.bg = c(NA, NA, "purple", "pink")
                }
            }
        } else {
            if(is.der.line) {
                if(is.true.fcn) {
                    legend_text <- c("observed", "true f", "mean predictive f*",
                                     "CI for f*", "deriv. 0")
                    col <- c("black", "red", "blue", rgb(0, 0, 1, 0.2), "dark green")
                    lwd = c(0, 2, 2, 0, 2)
                    lty = c(0, 1, 1, 0, 2)
                    pch = c(19, NA, NA, 22, NA)
                    pt.bg = c(NA, NA, NA, rgb(0, 0, 1, 0.2), NA)
                } else {
                    legend_text <- c("observed", "mean predictive f*",
                                     "CI for f*", "deriv. 0")
                    col <- c("black", "blue", rgb(0, 0, 1, 0.2), "dark green")
                    lwd = c(0, 2, 0, 1)
                    lty = c(0, 1, 0, 2)
                    pch = c(19, NA, 22, NA)
                    pt.bg = c(NA, NA, rgb(0, 0, 1, 0.2), NA)
                }
            } else {
                if(is.true.fcn) {
                    legend_text <- c("observed", "true f", "mean predictive f*",
                                     "CI for f*")
                    col <- c("black", "red", "blue", rgb(0, 0, 1, 0.2))
                    lwd = c(0, 2, 2, 0)
                    lty = c(0, 1, 1, 0)
                    pch = c(19, NA, NA, 22)
                    pt.bg = c(NA, NA, NA, rgb(0, 0, 1, 0.2))
                } else {
                    legend_text <- c("observed", "mean predictive f*",
                                     "CI for f*")
                    col <- c("black", "blue", rgb(0, 0, 1, 0.2))
                    lwd = c(0, 2, 0)
                    lty = c(0, 1, 0)
                    pch = c(19, NA, 22)
                    pt.bg = c(NA, NA, rgb(0, 0, 1, 0.2))
                }
            }
        }
        if(!is.null(legend.loc)) {
            loc = legend.loc
        } else {
            loc = "bottomright"
        }
        legend(loc, legend_text, col = col, lwd = lwd, lty = lty, bty = "n",
               pch = pch, pt.bg = pt.bg, pt.cex = 1, text.font = 1, cex = 1)
    }
    if(!is.null(title)) { title(main = list(title, cex = 1.5)) }
    if(!is.null(t_est))  abline(v = t_est, col = "purple", lwd = 2)
}