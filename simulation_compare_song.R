# ################################################################################
# # Fit Song(2006) nonparametric kernel smoothing
# ################################################################################
# Step 1: choose bandwidth using direct plug-in method by minimizing AMISE
library(KernSmooth)
bw_vec <- rep(0, 100)
for (i in 1:100) {
    bw_vec[i] <- dpill(x = YY[[i]]$x, y = YY[[i]]$y)
}

# Step 2: Obtain the fitted curve and fitted first-derivative curve
LQfit <- function(x, y, h) {
    # x <- data[ , 1]
    # y <- data[ , 2]
    l <- length(x)
    beta0 <- rep(0, l)
    beta1 <- rep(0, l)
    for(i in 1:l) {
        x.reg <- x - x[i]
        w <- dnorm(x - x[i], 0, h)
        fit <- lm(y ~ x.reg + I(x.reg ^ 2), weights = w)
        beta0[i] <- fit$coe[1]
        beta1[i] <- fit$coe[2]
    }
    beta <- cbind(beta0, beta1)
    return(beta)
}



beta_lst <- list()
for (i in 1:100) {
    beta_lst[[i]] <- LQfit(x = YY[[i]]$x, y = YY[[i]]$y, h = bw_vec[i])
}

# Step 3: Add 95% confidence band on the fitted first-derivative curve and 
# obtain p-values
AddCI <- function(x, y, h, beta) {
    # x <- data[, 1]
    # y <- data[, 2]
    l <- length(x)
    beta0 <- beta[, 1]
    beta1 <- beta[, 2]
    ##Test for equilibruim points
    w <- rep(0, l)
    diag <- rep(0, l)
    upperCI <- rep(0, l)
    lowerCI <- rep(0, l)
    
    upperCI_f <- rep(0, l)
    lowerCI_f <- rep(0, l)
    
    se <- rep(0, l)
    se_f <- rep(0, l)
    
    z <- rep(0, l)
    p <- rep(0 ,l)
    options(object.size = 1000000000)
    
    ##Estimate sigma^2
    for(i in 1:l){
        Xi <- cbind(rep( 1, l), x - x[i], (x - x[i]) ^ 2)
        ker <- dnorm(Xi[, 2], 0, h)
        A <- matrix(0, ncol = 3, nrow = 3)
        A[1, 1] <- sum(ker)
        A[1, 2] <- ker %*% Xi[, 2]
        A[2, 1] <- A[1, 2]
        A[1, 3] <- ker %*% Xi[, 3]
        A[2, 2] <- A[1, 3]
        A[3, 1] <- A[1, 3]
        A[2, 3] <- ker %*% Xi[, 2] ^ 3
        A[3, 2] <- A[2, 3]
        A[3, 3] <- ker %*% Xi[, 3] ^ 2
        B <- solve(A)[1, ]
        C <- rbind(ker, ker * Xi[, 2], ker * Xi[,3])
        wi <- B %*% C
        diag[i] <- wi[i]
        w <- rbind(w, wi)
    }
    w <- w[2:( l + 1), ]
    second <- sum(w ^ 2)
    first <- 2 * sum(diag)
    v <- first - second
    vari <- 1 / ( l - v ) * sum((y - beta0) ^ 2)
    
    ##Calculate the 95% confidence band
    for(i in 1:l) {
        X <- cbind(rep(1, l), x - x[i], (x - x[i]) ^ 2)
        kernel <- dnorm(X[, 2], 0, h)
        An <- matrix(0, ncol = 3, nrow = 3)
        Bn <- matrix(0, ncol = 3, nrow = 3)
        An[1, 1] <- sum(kernel) / l
        An[1, 2] <- kernel %*% X[, 2] / l
        An[2, 1] <- An[1, 2]
        An[1, 3] <- kernel %*% X[, 3] / l
        An[2, 2] <- An[1, 3]
        An[3, 1] <- An[1, 3]
        An[2, 3] <- kernel %*% X[, 2] ^ 3 / l
        An[3, 2] <- An[2, 3]
        An[3, 3] <- kernel %*% X[, 3] ^ 2 / l
        kernel2 <- kernel ^ 2
        Bn[1, 1] <- sum(kernel2) / l / l
        Bn[1, 2] <- kernel2 %*% X[, 2] / l / l
        Bn[2, 1] <- Bn[1, 2]
        Bn[1, 3] <- kernel2 %*% X[, 3] / l / l
        Bn[2, 2] <- Bn[1, 3]
        Bn[3, 1] <- Bn[1, 3]
        Bn[2, 3] <- kernel2 %*% X[, 2] ^ 3 / l / l
        Bn[3, 2] <- Bn[2, 3]
        Bn[3, 3] <- kernel2 %*% X[, 3] ^ 2 / l / l
        sol <- solve(An)
        temp <- sol %*% Bn %*% sol
        temp2 <- temp[2, 2]
        temp1 <- temp[1, 1]
        se[i] <- sqrt(vari * temp2)
        se_f[i] <- sqrt(vari * temp1)
        
        z[i] <- abs(beta1[i] / se[i])
        p[i] <- (1 - pnorm(z[i])) * 2
        upperCI[i] <- beta1[i] + qnorm(0.975) * se[i]
        lowerCI[i] <- beta1[i] - qnorm(0.975) * se[i]
        
        upperCI_f[i] <- beta0[i] + qnorm(0.975) * se_f[i]
        lowerCI_f[i] <- beta0[i] - qnorm(0.975) * se_f[i]
    }
    upperCI <- round(upperCI, 5)
    lowerCI <- round(lowerCI, 5)
    
    upperCI_f <- round(upperCI_f, 5)
    lowerCI_f <- round(lowerCI_f, 5)
    
    
    p <- round(p, 5)
    CIp <- cbind(upperCI, lowerCI, p)
    return(CIp)
}

AddCI_f <- function(x, y, h, beta) {
  # x <- data[, 1]
  # y <- data[, 2]
  l <- length(x)
  beta0 <- beta[, 1]
  beta1 <- beta[, 2]
  ##Test for equilibruim points
  w <- rep(0, l)
  diag <- rep(0, l)
  upperCI <- rep(0, l)
  lowerCI <- rep(0, l)
  
  upperCI_f <- rep(0, l)
  lowerCI_f <- rep(0, l)
  
  se <- rep(0, l)
  se_f <- rep(0, l)
  
  z <- rep(0, l)
  p <- rep(0 ,l)
  options(object.size = 1000000000)
  
  ##Estimate sigma^2
  for(i in 1:l){
    Xi <- cbind(rep( 1, l), x - x[i], (x - x[i]) ^ 2)
    ker <- dnorm(Xi[, 2], 0, h)
    A <- matrix(0, ncol = 3, nrow = 3)
    A[1, 1] <- sum(ker)
    A[1, 2] <- ker %*% Xi[, 2]
    A[2, 1] <- A[1, 2]
    A[1, 3] <- ker %*% Xi[, 3]
    A[2, 2] <- A[1, 3]
    A[3, 1] <- A[1, 3]
    A[2, 3] <- ker %*% Xi[, 2] ^ 3
    A[3, 2] <- A[2, 3]
    A[3, 3] <- ker %*% Xi[, 3] ^ 2
    B <- solve(A)[1, ]
    C <- rbind(ker, ker * Xi[, 2], ker * Xi[,3])
    wi <- B %*% C
    diag[i] <- wi[i]
    w <- rbind(w, wi)
  }
  w <- w[2:( l + 1), ]
  second <- sum(w ^ 2)
  first <- 2 * sum(diag)
  v <- first - second
  vari <- 1 / ( l - v ) * sum((y - beta0) ^ 2)
  
  ##Calculate the 95% confidence band
  for(i in 1:l) {
    X <- cbind(rep(1, l), x - x[i], (x - x[i]) ^ 2)
    kernel <- dnorm(X[, 2], 0, h)
    An <- matrix(0, ncol = 3, nrow = 3)
    Bn <- matrix(0, ncol = 3, nrow = 3)
    An[1, 1] <- sum(kernel) / l
    An[1, 2] <- kernel %*% X[, 2] / l
    An[2, 1] <- An[1, 2]
    An[1, 3] <- kernel %*% X[, 3] / l
    An[2, 2] <- An[1, 3]
    An[3, 1] <- An[1, 3]
    An[2, 3] <- kernel %*% X[, 2] ^ 3 / l
    An[3, 2] <- An[2, 3]
    An[3, 3] <- kernel %*% X[, 3] ^ 2 / l
    kernel2 <- kernel ^ 2
    Bn[1, 1] <- sum(kernel2) / l / l
    Bn[1, 2] <- kernel2 %*% X[, 2] / l / l
    Bn[2, 1] <- Bn[1, 2]
    Bn[1, 3] <- kernel2 %*% X[, 3] / l / l
    Bn[2, 2] <- Bn[1, 3]
    Bn[3, 1] <- Bn[1, 3]
    Bn[2, 3] <- kernel2 %*% X[, 2] ^ 3 / l / l
    Bn[3, 2] <- Bn[2, 3]
    Bn[3, 3] <- kernel2 %*% X[, 3] ^ 2 / l / l
    sol <- solve(An)
    temp <- sol %*% Bn %*% sol
    temp2 <- temp[2, 2]
    temp1 <- temp[1, 1]
    se[i] <- sqrt(vari * temp2)
    se_f[i] <- sqrt(vari * temp1)
    
    z[i] <- abs(beta1[i] / se[i])
    p[i] <- (1 - pnorm(z[i])) * 2
    upperCI[i] <- beta1[i] + qnorm(0.975) * se[i]
    lowerCI[i] <- beta1[i] - qnorm(0.975) * se[i]
    
    upperCI_f[i] <- beta0[i] + qnorm(0.975) * se_f[i]
    lowerCI_f[i] <- beta0[i] - qnorm(0.975) * se_f[i]
  }
  upperCI <- round(upperCI, 5)
  lowerCI <- round(lowerCI, 5)
  
  upperCI_f <- round(upperCI_f, 5)
  lowerCI_f <- round(lowerCI_f, 5)
  
  
  p <- round(p, 5)
  CI_f <- cbind(upperCI_f, lowerCI_f, se_f)
  return(CI_f)
}


ci_song_lst <- list()
for (i in 1:100) {
    ci_song_lst[[i]] <- AddCI(x = YY[[i]]$x, y = YY[[i]]$y, h = bw_vec[i], beta = beta_lst[[i]])
}

ci_f_song_lst <- list()
for (i in 1:100) {
  ci_f_song_lst[[i]] <- AddCI_f(x = YY[[i]]$x, y = YY[[i]]$y, h = bw_vec[i], beta = beta_lst[[i]])
}

mean(sapply(ci_f_song_lst, function(x) x[, 1] - x[, 2]))

mean(sapply(ci_f_song_lst, function(x) x[, 3]))

# Step 4: Identify equilibrium points manually based on
# the results from Step 2 and Step 3. There are two scenarios:
all_info_lst <- list()
for (i in 1:100) {
    all_info_lst[[i]] <- cbind(pos = YY[[i]]$x, beta_lst[[i]], ci_song_lst[[i]])
}

song_der_zero_idx_lst <- list()
song_der_zero_raw_idx_lst <- list()
for(i in 1:100) {
    der_zero_idx <- which(diff(sign(all_info_lst[[i]][, 3]))!= 0)
    print(der_zero_idx)
    idx <- der_zero_idx %in% c(1, 2, 48, 49)
    raw_idx <- der_zero_idx %in% c(1, 49)
    song_der_zero_raw_idx_lst[[i]] <- der_zero_idx[!raw_idx]
    song_der_zero_idx_lst[[i]] <- der_zero_idx[!idx]
    # song_der_zero_pt[i, ] <- YY[[i]]$x[der_zero_idx]
    # if (length(der_zero_idx) == 0) {
    #     
    # } else {
    #     
    # }
}



lapply(song_der_zero_idx_lst, function(x) length(x))
table(unlist(lapply(song_der_zero_idx_lst, function(x) length(x))))



song_der_zero_pt_lst <- song_der_zero_idx_lst
# Step 5-1: Acquire approximate 95% confidence intervals for identified peak locations.
ci_pts_song_lst <- list()

for (i in 1:100) {
    stationary_pt <- all_info_lst[[i]][song_der_zero_pt_lst[[i]], ]
    if(length(song_der_zero_pt_lst[[i]]) == 1) {
        stationary_pt_upper <- stationary_pt[4]
        stationary_pt_lower <- stationary_pt[5]
    } else {
        stationary_pt_upper <- stationary_pt[, 4]
        stationary_pt_lower <- stationary_pt[, 5]
    }
    
    upper_vec <- c()
    lower_vec <- c()
    for(j in 1:length(song_der_zero_pt_lst[[i]])) {
        idx <- song_der_zero_pt_lst[[i]][j]
        if (all_info_lst[[i]][idx, 3] - 
            all_info_lst[[i]][idx-1, 3] < 0) {
            
            
            upper_idx_search <- which(diff(rev(all_info_lst[[i]][1:(idx-1), 3])) < 0)
            
            if (length(upper_idx_search) != 0) {
                upper_idx <- all_info_lst[[i]][(idx-1):(idx-1-upper_idx_search[1]), 3][all_info_lst[[i]][(idx-1):(idx-1-upper_idx_search[1]), 3] > 
                                                                                           stationary_pt_upper[j]][1]
            } else {
                upper_idx <- rev(all_info_lst[[i]][1:(idx-1), 3][all_info_lst[[i]][1:(idx-1), 3] > 
                                                                     stationary_pt_upper[j]])[1]
            }
            
            if (!is.na(upper_idx)) {
                lower_pt <- all_info_lst[[i]][which(all_info_lst[[i]][, 3] == upper_idx), 1]
            }
            
            
            
            lower_idx_search <-  which(diff(all_info_lst[[i]][(idx+1):n, 3]) > 0)
            
            
            if(length(lower_idx_search) != 0) {
                lower_idx <- all_info_lst[[i]][(idx+1):(idx+lower_idx_search[1]), 3][all_info_lst[[i]][(idx+1):(idx+lower_idx_search[1]), 3] < 
                                                                                         stationary_pt_lower[j]][1]
            } else{
                lower_idx <- all_info_lst[[i]][(idx+1):n, 3][all_info_lst[[i]][(idx+1):n, 3] < 
                                                                 stationary_pt_lower[j]][1]
            }
            
            
            if (!is.na(lower_idx)) {
                upper_pt <- all_info_lst[[i]][which(all_info_lst[[i]][, 3] == lower_idx), 1]
            }
            
            
            if(is.na(upper_idx) & !is.na(lower_idx)) {
                
                lower_pt <- 2 * stationary_pt[j, 1] - upper_pt
            }
            
            if(!is.na(upper_idx) & is.na(lower_idx)) {
                upper_pt <- 2 * stationary_pt[j, 1] - lower_pt
            }
            
            if(is.na(upper_idx) & is.na(lower_idx)) {
                lower_pt <- NA
                upper_pt <- NA
            }
            
            
            
            # lower_idx <- all_info_lst[[i]][(idx+1):n, 3][all_info_lst[[i]][(idx+1):n, 3] < 
            #                                                      stationary_pt_lower[j]][1]
            # upper_pt <- all_info_lst[[i]][which(all_info_lst[[i]][, 3] == lower_idx), 1]
            
            
            
        } else {
            upper_idx_search <- which(diff(all_info_lst[[i]][(idx+1):n, 3]) < 0)
            if (length(upper_idx_search) != 0) {
                upper_idx <- all_info_lst[[i]][(idx+1):(idx+upper_idx_search[1]), 3][all_info_lst[[i]][(idx+1):(idx+upper_idx_search[1]), 3] > 
                                                                                         stationary_pt_upper[j]][1]
            } else {
                upper_idx <- all_info_lst[[i]][(idx+1):n, 3][all_info_lst[[i]][(idx+1):n, 3] >
                                                                 stationary_pt_upper[j]][1]
                # upper_idx <- NA
            }
            
            if (!is.na(upper_idx)) {
                upper_pt <- all_info_lst[[i]][which(all_info_lst[[i]][, 3] == upper_idx), 1]
            }
            
            lower_idx_search <- which(diff(rev(all_info_lst[[i]][1:(idx-1), 3])) > 0)
            if(length(lower_idx_search) != 0) {
                lower_idx <- all_info_lst[[i]][(idx-1):(idx-lower_idx_search[1]), 3][all_info_lst[[i]][(idx-1):(idx-lower_idx_search[1]), 3] < 
                                                                                         stationary_pt_lower[j]][1]
            } else{
                # lower_idx <- NA
                lower_idx <- rev(all_info_lst[[i]][1:(idx-1), 3][all_info_lst[[i]][1:(idx-1), 3] <
                                                                     stationary_pt_lower[j]])[1]
            }
            
            # lower_idx <- rev(all_info_lst[[i]][1:(idx-1), 3][all_info_lst[[i]][1:(idx-1), 3] < 
            #                                                  stationary_pt_lower[j]])[1]
            if (!is.na(lower_idx)) {
                lower_pt <- all_info_lst[[i]][which(all_info_lst[[i]][, 3] == lower_idx), 1]
            }
            # lower_pt <- all_info_lst[[i]][which(all_info_lst[[i]][, 3] == lower_idx), 1]
            
            if(is.na(upper_idx) & !is.na(lower_idx)) {
                upper_pt <- 2 * stationary_pt[j, 1] - lower_pt
            }
            
            if(!is.na(upper_idx) & is.na(lower_idx)) {
                lower_pt <- 2 * stationary_pt[j, 1] - upper_pt
            }
            
            if(is.na(upper_idx) & is.na(lower_idx)) {
                lower_pt <- NA
                upper_pt <- NA
            }
        }
        upper_vec <- c(upper_vec, upper_pt)
        lower_vec <- c(lower_vec, lower_pt)
    }
    
    
    ci_pts_song_lst[[i]] <- cbind(lower_vec, upper_vec)
}

# # Step 5-2: Bonferroni correction and FDR
# p_bonferroni_lst <- list()
# p_bh_lst <- list()
# p_lst <- list()
# for (i in 1:100) {
#     p_lst[[i]] <- all_info_lst[[i]][song_der_zero_pt_lst[[i]], 6]
#     p_bonferroni_lst[[i]] <- p.adjust(p_lst[[i]], method = "bonferroni")
#     p_bh_lst[[i]] <- p.adjust(p_lst[[i]], method = "BH")
# }
# 
# sapply(p_lst, function(x) sum(x < 0.05))
# sapply(p_bonferroni_lst, function(x) sum(x < 0.05))
# sapply(p_bh_lst, function(x) sum(x < 0.05))


# Step 6: plotting
imp_lst <- list()

for(i in 1:100) {
    position <- song_der_zero_idx_lst[[i]]
    x_val <- YY[[i]]$x[position]
    lCI <- ci_pts_song_lst[[i]][, 1]
    uCI <- ci_pts_song_lst[[i]][, 2]
    imp_lst[[i]] <- list(position = position,
                         x = x_val,
                         lCI = lCI,
                         uCI = uCI)
}



Curve <- function(x, y, beta, CIp, imp) {
    # x <- data[ , 1]
    # y <- data[ , 2]
    l <- length(x)
    ma <- max(x)
    beta0 <- beta[ , 1]
    beta1 <- beta[ , 2]
    upperCI <- CIp[ , 1]
    lowerCI <- CIp[ , 2]
    # par(mfrow = c( 2, 1))
    ## nonparametric fitted curve with 95% CI
    plot(x, y, xlim = c(0, ma), type = 'p', ylim = c(0, 2.5),
         xlab = 'x', ylab = 'y', pch = 19)
    lines(x, beta0, lty = 4, lwd = 3, col = "green4")
    # title('Nonparametric Fitted Curve')
    pos <-as.vector(imp$position)
    vx <-as.vector(imp$x)
    ll <- length(vx)
    vx1 <-as.vector(imp$lCI)
    vx2 <-as.vector(imp$uCI)
    vy <- beta0[pos]
    
    for (i in 1:ll){
        x1 <- vx1[i]
        x2 <- vx2[i]
        y1 <- vy[i]
        polygon(c(x1, x1, x2, x2), c(0, 0.04, 0.04, 0) + 0.02*i, col = i,
                border = "white")
        points(vx[i], 0.02*i, pch = 18, col = i, cex = 1.5)
    }

    # ##first-derivative curve with 95% CI
    # plot(x, beta1, type = 'l', xlim = c(0, ma), xlab = 'x', 
    #      ylab = 'First derivative')
    # title( 'Estimated First Derivative Curve with 95% Confidence Interval')
    # for (i in 1:( l - 1)){
    #     x1 <- x[i]
    #     x2 <- x[i + 1]
    #     y1 <- lowerCI[i]
    #     y2 <- lowerCI[i + 1]
    #     y3 <- upperCI[i]
    #     y4 <- upperCI[i + 1]
    #     polygon(c(x1, x1, x2, x2), c(y1, y3, y2, y4), col = "lightgrey")
    # }
    # cross <- rep(0, ll)
    # par(new = T)
    # plot(vx, cross, pch = "x", xlim = c( 0, ma), ylim = c( -2, 2), 
    #      xlab = '', ylab = '')
    # abline(0)
    # par(new = T)
    # plot(x, beta1, type = 'l', xlim = c( 0, ma), ylim = c( -2, 2),
    #      xlab = 'Chromosomal coordinates (in kb)', ylab = 'First derivative'
}

for(i in 1:100) {
    print(length(imp_lst[[i]]$lCI))
}


ci_no_na_lst <- list()

for(i in 1:100) {
    lCI <- ci_pts_song_lst[[i]][, 1]
    uCI <- ci_pts_song_lst[[i]][, 2]
    ci_no_na_lst[[i]] <- list(lCI = lCI[!is.na(lCI)],
                         uCI = uCI[!is.na(uCI)])
}



c_1_vec <- rep(0, 100)
c_2_vec <- rep(0, 100)
for (i in 1:100) {
    ## first s.p.
    lower_vec <- imp_lst[[i]]$lCI[!is.na(imp_lst[[i]]$lCI)]
    upper_vec <- imp_lst[[i]]$uCI[!is.na(imp_lst[[i]]$uCI)]
    c_1 <- 0
    for (j in 1:length(lower_vec)) {
        if (!(cri_pts[1] > lower_vec[j] && cri_pts[1] < upper_vec[j])) {
            c_1 <- c_1 + 1
        }
    }
    c_1_vec[i] <- c_1
    ## second s.p.
    c_2 <- 0
    for (j in 1:length(lower_vec)) {
        if (!(cri_pts[2] > lower_vec[j] && cri_pts[2] < upper_vec[j])) {
            c_2 <- c_2 + 1
        }
    }
    c_2_vec[i] <- c_2
}


# number of missed sp
sapply(1:100, function(x) length(ci_no_na_lst[[x]]$lCI))
sum(sapply(1:100, function(x) length(ci_no_na_lst[[x]]$lCI)) == c_1_vec)
sum(sapply(1:100, function(x) length(ci_no_na_lst[[x]]$lCI)) == c_2_vec)


# idxx <- 1
# x = YY[[idxx]]$x
# y = YY[[idxx]]$y
# beta = beta_lst[[idxx]]
# CIp = ci_song_lst[[idxx]] 
# imp = imp_lst[[idxx]]

Curve(x = YY[[idxx]]$x, y = YY[[idxx]]$y, beta = beta_lst[[idxx]], 
      CIp = ci_song_lst[[idxx]], 
      imp = imp_lst[[idxx]])
lines(idx, regfcn(idx),
      col = "red", lwd = 2)


par(mfrow = c(2, 5))
for (i in 1:100) {
    Curve(x = YY[[i]]$x, y = YY[[i]]$y, beta = beta_lst[[i]], 
          CIp = ci_song_lst[[i]], 
          imp = imp_lst[[i]])
    title(main = list(paste('NKS: data', i), cex = 1.4))
    lines(idx, regfcn(idx),
          col = "red", lwd = 2)
    abline(v = cri_pts, col = "purple")
}

for (k in 1:5) {
    plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
                     mu_test = pred_dgp_single_lst[[k]]$mu_test,
                     CI_Low_f = pred_dgp_single_lst[[k]]$ci_low,
                     CI_High_f = pred_dgp_single_lst[[k]]$ci_high,
                     xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
                     ylim = c(0, 2.5), is.true.fcn = TRUE, cex = 1, 
                     plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
                     pred_lwd = 2, title = paste("single-DGP data", k), 
                     is.legend = FALSE)
    # hist_t(sample_t_EM_single_lst[[k]], 0, 2, ylim = c(0, 4), col = "blue",
    #        den.line = FALSE, cri_pt = cri_pts)
    # title(list(paste("Distribution of t: data", k), cex = 1.8))
    for ( i in 1:hpd_interval_hist_lst[[k]]$no_cluster) {
        segments(x0 = hpd_interval_hist_lst[[k]]$ci_lower[[i]], 
                 y0 = 0.04 + 0.02*i, 
                 x1 = hpd_interval_hist_lst[[k]]$ci_upper[[i]], 
                 y1 = 0.04 + 0.02*i, col = i,
                 lwd = 4)
        points(mean(hpd_interval_hist_lst[[k]]$sample_cluster_lst[[i]]), 
               0.02*i, pch = 18, col = i, cex = 1.5)
    }
} 

# c_1_vec <- rep(0, 100)
# c_2_vec <- rep(0, 100)
# for (i in 1:100) {
#     ## first s.p.
#     c_1 <- 0
#     for (j in 1:hpd_interval_hist_lst[[i]]$no_cluster) {
#         if (!(cri_pts[1] > hpd_interval_hist_lst[[i]]$ci_lower[[j]] && 
#               cri_pts[1] < hpd_interval_hist_lst[[i]]$ci_upper[[j]])) {
#             c_1 <- c_1 + 1
#         }
#     }
#     c_1_vec[i] <- c_1
#     ## second s.p.
#     c_2 <- 0
#     for (j in 1:hpd_interval_hist_lst[[i]]$no_cluster) {
#         if (!(cri_pts[2] > (hpd_interval_hist_lst[[i]]$ci_lower[[j]]-0.009) && 
#               cri_pts[2] < (hpd_interval_hist_lst[[i]]$ci_upper[[j]]+0.009))) {
#             c_2 <- c_2 + 1
#         }
#     }
#     c_2_vec[i] <- c_2
# }
# 
# # number of missed sp
# sapply(1:100, function(x) hpd_interval_hist_lst[[x]]$no_cluster)
# sum(sapply(1:100, function(x) hpd_interval_hist_lst[[x]]$no_cluster) == c_1_vec)
# sum(sapply(1:100, function(x) hpd_interval_hist_lst[[x]]$no_cluster) == c_2_vec)
# 
# for ( i in 1:100) {
#     
# }
# 
# sapply(1:100, function(x) mean(unlist(hpd_interval_hist_lst[[x]]$ci_upper) - 
#                                    unlist(hpd_interval_hist_lst[[x]]$ci_lower)))
# sapply(1:100, function(x) unlist(hpd_interval_hist_lst[[x]]$ci_upper))



table(sapply(1:100, function(x) {
    hpd_interval_hist_lst[[x]]$no_cluster}))

table(sapply(song_der_zero_idx_lst, length))

table(sapply(song_der_zero_raw_idx_lst, length))


## RMSE

# rmse_song <- sapply(beta_lst, function(x){
#   pred_f <- x[, 1]
#   rmse_f(x, true_fcn_val)})
# mean(rmse_dgp_single)

# rmse_song <- rep(0, 100)
diff_sq_song <- matrix(0, 50, 100)
for (i in 1:100) {
  true_fcn_val <- regfcn(YY[[i]]$x)
  pred_f <- beta_lst[[i]][, 1]
  # rmse_song[i] <- rmse_f(pred_f, true_fcn_val)
  diff_sq_song[, i] <- (pred_f - true_fcn_val) ^ 2
}
# mean(rmse_song)
rmse_song_x <- sqrt(apply(diff_sq_song, 1, mean))
plot(YY[[1]]$x, diff_sq_song[, 1], type = "l")



################################################################################


# hpd_interval_hist_lst[[1]]$ci_upper
# 
# table(sapply(hpd_interval_hist_lst, function(x)x$no_cluster))
post_mean_lst <- list()
for (k in 1:100) {
  post_mean_lst[[k]] <- sapply(hpd_interval_hist_lst[[k]]$sample_cluster_lst, function(x) mean(x))
}

crit_1_mean <- sapply(post_mean_lst, function(x) x[x < 1])
crit_2_mean <- sapply(post_mean_lst, function(x) x[x > 1])

mean( (unlist(crit_1_mean) - cri_pts[1]) ^ 2) 
mean( (unlist(crit_2_mean) - cri_pts[2]) ^ 2)



est_mean_lst_song <- lapply(imp_lst, function(x) x$x)


crit_1_mean_song <- sapply(est_mean_lst_song, function(x) x[x < 1])
crit_2_mean_song <- sapply(est_mean_lst_song, function(x) x[x > 1])

sqrt(mean( (unlist(crit_1_mean_song) - cri_pts[1]) ^ 2) )
sqrt(mean( (unlist(crit_2_mean_song) - cri_pts[2]) ^ 2))








