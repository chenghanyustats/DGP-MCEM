## original Stochastic_EM_multi()

Stochastic_EM_multi <- function(y, x, H0, theta_init = c(1, 1, 1), epsilon = 1e-4, 
                                D = 100, a_vec = c(0, 1), b_vec = c(1, 2), 
                                n_sample = 4000, max_iter = 100,
                                lower = c(0.001, 0.001, 0.001),
                                upper = c(1/0.001, 1/0.001, 1/0.001),
                                ctrl = list(TOL = 1e-5, trace = 0),
                                ga_shape = 5, ga_rate = 5,
                                is.sig.par = TRUE) {
    ################################
    # Arguments
    
    ################################
    thetas <- matrix(0, max_iter, length(theta_init))
    
    print("Stochastic EM (multiple t) Begins")
    
    eps <- epsilon + 1
    count <- 1
    # print(paste("epsilon =", round(eps, decimalplaces(epsilon) * 2), 
    #             " count", count))
    
    nn <- length(a_vec)  ## number of critical points
    
    # Initialization of parameters
    thetas[1, ] <- theta_k <- theta_init
    theta_new <- theta_k
    log_prior_den <-  log(prod(1 / (b_vec - a_vec)))
    
    while(eps > epsilon) {
        
        if (count == 100) {
            print("Does not converge")
            break
        }
        # ******** Stochastic E-step ************ #
        # Sampling posterior distributuion of t given theta = (sig, tau, h)
        
        if (count == 1) {
            # sample_t_mat <- matrix(c(sort(runif(D, min = a_vec[1], max = b_vec[1])), 
            #                      sort(runif(D, min = a_vec[2], max = b_vec[2]))), n_mc, 2)
            
            sample_t_mat <- matrix(0, D, nn)
            for (i in 1:nn) {
                # sample_t_mat[, i] <- sort(runif(n_mc, min = a_vec[i], max = b_vec[i]))
                sample_t_mat[, i] <- runif(D, min = a_vec[i], max = b_vec[i])
            }
        } else {
            opt_res_t_new_h <- solnp(pars = c(theta_k, (a_vec + b_vec) / 2), 
                                     fun = log_marginal_lik_gp_der_t_new_h,
                                     LB = c(lower, a_vec), UB = c(upper, b_vec),
                                     H0 = H0, control = ctrl, y = y, x_vec = x,
                                     ga_shape = ga_shape, ga_rate = ga_rate,
                                     is.sig.par = is.sig.par)
            M_const <- max(exp(-opt_res_t_new_h$values[length(opt_res_t_new_h$values)]))
            
            sample_t_mat <- matrix(0, nrow = D, ncol = nn)
            s <- 1
            while (s <= D) {
                # t1_star <- runif(1, min = a_vec[1], max = b_vec[1])
                # t2_star <- runif(1, min = a_vec[2], max = b_vec[2])
                t_star <- runif(nn, min = a_vec, max = b_vec)
                # candi_den <- 1 / (b_vec[1] - a_vec[1]) * 1 / (b_vec[2] - a_vec[2])
                candi_den <- prod(1 / (b_vec - a_vec))
                
                log_den <- -log_marginal_lik_gp_der_new_h(theta = theta_k, y = y,
                                                          x_vec = x, 
                                                          der_vec = t_star, 
                                                          H0 = H0,
                                                          ga_shape = ga_shape, 
                                                          ga_rate = ga_rate,
                                                          is.sig.par = is.sig.par)
                den <- exp(log_den + log_prior_den)
                if (runif(1) * (M_const / 10) < (den / candi_den)) {
                    sample_t_mat[s, ] <- t_star
                    cat("no. draw:", s, "\r")
                    s <- s + 1
                }
            }
        }
        
        # ******** M-step for theta and sigma ************ #
        res <- solnp(pars = theta_k, fun = marg_lik_gp_der_mc_new_h,
                     LB = lower, control = ctrl, y = y, x_vec = x, 
                     H0 = H0, der_mc_mat = sample_t_mat,
                     ga_shape = ga_shape, ga_rate = ga_rate,
                     is.sig.par = is.sig.par)
        theta_k <- res$par
        
        # ******** Update epsilon and parameters ************ #
        eps <- sum((theta_new - theta_k) ^ 2)
        theta_new <- theta_k
        count <- count + 1
        print(paste("epsilon =", round(eps, 4), " count", count))
        thetas[count, ] <- theta_new
    }
    
    sample_t_mat <- matrix(0, n_sample, nn)
    opt_res_t_new_h <- solnp(pars = c(theta_k, (a_vec + b_vec) / 2), 
                             fun = log_marginal_lik_gp_der_t_new_h,
                             LB = c(lower, a_vec), UB = c(upper, b_vec),
                             H0 = H0, control = ctrl, y = y, x_vec = x,
                             ga_shape = ga_shape, 
                             ga_rate = ga_rate,
                             is.sig.par = is.sig.par)
    M_const <- max(exp(-opt_res_t_new_h$values[length(opt_res_t_new_h$values)]))
    
    s <- 1
    while (s <= n_sample) {
        t_star <- runif(nn, min = a_vec, max = b_vec)
        candi_den <- prod(1 / (b_vec - a_vec))
        
        log_den <- -log_marginal_lik_gp_der_new_h(theta = theta_new, y = y,
                                                  x_vec = x, 
                                                  der_vec = t_star, 
                                                  H0 = H0,
                                                  ga_shape = ga_shape, 
                                                  ga_rate = ga_rate,
                                                  is.sig.par = is.sig.par)
        den <- exp(log_den + log_prior_den)
        if (runif(1) * (M_const / 10) < (den / candi_den)) {
            sample_t_mat[s, ] <- t_star
            cat("final draw:", s, "\r")
            s <- s + 1
        }
    }
    
    cat("\n")
    print("Done!")
    colnames(thetas) <- c("sigma", "tau", "h")
    
    return(list(sample_t_mat = sample_t_mat, thetas = thetas[1:count, ]))
}


stochastic_em_dgp_multi <- function(y, x, H0, theta_init = c(1, 1, 1), epsilon = 1e-4, 
                                D = 100, a_vec = c(0, 1), b_vec = c(1, 2), 
                                n_sample = 4000, max_iter = 100,
                                lower = c(0.001, 0.001, 0.001),
                                upper = c(1/0.001, 1/0.001, 1/0.001),
                                ctrl = list(TOL = 1e-5, trace = 0),
                                ga_shape = 5, ga_rate = 5,
                                is.sig.par = TRUE, tune = 1) {
    ################################
    # Arguments
    
    ################################
    thetas <- matrix(0, max_iter, length(theta_init))
    
    print("Stochastic EM (multiple t) Begins")
    
    eps <- epsilon + 1
    count <- 1
    # print(paste("epsilon =", round(eps, decimalplaces(epsilon) * 2), 
    #             " count", count))
    
    nn <- length(a_vec)  ## number of critical points
    
    # Initialization of parameters
    thetas[1, ] <- theta_k <- theta_init
    theta_new <- theta_k
    log_prior_den <-  log(prod(1 / (b_vec - a_vec)))
    
    while(eps > epsilon) {
        
        if (count == 100) {
            print("Does not converge")
            break
        }
        # ******** Stochastic E-step ************ #
        # Sampling posterior distributuion of t given theta = (sig, tau, h)
        
        if (count == 1) {
            
            
            # sample_t_mat <- matrix(0, D, nn)
            # for (i in 1:nn) {
            #     # sample_t_mat[, i] <- sort(runif(n_mc, min = a_vec[i], max = b_vec[i]))
            #     sample_t_mat[, i] <- runif(D, min = a_vec[i], max = b_vec[i])
            # }
            # 
            
            # sample_t_mat <- matrix(0, nn, D)
            # for (i in 1:nn) {
            #     sample_t_mat[i, ] <- runif(D, min = a_vec[i], max = b_vec[i])
            # }
            
            sample_t_mat <- matrix(stats::runif(D * nn, min = a_vec,
                                                max = b_vec), nn, D)
            # sample_t_mat <- matrix(stats::runif(D * nn, min = a_vec,
            #                                     max = b_vec), D, nn)
            # print("1")
        } else {
            # print("2")
            opt_res_t <- Rsolnp::solnp(pars = c(theta_k, (a_vec + b_vec) / 2), 
                                       fun = log_mar_lik_gp_der_t,
                                       LB = c(lower, a_vec), 
                                       UB = c(upper, b_vec),
                                       H0 = H0, control = ctrl, 
                                       y = y, x_vec = x,
                                       ga_shape = ga_shape, ga_rate = ga_rate,
                                       is.sig.par = is.sig.par)
            M_const <- max(exp(-opt_res_t$values[length(opt_res_t$values)]))
            # print("3")
            
            sample_t_mat <- matrix(0, nrow = nn, ncol = D)
            # sample_t_mat <- matrix(0, nrow = D, ncol = nn)
            s <- 1
            print("Sampling begins")
            iter <- 0
            while (s <= D) {
                t_star <- stats::runif(nn, min = a_vec, max = b_vec)
                candi_den <- 1 / prod(b_vec - a_vec)
                
                log_den <- -log_mar_lik_gp_der(theta = theta_k, y = y,
                                               x_vec = x, der_vec = t_star, 
                                               H0 = H0, ga_shape = ga_shape, 
                                               ga_rate = ga_rate,
                                               is.sig.par = is.sig.par)
                den <- exp(log_den + log_prior_den)
                # print(stats::runif(1) * (M_const / 100))
                # print(den / candi_den)
                
                if (stats::runif(1) * (M_const / tune) < (den / candi_den)) {
                    sample_t_mat[, s] <- t_star
                    cat("no. draw:", s, "\r")
                    s <- s + 1
                    iter <- 0
                }
                # iter <- iter + 1
                # if (iter %% 20 == 0) {
                #     tune <- 100 + tune
                # }
                # print(iter)
                # print(tune)
            }
            
        }
        # print("4")
        # ******** M-step for theta and sigma ************ #
        # marg_lik_gp_der_mc_new_h_1
        # log_mar_lik_gp_der_mc_multi
        res <- solnp(pars = theta_k, fun = log_mar_lik_gp_der_mc_multi,
                     LB = lower, control = ctrl, y = y, x_vec = x,
                     H0 = H0, der_mc_mat = sample_t_mat,
                     ga_shape = ga_shape, ga_rate = ga_rate,
                     is.sig.par = is.sig.par)
        theta_k <- res$par
        
        # ******** M-step for theta and sigma ************ #
        # res <- solnp(pars = theta_k, fun = marg_lik_gp_der_mc_new_h,
        #              LB = lower, control = ctrl, y = y, x_vec = x,
        #              H0 = H0, der_mc_mat = sample_t_mat,
        #              ga_shape = ga_shape, ga_rate = ga_rate,
        #              is.sig.par = is.sig.par)
        # theta_k <- res$par
        
        # res <- Rsolnp::solnp(pars = theta_k, fun = log_mar_lik_gp_der_mc,
        #                      LB = lower, control = ctrl, y = y, x_vec = x,
        #                      H0 = H0, der_mc_mat = sample_t_mat,
        #                      ga_shape = ga_shape, ga_rate = ga_rate,
        #                      is.sig.par = is.sig.par)
        # theta_k <- res$par
        # print("5")
        # ******** Update epsilon and parameters ************ #
        eps <- sum((theta_new - theta_k) ^ 2)
        theta_new <- theta_k
        count <- count + 1
        print(paste("epsilon =", round(eps, 4), " count", count))
        thetas[count, ] <- theta_new
    }
    

    opt_res_t <- Rsolnp::solnp(pars = c(theta_k, apply(sample_t_mat, 1, mean)), 
                               fun = log_mar_lik_gp_der_t,
                               LB = c(lower, a_vec), UB = c(upper, b_vec),
                               H0 = H0, control = ctrl, y = y, x_vec = x,
                               ga_shape = ga_shape, 
                               ga_rate = ga_rate,
                               is.sig.par = is.sig.par)
    M_const <- max(exp(-opt_res_t$values[length(opt_res_t$values)]))
    
    sample_t_mat <- matrix(0L, nrow = n_sample, ncol = nn)
    # sample_t_mat <- matrix(0L, nn, n_sample)
    s <- 1
    # tune <- tune + 100
    iter <- 0
    while (s <= n_sample) {
        t_star <- stats::runif(nn, min = a_vec, max = b_vec)
        candi_den <- 1 / prod(b_vec - a_vec)
        
        log_den <- -log_mar_lik_gp_der(theta = theta_new, y = y,
                                       x_vec = x, 
                                       der_vec = t_star, 
                                       H0 = H0,
                                       ga_shape = ga_shape, 
                                       ga_rate = ga_rate,
                                       is.sig.par = is.sig.par)
        den <- exp(log_den + log_prior_den)
        if (stats::runif(1) * (M_const / tune) < (den / candi_den)) {
            sample_t_mat[s, ] <- t_star
            cat("final draw:", s, "\r")
            s <- s + 1
            iter <- 0
        }
        # iter <- iter + 1
        # if (iter %% 20 == 0) {
        #     tune <- 100 + tune
        # }
    }
    
    cat("\n")
    print("Done!")
    colnames(thetas) <- c("sigma", "tau", "h")
    
    return(list(sample_t_mat = sample_t_mat, thetas = thetas[1:count, ],
                tune = tune))
}




stochastic_em_dgp_multi_matern <- function(y, x, H0, theta_init = c(1, 1, 1), epsilon = 1e-4, 
                                    D = 100, a_vec = c(0, 1), b_vec = c(1, 2), 
                                    n_sample = 4000, max_iter = 100,
                                    lower = c(0.001, 0.001, 0.001),
                                    upper = c(1/0.001, 1/0.001, 1/0.001),
                                    ctrl = list(TOL = 1e-5, trace = 0),
                                    ga_shape = 5, ga_rate = 5, nu = 1.5,
                                    is.sig.par = TRUE, tune = 1) {
    ################################
    # Arguments
    
    ################################
    thetas <- matrix(0, max_iter, length(theta_init))
    
    print("Stochastic EM (multiple t) Begins")
    
    eps <- epsilon + 1
    count <- 1
    # print(paste("epsilon =", round(eps, decimalplaces(epsilon) * 2), 
    #             " count", count))
    
    nn <- length(a_vec)  ## number of critical points
    
    # Initialization of parameters
    thetas[1, ] <- theta_k <- theta_init
    theta_new <- theta_k
    log_prior_den <-  log(prod(1 / (b_vec - a_vec)))
    
    while(eps > epsilon) {
        
        if (count == 100) {
            print("Does not converge")
            break
        }
        # ******** Stochastic E-step ************ #
        # Sampling posterior distributuion of t given theta = (sig, tau, h)
        
        if (count == 1) {
            
            
            # sample_t_mat <- matrix(0, D, nn)
            # for (i in 1:nn) {
            #     # sample_t_mat[, i] <- sort(runif(n_mc, min = a_vec[i], max = b_vec[i]))
            #     sample_t_mat[, i] <- runif(D, min = a_vec[i], max = b_vec[i])
            # }
            # 
            
            # sample_t_mat <- matrix(0, nn, D)
            # for (i in 1:nn) {
            #     sample_t_mat[i, ] <- runif(D, min = a_vec[i], max = b_vec[i])
            # }
            
            sample_t_mat <- matrix(stats::runif(D * nn, min = a_vec,
                                                max = b_vec), nn, D)
            # sample_t_mat <- matrix(stats::runif(D * nn, min = a_vec,
            #                                     max = b_vec), D, nn)
            # print("1")
        } else {
            # print("2")
            opt_res_t <- Rsolnp::solnp(pars = c(theta_k, (a_vec + b_vec) / 2), 
                                       fun = log_mar_lik_gp_der_t_matern,
                                       LB = c(lower, a_vec), 
                                       UB = c(upper, b_vec),
                                       H0 = H0, control = ctrl, 
                                       y = y, x_vec = x, nu = nu,
                                       ga_shape = ga_shape, ga_rate = ga_rate,
                                       is.sig.par = is.sig.par)
            M_const <- max(exp(-opt_res_t$values[length(opt_res_t$values)]))
            # print("3")
            
            sample_t_mat <- matrix(0, nrow = nn, ncol = D)
            # sample_t_mat <- matrix(0, nrow = D, ncol = nn)
            s <- 1
            print("Sampling begins")
            iter <- 0
            while (s <= D) {
                t_star <- stats::runif(nn, min = a_vec, max = b_vec)
                candi_den <- 1 / prod(b_vec - a_vec)
                
                log_den <- -log_mar_lik_gp_der_matern(theta = theta_k, y = y,
                                               x_vec = x, der_vec = t_star, 
                                               H0 = H0, ga_shape = ga_shape, 
                                               ga_rate = ga_rate, nu = nu,
                                               is.sig.par = is.sig.par)
                den <- exp(log_den + log_prior_den)
                # print(stats::runif(1) * (M_const / 100))
                # print(den / candi_den)
                
                if (stats::runif(1) * (M_const / tune) < (den / candi_den)) {
                    sample_t_mat[, s] <- t_star
                    cat("no. draw:", s, "\r")
                    s <- s + 1
                    iter <- 0
                }
                # iter <- iter + 1
                # if (iter %% 20 == 0) {
                #     tune <- 100 + tune
                # }
                # print(iter)
                # print(tune)
            }
            
        }
        # print("4")
        # ******** M-step for theta and sigma ************ #
        # marg_lik_gp_der_mc_new_h_1
        # log_mar_lik_gp_der_mc_multi
        res <- solnp(pars = theta_k, fun = log_mar_lik_gp_der_mc_multi_matern,
                     LB = lower, control = ctrl, y = y, x_vec = x,
                     H0 = H0, der_mc_mat = sample_t_mat, nu = nu,
                     ga_shape = ga_shape, ga_rate = ga_rate,
                     is.sig.par = is.sig.par)
        theta_k <- res$par
        
        # ******** M-step for theta and sigma ************ #
        # res <- solnp(pars = theta_k, fun = marg_lik_gp_der_mc_new_h,
        #              LB = lower, control = ctrl, y = y, x_vec = x,
        #              H0 = H0, der_mc_mat = sample_t_mat,
        #              ga_shape = ga_shape, ga_rate = ga_rate,
        #              is.sig.par = is.sig.par)
        # theta_k <- res$par
        
        # res <- Rsolnp::solnp(pars = theta_k, fun = log_mar_lik_gp_der_mc,
        #                      LB = lower, control = ctrl, y = y, x_vec = x,
        #                      H0 = H0, der_mc_mat = sample_t_mat,
        #                      ga_shape = ga_shape, ga_rate = ga_rate,
        #                      is.sig.par = is.sig.par)
        # theta_k <- res$par
        # print("5")
        # ******** Update epsilon and parameters ************ #
        eps <- sum((theta_new - theta_k) ^ 2)
        theta_new <- theta_k
        count <- count + 1
        print(paste("epsilon =", round(eps, 4), " count", count))
        thetas[count, ] <- theta_new
    }
    
    
    opt_res_t <- Rsolnp::solnp(pars = c(theta_k, apply(sample_t_mat, 1, mean)), 
                               fun = log_mar_lik_gp_der_t_matern,
                               LB = c(lower, a_vec), UB = c(upper, b_vec),
                               H0 = H0, control = ctrl, y = y, x_vec = x,
                               ga_shape = ga_shape, 
                               ga_rate = ga_rate, nu = nu,
                               is.sig.par = is.sig.par)
    M_const <- max(exp(-opt_res_t$values[length(opt_res_t$values)]))
    
    sample_t_mat <- matrix(0L, nrow = n_sample, ncol = nn)
    # sample_t_mat <- matrix(0L, nn, n_sample)
    s <- 1
    # tune <- tune + 100
    iter <- 0
    while (s <= n_sample) {
        t_star <- stats::runif(nn, min = a_vec, max = b_vec)
        candi_den <- 1 / prod(b_vec - a_vec)
        
        log_den <- -log_mar_lik_gp_der_matern(theta = theta_new, y = y,
                                       x_vec = x, 
                                       der_vec = t_star, 
                                       H0 = H0,
                                       ga_shape = ga_shape, 
                                       ga_rate = ga_rate, nu = nu, 
                                       is.sig.par = is.sig.par)
        den <- exp(log_den + log_prior_den)
        if (stats::runif(1) * (M_const / tune) < (den / candi_den)) {
            sample_t_mat[s, ] <- t_star
            cat("final draw:", s, "\r")
            s <- s + 1
            iter <- 0
        }
        # iter <- iter + 1
        # if (iter %% 20 == 0) {
        #     tune <- 100 + tune
        # }
    }
    
    cat("\n")
    print("Done!")
    colnames(thetas) <- c("sigma", "tau", "h")
    
    return(list(sample_t_mat = sample_t_mat, thetas = thetas[1:count, ],
                tune = tune))
}



# 
# log_mar_lik_gp_der_t
# log_mar_lik_gp_der_t_matern
# 
# 
# log_mar_lik_gp_der
# log_mar_lik_gp_der_matern
# 
# log_mar_lik_gp_der_mc_multi
# log_mar_lik_gp_der_mc_matern

# a_vec <- c(0, 1)
# b_vec <- c(1, 2)
# 
# t_mat <- matrix(0, 100, 2)
# 
# 
# runif(10, min = a_vec, max = b_vec)
# microbenchmark(matrix(runif(200, min = a_vec, max = b_vec), 100, 2),
#                for (i in 1:2) {
#                    t_mat[, i] <- runif(100, min = a_vec[i], max = b_vec[i])
#                })

# microbenchmark(prod(1 / (b_vec - a_vec)), 1/prod(b_vec - a_vec), times = 10000)
