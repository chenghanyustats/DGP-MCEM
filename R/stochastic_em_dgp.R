# original Stochastic_EM()
# y = YY[[k]]$y
# x = YY[[k]]$x
# H0 = H0_diff[[k]]
# theta_init = c(1, 1, 1)
# epsilon = 1e-4
# D = 100
# a = 0
# b = 2
# n_sample = 4000
# max_iter = 100
# lower = c(0.001, 0.001, 0.001)
# upper = c(1/0.001, 1/0.001, 1/0.001)
# init_t = 0.5
# shape1 = 1
# shape2 = 1
# ctrl = list(TOL = 1e-5, trace = 0)
# ga_shape = 5
# ga_rate = 5
# is.sig.par = FALSE

stochastic_em_dgp <- function(y, x, H0, theta_init = c(1, 1, 1), epsilon = 1e-4, 
                          D = 100, a = 0, b = 2, n_sample = 4000, max_iter = 100,
                          lower = c(0.001, 0.001, 0.001),
                          upper = c(1/0.001, 1/0.001, 1/0.001), init_t = 0.5, 
                          shape1 = 1, shape2 = 1,
                          ctrl = list(TOL = 1e-5, trace = 0),
                          ga_shape = 5, ga_rate = 5,
                          is.sig.par = TRUE) {
    ################################
    # Arguments
    
    ################################
    thetas <- matrix(0L, max_iter, length(theta_init))
    
    print("Stochastic EM new (one t) Begins")
    
    eps <- epsilon + 1
    count <- 1
    
    # Initialization of parameters
    thetas[1, ] <- theta_k <- theta_init
    theta_new <- theta_k
    if (shape1 == 1 && shape2 == 1) {
        log_prior_den <- log_dgbeta((a + b) / 2, shape1 = shape1, 
                                    shape2 = shape2, a = a, b = b)
    }
    
    
    while(eps > epsilon) {
        
        if (count == 100) {
            print("Does not converge")
            break
        }
        
        # ******** Stochastic E-step ************ #
        # Sampling posterior distributuion of t given theta = (sig, tau, h)
        # print("1")
        if (count == 1) {
            sample_t <- matrix(stats::runif(D, min = a, max = b), 1, D)
        } else {
            sample_t <- matrix(0L, 1, D)
            # print("2")
            opt_res_t <- Rsolnp::solnp(pars = c(theta_k, init_t),
                                       fun = log_mar_lik_gp_der_t,
                                       LB = c(lower, a),
                                       UB = c(upper, b), H0 = H0, 
                                       control = ctrl, y = y, x_vec = x,
                                       ga_shape = ga_shape, ga_rate = ga_rate,
                                       is.sig.par = is.sig.par)
            
            M_const <- max(exp(-opt_res_t$values[length(opt_res_t$values)]))
            
            s <- 1
            print("Sampling begins")
            while (s <= D) {
                t1_star <- stats::runif(1, min = a, max = b)
                candi_den <- 1 / (b - a)
                
                log_den <- -log_mar_lik_gp_der(theta = theta_k, y = y, 
                                               x_vec = x, 
                                               der_vec = t1_star,
                                               H0 = H0,
                                               ga_shape = ga_shape, 
                                               ga_rate = ga_rate,
                                               is.sig.par = is.sig.par)   
                if (shape1 != 1 || shape2 != 1) {
                    log_prior_den <- log_dgbeta(t1_star, shape1 = shape1,
                                                shape2 = shape2, a = a, b = b)
                }
                
                den <- exp(log_den + log_prior_den)
                # (den / candi_den)
                if (stats::runif(1) * (M_const / 1) < (den / candi_den)) {
                    sample_t[[s]] <- t1_star
                    cat("no. draw:", s, "\r")
                    s <- s + 1
                }
            }
        }
        # print("3")
        # ******** M-step for theta and sigma ************ #
        res <- Rsolnp::solnp(pars = theta_k, fun = log_mar_lik_gp_der_mc,
                             LB = lower, control = ctrl, y = y, x_vec = x, 
                             H0 = H0, der_mc_mat = sample_t,
                             ga_shape = ga_shape, ga_rate = ga_rate,
                             is.sig.par = is.sig.par)
        theta_k <- res$par
        # print("4")
        # ******** Update epsilon and parameters ************ #
        eps <- sum((theta_new - theta_k) ^ 2)
        theta_new <- theta_k
        count <- count + 1
        print(paste("epsilon =", round(eps, 4), " count", count))
        thetas[count, ] <- theta_new
    }
    # print("5")
    opt_res_t <- Rsolnp::solnp(pars = c(theta_new, sum(sample_t) / D),
                               fun = log_mar_lik_gp_der_t,
                               LB = c(lower, a),
                               UB = c(upper, b), H0 = H0,
                               control = ctrl, y = y, x_vec = x,
                               ga_shape = ga_shape,
                               ga_rate = ga_rate,
                               is.sig.par = is.sig.par)
    
    # opt_res_t <- Rsolnp::solnp(pars = c(theta_new, init_t),
    #                            fun = log_mar_lik_gp_der_t,
    #                            LB = c(lower, a),
    #                            UB = c(upper, b), H0 = H0,
    #                            control = ctrl, y = y, x_vec = x,
    #                            ga_shape = ga_shape,
    #                            ga_rate = ga_rate,
    #                            is.sig.par = is.sig.par)
    
    M_const <- max(exp(-opt_res_t$values[length(opt_res_t$values)]))
    
    
    # sample_t <- matrix(0L, 1, D)
    sample_t <- matrix(0L, n_sample, 1)
    
    s <- 1
    while (s <= n_sample) {
        t1_star <- stats::runif(1, min = a, max = b)
        candi_den <- 1 / (b - a)
        
        log_den <- -log_mar_lik_gp_der(theta = theta_new, y = y, 
                                       x_vec = x, der_vec = t1_star,
                                       H0 = H0,
                                       ga_shape = ga_shape, 
                                       ga_rate = ga_rate,
                                       is.sig.par = is.sig.par)
        if (shape1 != 1 || shape2 != 1) {
            log_prior_den <-  log_dgbeta(t1_star, shape1 = shape1,
                                         shape2 = shape2, a = a, b = b)
        }
        
        den <- exp(log_den + log_prior_den)
        if (stats::runif(1) * (M_const / 1) < (den / candi_den)) {
            sample_t[s, ] <- t1_star
            cat("final draw:", s, "\r")
            s <- s + 1
        }
    }
    cat("\n")
    print("Done!")
    colnames(thetas) <- c("sigma", "tau", "h")
    
    return(list(sample_t = sample_t, thetas = thetas[1:count, ]))
}


stochastic_em_dgp <- function(y, x, H0, theta_init = c(1, 1, 1), epsilon = 1e-4, 
                              D = 100, a = 0, b = 2, n_sample = 4000, max_iter = 100,
                              lower = c(0.001, 0.001, 0.001),
                              upper = c(1/0.001, 1/0.001, 1/0.001), init_t = 0.5, 
                              shape1 = 1, shape2 = 1,
                              ctrl = list(TOL = 1e-5, trace = 0),
                              ga_shape = 5, ga_rate = 5,
                              a_h = 1, b_h = 1, mixture_prob = c(0.5, 0.5),
                              is.sig.par = TRUE,
                              is.h.par = FALSE, tune = 1) {
    ################################
    # Arguments
    
    ################################
    thetas <- matrix(0L, max_iter, length(theta_init))
    
    print("Stochastic EM new (one t) Begins")
    
    eps <- epsilon + 1
    count <- 1
    
    n_mixture <- length(shape1)
    
    
    # Initialization of parameters
    thetas[1, ] <- theta_k <- theta_init
    theta_new <- theta_k
    # if (shape1 == 1 && shape2 == 1) {
    #     log_prior_den <- log_dgbeta((a + b) / 2, shape1 = shape1, 
    #                                 shape2 = shape2, a = a, b = b)
    # }
    if (n_mixture == 1) {
        if (shape1 == 1 && shape2 == 1) {
            log_prior_den <- log_dgbeta((a + b) / 2, shape1 = shape1,
                                               shape2 = shape2, 
                                               a = a, b = b)
        }
    }
    
    while(eps > epsilon) {
        
        if (count == 100) {
            print("Does not converge")
            break
        }
        
        # ******** Stochastic E-step ************ #
        # Sampling posterior distributuion of t given theta = (sig, tau, h)
        # print("1")
        if (count == 1) {
            # if(n_mixture == 1) {
            #     sample_t <- rbeta(D, shape1 = shape1, shape2 = shape2) * (b - a) + a
            # } else {
            #     sample_t <- c()
            #     for (m in 1:n_mixture) {
            #         D_m <- ceiling(D * mixture_prob[m])
            #         sample_t_m <- rbeta(D_m, shape1 = shape1[m],
            #                             shape2 = shape2[m])
            #         sample_t <- c(sample_t, sample_t_m)
            #     }
            #     sample_t <- sample_t[1:D]
            # }
            # 
            # 

            sample_t <- matrix(stats::runif(D, min = a, max = b), 1, D)
        } else {
            sample_t <- matrix(0L, 1, D)
            # print("2")
            opt_res_t <- Rsolnp::solnp(pars = c(theta_k, init_t),
                                       fun = log_mar_lik_gp_der_t,
                                       LB = c(lower, a),
                                       UB = c(upper, b), H0 = H0, 
                                       control = ctrl, y = y, x_vec = x,
                                       ga_shape = ga_shape, ga_rate = ga_rate,
                                       is.sig.par = is.sig.par)
            
            M_const <- max(exp(-opt_res_t$values[length(opt_res_t$values)]))
            
            s <- 1
            print("Sampling begins")
            while (s <= D) {
                t1_star <- stats::runif(1, min = a, max = b)
                candi_den <- 1 / (b - a)
                
                log_den <- -log_mar_lik_gp_der(theta = theta_k, y = y, 
                                               x_vec = x, 
                                               der_vec = t1_star,
                                               H0 = H0,
                                               ga_shape = ga_shape, 
                                               ga_rate = ga_rate,
                                               is.sig.par = is.sig.par)   
                # if (shape1 != 1 || shape2 != 1) {
                #     log_prior_den <- log_dgbeta(t1_star, shape1 = shape1,
                #                                 shape2 = shape2, a = a, b = b)
                # }
                
                if (n_mixture == 1) {
                    if (shape1 != 1 || shape2 != 1) {
                        log_prior_den <- log_dgbeta(t1_star, shape1 = shape1,
                                                           shape2 = shape2,
                                                           a = a, b = b)
                    }
                } else {
                    prior_den <- 0
                    for (m in 1:n_mixture) {
                        prior_den <- prior_den +
                            mixture_prob[m] * exp(log_dgbeta(t1_star, shape1 = shape1[m],
                                                                    shape2 = shape2[m],
                                                                    a = a, b = b))
                    }
                    log_prior_den <- log(prior_den)
                }
                
                
                den <- exp(log_den + log_prior_den)
                # (den / candi_den)
                if (stats::runif(1) * (M_const / tune) < (den / candi_den)) {
                    sample_t[[s]] <- t1_star
                    cat("no. draw:", s, "\r")
                    s <- s + 1
                }
            }
        }
        # print("3")
        # ******** M-step for theta and sigma ************ #
        res <- Rsolnp::solnp(pars = theta_k, fun = log_mar_lik_gp_der_mc,
                             LB = lower, control = ctrl, y = y, x_vec = x, 
                             H0 = H0, der_mc_mat = sample_t,
                             ga_shape = ga_shape, ga_rate = ga_rate,
                             is.sig.par = is.sig.par,
                             is.h.par = is.h.par,
                             a_h = a_h, b_h = b_h)
        theta_k <- res$par
        mar_post_k <- res$values[length(res$values)]
        # print("4")
        # ******** Update epsilon and parameters ************ #
        eps <- sum((theta_new - theta_k) ^ 2)
        theta_new <- theta_k
        count <- count + 1
        print(paste("epsilon =", round(eps, 4), " count", count))
        print(paste("(increasing) marginal posterior =", -mar_post_k))
        thetas[count, ] <- theta_new
    }
    # print("5")
    opt_res_t <- Rsolnp::solnp(pars = c(theta_new, sum(sample_t) / D),
                               fun = log_mar_lik_gp_der_t,
                               LB = c(lower, a),
                               UB = c(upper, b), H0 = H0,
                               control = ctrl, y = y, x_vec = x,
                               ga_shape = ga_shape,
                               ga_rate = ga_rate,
                               is.sig.par = is.sig.par)
    
    # opt_res_t <- Rsolnp::solnp(pars = c(theta_new, init_t),
    #                            fun = log_mar_lik_gp_der_t,
    #                            LB = c(lower, a),
    #                            UB = c(upper, b), H0 = H0,
    #                            control = ctrl, y = y, x_vec = x,
    #                            ga_shape = ga_shape,
    #                            ga_rate = ga_rate,
    #                            is.sig.par = is.sig.par)
    
    M_const <- max(exp(-opt_res_t$values[length(opt_res_t$values)]))
    
    
    # sample_t <- matrix(0L, 1, D)
    sample_t <- matrix(0L, n_sample, 1)
    
    s <- 1
    while (s <= n_sample) {
        t1_star <- stats::runif(1, min = a, max = b)
        candi_den <- 1 / (b - a)
        
        log_den <- -log_mar_lik_gp_der(theta = theta_new, y = y, 
                                       x_vec = x, der_vec = t1_star,
                                       H0 = H0,
                                       ga_shape = ga_shape, 
                                       ga_rate = ga_rate,
                                       is.sig.par = is.sig.par)
        # if (shape1 != 1 || shape2 != 1) {
        #     log_prior_den <-  log_dgbeta(t1_star, shape1 = shape1,
        #                                  shape2 = shape2, a = a, b = b)
        # }
        
        if (n_mixture == 1) {
            if (shape1 != 1 || shape2 != 1) {
                log_prior_den <- log_dgbeta(t1_star, shape1 = shape1,
                                                   shape2 = shape2,
                                                   a = a, b = b)
            }
        } else {
            prior_den <- 0
            for (m in 1:n_mixture) {
                prior_den <- prior_den +
                    mixture_prob[m] * exp(log_dgbeta(t1_star, 
                                                            shape1 = shape1[m],
                                                            shape2 = shape2[m],
                                                            a = a, b = b))
            }
            log_prior_den <- log(prior_den)
        }
        
        
        den <- exp(log_den + log_prior_den)
        if (stats::runif(1) * (M_const / tune) < (den / candi_den)) {
            sample_t[s, ] <- t1_star
            cat("final draw:", s, "\r")
            s <- s + 1
        }
    }
    cat("\n")
    print("Done!")
    colnames(thetas) <- c("sigma", "tau", "h")
    
    return(list(sample_t = sample_t, thetas = thetas[1:count, ]))
}

# k <- 1
# y = YY[[k]]$y
# x = YY[[k]]$x
# H0 = H0_diff[[k]]
# theta_init = c(1, 1, 1)
# epsilon = 1e-4
# D = 100
# a = 0
# b = 2
# n_sample = 1000
# max_iter = 100
# lower = c(0.0001, 0.0001, 0.0001)
# upper = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001)
# shape1 = 1
# shape2 = 1
# ctrl = list(TOL = 1e-5, trace = 0)
# is.sig.par = FALSE
# nu = 1.5
# is.h.par = FALSE
# init_t = 0.5

stochastic_em_dgp_matern <- function(y, x, H0, theta_init = c(1, 1, 1), epsilon = 1e-4, 
                              D = 100, a = 0, b = 2, n_sample = 4000, max_iter = 100,
                              lower = c(0.001, 0.001, 0.001),
                              upper = c(1/0.001, 1/0.001, 1/0.001), init_t = 0.5, 
                              shape1 = 1, shape2 = 1,
                              ctrl = list(TOL = 1e-5, trace = 0),
                              ga_shape = 5, ga_rate = 5,
                              a_h = 1, b_h = 1, 
                              nu = 1.5,
                              is.sig.par = TRUE,
                              is.h.par = TRUE) {
    ################################
    # Arguments
    
    ################################
    thetas <- matrix(0L, max_iter, length(theta_init))
    
    print("Stochastic EM new (one t) Begins")
    
    eps <- epsilon + 1
    count <- 1
    
    # Initialization of parameters
    thetas[1, ] <- theta_k <- theta_init
    theta_new <- theta_k
    if (shape1 == 1 && shape2 == 1) {
        log_prior_den <- log_dgbeta((a + b) / 2, shape1 = shape1, 
                                    shape2 = shape2, a = a, b = b)
    }
    
    
    while(eps > epsilon) {
        
        if (count == 100) {
            print("Does not converge")
            break
        }
        
        # ******** Stochastic E-step ************ #
        # Sampling posterior distributuion of t given theta = (sig, tau, h)
        # print("1")
        if (count == 1) {
            sample_t <- matrix(stats::runif(D, min = a, max = b), 1, D)
        } else {
            sample_t <- matrix(0L, 1, D)
            # print("2")
            opt_res_t <- Rsolnp::solnp(pars = c(theta_k, init_t),
                                       fun = log_mar_lik_gp_der_t_matern,
                                       LB = c(lower, a),
                                       UB = c(upper, b), H0 = H0, nu = nu,
                                       control = ctrl, y = y, x_vec = x,
                                       ga_shape = ga_shape, ga_rate = ga_rate,
                                       is.sig.par = is.sig.par)
            
            M_const <- max(exp(-opt_res_t$values[length(opt_res_t$values)]))
            
            s <- 1
            print("Sampling begins")
            while (s <= D) {
                t1_star <- stats::runif(1, min = a, max = b)
                candi_den <- 1 / (b - a)
                
                log_den <- -log_mar_lik_gp_der_matern(theta = theta_k, y = y, 
                                               x_vec = x, 
                                               der_vec = t1_star,
                                               H0 = H0, nu = nu,
                                               ga_shape = ga_shape, 
                                               ga_rate = ga_rate,
                                               is.sig.par = is.sig.par)   
                if (shape1 != 1 || shape2 != 1) {
                    log_prior_den <- log_dgbeta(t1_star, shape1 = shape1,
                                                shape2 = shape2, a = a, b = b)
                }
                
                den <- exp(log_den + log_prior_den)
                # (den / candi_den)
                if (stats::runif(1) * (M_const / 1) < (den / candi_den)) {
                    sample_t[[s]] <- t1_star
                    cat("no. draw:", s, "\r")
                    s <- s + 1
                }
            }
        }
        # print("3")
        # ******** M-step for theta and sigma ************ #
        res <- Rsolnp::solnp(pars = theta_k, fun = log_mar_lik_gp_der_mc_matern,
                             LB = lower, control = ctrl, y = y, x_vec = x, 
                             H0 = H0, der_mc_mat = sample_t, nu = nu,
                             ga_shape = ga_shape, ga_rate = ga_rate,
                             is.sig.par = is.sig.par,
                             is.h.par = is.h.par,
                             a_h = a_h, b_h = b_h)
        theta_k <- res$par
        mar_post_k <- res$values[length(res$values)]
        # print("4")
        # ******** Update epsilon and parameters ************ #
        eps <- sum((theta_new - theta_k) ^ 2)
        theta_new <- theta_k
        count <- count + 1
        print(paste("epsilon =", round(eps, 4), " count", count))
        print(paste("(increasing) marginal posterior =", -mar_post_k))
        thetas[count, ] <- theta_new
    }
    # print("5")
    opt_res_t <- Rsolnp::solnp(pars = c(theta_new, sum(sample_t) / D),
                               fun = log_mar_lik_gp_der_t_matern,
                               LB = c(lower, a),
                               UB = c(upper, b), H0 = H0, nu = nu,
                               control = ctrl, y = y, x_vec = x,
                               ga_shape = ga_shape,
                               ga_rate = ga_rate,
                               is.sig.par = is.sig.par)
    
    # opt_res_t <- Rsolnp::solnp(pars = c(theta_new, init_t),
    #                            fun = log_mar_lik_gp_der_t,
    #                            LB = c(lower, a),
    #                            UB = c(upper, b), H0 = H0,
    #                            control = ctrl, y = y, x_vec = x,
    #                            ga_shape = ga_shape,
    #                            ga_rate = ga_rate,
    #                            is.sig.par = is.sig.par)
    
    M_const <- max(exp(-opt_res_t$values[length(opt_res_t$values)]))
    
    
    # sample_t <- matrix(0L, 1, D)
    sample_t <- matrix(0L, n_sample, 1)
    
    s <- 1
    while (s <= n_sample) {
        t1_star <- stats::runif(1, min = a, max = b)
        candi_den <- 1 / (b - a)
        
        log_den <- -log_mar_lik_gp_der_matern(theta = theta_new, y = y, 
                                       x_vec = x, der_vec = t1_star,
                                       H0 = H0, nu = nu,
                                       ga_shape = ga_shape, 
                                       ga_rate = ga_rate,
                                       is.sig.par = is.sig.par)
        if (shape1 != 1 || shape2 != 1) {
            log_prior_den <-  log_dgbeta(t1_star, shape1 = shape1,
                                         shape2 = shape2, a = a, b = b)
        }
        
        den <- exp(log_den + log_prior_den)
        if (stats::runif(1) * (M_const / 1) < (den / candi_den)) {
            sample_t[s, ] <- t1_star
            cat("final draw:", s, "\r")
            s <- s + 1
        }
    }
    cat("\n")
    print("Done!")
    colnames(thetas) <- c("sigma", "tau", "h")
    
    return(list(sample_t = sample_t, thetas = thetas[1:count, ]))
}
# y = YY[[k]]$y
# x = YY[[k]]$x
# H0 = H0_diff[[k]]
# theta_init = c(1, 1, 1)
# epsilon = 1e-4
# D = 100
# a = 0
# b = 2
# n_sample = 1000
# max_iter = 100
# lower = c(0.0001, 0.0001, 0.0001)
# upper = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001)
# shape1 = 1
# shape2 = 1
# ctrl = list(TOL = 1e-5, trace = 0)
# is.sig.par = FALSE
# is.h.par = FALSE

stochastic_em_dgp_new <- function(y, x, H0, theta_init = c(1, 1, 1), epsilon = 1e-4, 
                              D = 100, a = 0, b = 2, n_sample = 4000, max_iter = 100,
                              lower = c(0.001, 0.001, 0.001),
                              upper = c(1/0.001, 1/0.001, 1/0.001), init_t = 0.5, 
                              shape1 = 1, shape2 = 1,
                              ctrl = list(TOL = 1e-5, trace = 0),
                              ga_shape = 5, ga_rate = 5,
                              a_h = 1, b_h = 1, 
                              is.sig.par = TRUE,
                              is.h.par = TRUE) {
    ################################
    # Arguments
    
    ################################
    thetas <- matrix(0L, max_iter, length(theta_init))
    
    print("Stochastic EM new (one t) Begins")
    
    eps <- epsilon + 1
    count <- 1
    
    # Initialization of parameters
    thetas[1, ] <- theta_k <- theta_init
    theta_new <- theta_k
    if (shape1 == 1 && shape2 == 1) {
        log_prior_den <- log_dgbeta((a + b) / 2, shape1 = shape1, 
                                    shape2 = shape2, a = a, b = b)
    }
    
    
    while(eps > epsilon) {
        
        if (count == 100) {
            print("Does not converge")
            break
        }
        
        # ******** Stochastic E-step ************ #
        # Sampling posterior distributuion of t given theta = (sig, tau, h)
        # print("1")
        if (count == 1) {
            sample_t <- matrix(stats::runif(D, min = a, max = b), 1, D)
        } else {
            sample_t <- matrix(0L, 1, D)
            # print("2")
            opt_res_t <- Rsolnp::solnp(pars = c(theta_k, init_t),
                                       fun = log_mar_lik_gp_der_t,
                                       LB = c(lower, a),
                                       UB = c(upper, b), H0 = H0, 
                                       control = ctrl, y = y, x_vec = x,
                                       ga_shape = ga_shape, ga_rate = ga_rate,
                                       is.sig.par = is.sig.par)
            
            M_const <- max(exp(-opt_res_t$values[length(opt_res_t$values)]))
            
            s <- 1
            print("Sampling begins")
            while (s <= D) {
                t1_star <- stats::runif(1, min = a, max = b)
                candi_den <- 1 / (b - a)
                
                log_den <- -log_mar_lik_gp_der(theta = theta_k, y = y, 
                                               x_vec = x, 
                                               der_vec = t1_star,
                                               H0 = H0,
                                               ga_shape = ga_shape, 
                                               ga_rate = ga_rate,
                                               is.sig.par = is.sig.par)   
                if (shape1 != 1 || shape2 != 1) {
                    log_prior_den <- log_dgbeta(t1_star, shape1 = shape1,
                                                shape2 = shape2, a = a, b = b)
                }
                
                den <- exp(log_den + log_prior_den)
                # (den / candi_den)
                if (stats::runif(1) * (M_const / 1) < (den / candi_den)) {
                    sample_t[[s]] <- t1_star
                    cat("no. draw:", s, "\r")
                    s <- s + 1
                }
            }
        }
        # print("3")
        # ******** M-step for theta and sigma ************ #
        res <- Rsolnp::solnp(pars = theta_k, fun = log_mar_lik_gp_der_mc_new2,
                             LB = lower, control = ctrl, y = y, x_vec = x, 
                             H0 = H0, der_mc_mat = sample_t,
                             ga_shape = ga_shape, ga_rate = ga_rate,
                             is.sig.par = is.sig.par,
                             is.h.par = is.h.par,
                             a_h = a_h, b_h = b_h)
        theta_k <- res$par
        mar_post_k <- res$values[length(res$values)]
        # print("4")
        # ******** Update epsilon and parameters ************ #
        eps <- sum((theta_new - theta_k) ^ 2)
        theta_new <- theta_k
        count <- count + 1
        print(paste("epsilon =", round(eps, 4), " count", count))
        print(paste("(increasing) marginal posterior =", -mar_post_k))
        thetas[count, ] <- theta_new
    }
    # print("5")
    opt_res_t <- Rsolnp::solnp(pars = c(theta_new, sum(sample_t) / D),
                               fun = log_mar_lik_gp_der_t,
                               LB = c(lower, a),
                               UB = c(upper, b), H0 = H0,
                               control = ctrl, y = y, x_vec = x,
                               ga_shape = ga_shape,
                               ga_rate = ga_rate,
                               is.sig.par = is.sig.par)
    
    # opt_res_t <- Rsolnp::solnp(pars = c(theta_new, init_t),
    #                            fun = log_mar_lik_gp_der_t,
    #                            LB = c(lower, a),
    #                            UB = c(upper, b), H0 = H0,
    #                            control = ctrl, y = y, x_vec = x,
    #                            ga_shape = ga_shape,
    #                            ga_rate = ga_rate,
    #                            is.sig.par = is.sig.par)
    
    M_const <- max(exp(-opt_res_t$values[length(opt_res_t$values)]))
    
    
    # sample_t <- matrix(0L, 1, D)
    sample_t <- matrix(0L, n_sample, 1)
    
    s <- 1
    while (s <= n_sample) {
        t1_star <- stats::runif(1, min = a, max = b)
        candi_den <- 1 / (b - a)
        
        log_den <- -log_mar_lik_gp_der(theta = theta_new, y = y, 
                                       x_vec = x, der_vec = t1_star,
                                       H0 = H0,
                                       ga_shape = ga_shape, 
                                       ga_rate = ga_rate,
                                       is.sig.par = is.sig.par)
        if (shape1 != 1 || shape2 != 1) {
            log_prior_den <-  log_dgbeta(t1_star, shape1 = shape1,
                                         shape2 = shape2, a = a, b = b)
        }
        
        den <- exp(log_den + log_prior_den)
        if (stats::runif(1) * (M_const / 1) < (den / candi_den)) {
            sample_t[s, ] <- t1_star
            cat("final draw:", s, "\r")
            s <- s + 1
        }
    }
    cat("\n")
    print("Done!")
    colnames(thetas) <- c("sigma", "tau", "h")
    
    return(list(sample_t = sample_t, thetas = thetas[1:count, ]))
}





stochastic_em_dgp_2_par <- function(y, x, H0, 
                                    theta_init = c(1, 1), 
                                    sig_init = 1,
                                    epsilon = 1e-4, 
                                    D = 100, a = 0, b = 2, n_sample = 4000, max_iter = 100,
                                    lower = c(0.001, 0.001, 0.001),
                                    upper = c(1/0.001, 1/0.001, 1/0.001), init_t = 0.5, 
                                    shape1 = 1, shape2 = 1,
                                    ctrl = list(TOL = 1e-5, trace = 0),
                                    ga_shape = 5, ga_rate = 5,
                                    a_h = 1, b_h = 1, 
                                    is.sig.par = TRUE,
                                    is.h.par = TRUE) {
    ################################
    # Arguments
    
    ################################
    thetas <- matrix(0L, max_iter, length(theta_init))
    sigs <- rep(0, max_iter)
    taus <- rep(0, max_iter)
    n <- length(y)
    print("Stochastic EM new (one t) Begins")
    
    eps <- epsilon + 1
    count <- 1
    
    # Initialization of parameters
    thetas[1, ] <- theta_k <- theta_init
    theta_new <- theta_k
    
    sigs[1] <- sig_k <- sig_init
    sig_new <- sig_k
    
    taus[1] <- tau_k <- sigs[1] * thetas[1, 1]
    tau_new <- tau_k
    
    
    if (shape1 == 1 && shape2 == 1) {
        log_prior_den <- log_dgbeta((a + b) / 2, shape1 = shape1, 
                                    shape2 = shape2, a = a, b = b)
    }
    
    
    while(eps > epsilon) {
        
        if (count == 100) {
            print("Does not converge")
            break
        }
        
        # ******** Stochastic E-step ************ #
        # Sampling posterior distributuion of t given theta = (sig, tau, h)
        # print("1")
        if (count == 1) {
            sample_t <- matrix(stats::runif(D, min = a, max = b), 1, D)
        } else {
            sample_t <- matrix(0L, 1, D)
            # print("2")
            opt_res_t <- Rsolnp::solnp(pars = c(sig_k, tau_k, theta_k[2], init_t),
                                       fun = log_mar_lik_gp_der_t,
                                       LB = c(lower, a),
                                       UB = c(upper, b), H0 = H0, 
                                       control = ctrl, y = y, x_vec = x,
                                       ga_shape = ga_shape, ga_rate = ga_rate,
                                       is.sig.par = is.sig.par)
            
            M_const <- max(exp(-opt_res_t$values[length(opt_res_t$values)]))
            
            s <- 1
            print("Sampling begins")
            while (s <= D) {
                t1_star <- stats::runif(1, min = a, max = b)
                candi_den <- 1 / (b - a)
                
                log_den <- -log_mar_lik_gp_der(theta = c(sig_k, tau_k, theta_k[2]), y = y, 
                                               x_vec = x, 
                                               der_vec = t1_star,
                                               H0 = H0,
                                               ga_shape = ga_shape, 
                                               ga_rate = ga_rate,
                                               is.sig.par = is.sig.par)   
                if (shape1 != 1 || shape2 != 1) {
                    log_prior_den <- log_dgbeta(t1_star, shape1 = shape1,
                                                shape2 = shape2, a = a, b = b)
                }
                
                den <- exp(log_den + log_prior_den)
                # (den / candi_den)
                if (stats::runif(1) * (M_const / 1) < (den / candi_den)) {
                    sample_t[[s]] <- t1_star
                    cat("no. draw:", s, "\r")
                    s <- s + 1
                }
            }
        }
        # print("3")
        # ******** M-step for theta and sigma ************ #
        
        ## (lambda, h) (tau2 = sig2 * lambda)
        
        # theta, y, x_vec, der_mc_mat, H0,
        # a_h = 1, b_h = 1, sig2, 
        # is.h.par = TRUE
        res <- Rsolnp::solnp(pars = theta_k, 
                             fun = log_mar_lik_gp_der_mc_2_par,
                             LB = lower[1:2], control = ctrl, 
                             y = y, x_vec = x, 
                             H0 = H0, der_mc_mat = sample_t,
                             is.h.par = is.h.par, sig2 = sig2,
                             a_h = a_h, b_h = b_h)
        theta_k <- res$par
        mar_post_k <- res$values[length(res$values)]
        

        
        
        ## sig
        Kff_no_tau <- se_ker_no_tau(H0 = H0, h = theta_k[2])
        A <- Kff_no_tau + diag(1 / theta_k[1], n)
        
        Kdf_lst_no_tau <- as.list(data.frame(apply(sample_t, 2, computeCovDer1, 
                                                   idx2 = x,
                                                   ker_der1_fcn = se_ker_der1_no_tau,
                                                   h = theta_k[2])))
        
        Kdf_lst_no_tau <- lapply(Kdf_lst_no_tau, function(x){
            matrix(x, nrow = nrow(sample_t), ncol = n)
        })
        
        B_arry_no_tau <- sapply(Kdf_lst_no_tau, create_B, h = theta_k[2], 
                                simplify = "array")
        
        avgB <- apply(B_arry_no_tau, c(1, 2), mean)
        sig_k <- sqrt(update_sig2(y = y, lambda = theta_k[1], A = A, avgB = avgB, 
                                  is.sig.par = is.sig.par, ga_shape = 5, 
                                  ga_rate = 5))
        
        tau_k <- sig_k * sqrt(theta_k[1])
        
        
        # print("4")
        # ******** Update epsilon and parameters ************ #
        eps <- sum((theta_new - theta_k) ^ 2) + (sig_new - sig_k) ^ 2
        theta_new <- theta_k
        sig_new <- sig_k
        tau_new <- tau_k
        # h_new <- theta_new[2]
        count <- count + 1
        print(paste("epsilon =", round(eps, 4), " count", count))
        print(paste("(increasing) marginal posterior =", -mar_post_k))
        thetas[count, ] <- theta_new
        sigs[count] <- sig_new
        taus[count] <- tau_new
    }
    # print("5")
    opt_res_t <- Rsolnp::solnp(pars = c(sig_new, tau_new, h_new, sum(sample_t) / D),
                               fun = log_mar_lik_gp_der_t,
                               LB = c(lower, a),
                               UB = c(upper, b), H0 = H0,
                               control = ctrl, y = y, x_vec = x,
                               ga_shape = ga_shape,
                               ga_rate = ga_rate,
                               is.sig.par = is.sig.par)
    
    # opt_res_t <- Rsolnp::solnp(pars = c(theta_new, init_t),
    #                            fun = log_mar_lik_gp_der_t,
    #                            LB = c(lower, a),
    #                            UB = c(upper, b), H0 = H0,
    #                            control = ctrl, y = y, x_vec = x,
    #                            ga_shape = ga_shape,
    #                            ga_rate = ga_rate,
    #                            is.sig.par = is.sig.par)
    
    M_const <- max(exp(-opt_res_t$values[length(opt_res_t$values)]))
    
    
    # sample_t <- matrix(0L, 1, D)
    sample_t <- matrix(0L, n_sample, 1)
    
    s <- 1
    while (s <= n_sample) {
        t1_star <- stats::runif(1, min = a, max = b)
        candi_den <- 1 / (b - a)
        
        log_den <- -log_mar_lik_gp_der(theta = c(sig_new, tau_new, h_new), y = y, 
                                       x_vec = x, der_vec = t1_star,
                                       H0 = H0,
                                       ga_shape = ga_shape, 
                                       ga_rate = ga_rate,
                                       is.sig.par = is.sig.par)
        if (shape1 != 1 || shape2 != 1) {
            log_prior_den <-  log_dgbeta(t1_star, shape1 = shape1,
                                         shape2 = shape2, a = a, b = b)
        }
        
        den <- exp(log_den + log_prior_den)
        if (stats::runif(1) * (M_const / 1) < (den / candi_den)) {
            sample_t[s, ] <- t1_star
            cat("final draw:", s, "\r")
            s <- s + 1
        }
    }
    cat("\n")
    print("Done!")
    # colnames(thetas) <- c("sigma", "tau", "h")
    
    return(list(sample_t = sample_t, 
                thetas = thetas[1:count, ],
                sigs = sigs[1:count],
                taus = taus[1:count]))
}






# microbenchmark(rep(0, 100), 
#                matrix(rnorm(100), 1, 100),
#                matrix(rep(0, 100), 1, 100), 
#                as.matrix(rep(0, 100)))
# mat <- matrix(rnorm(100), 1, 100)

# system.time(for(i in 1:100) {
#     res <- Rsolnp::solnp(pars = theta_k, fun = log_mar_lik_gp_der_mc_new,
#                          LB = lower, control = ctrl, y = y, x_vec = x, 
#                          H0 = H0, der_mc_mat = sample_t,
#                          ga_shape = ga_shape, ga_rate = ga_rate,
#                          is.sig.par = is.sig.par,
#                          is.h.par = is.h.par,
#                          a_h = a_h, b_h = b_h)
# })
# 
# system.time(for(i in 1:100) {
#     res <- Rsolnp::solnp(pars = theta_k, fun = log_mar_lik_gp_der_mc_new2,
#                          LB = lower, control = ctrl, y = y, x_vec = x, 
#                          H0 = H0, der_mc_mat = sample_t,
#                          ga_shape = ga_shape, ga_rate = ga_rate,
#                          is.sig.par = is.sig.par,
#                          is.h.par = is.h.par,
#                          a_h = a_h, b_h = b_h)
# })
# 
# system.time(for(i in 1:100) {
#     res <- Rsolnp::solnp(pars = theta_k, fun = log_mar_lik_gp_der_mc,
#                          LB = lower, control = ctrl, y = y, x_vec = x, 
#                          H0 = H0, der_mc_mat = sample_t,
#                          ga_shape = ga_shape, ga_rate = ga_rate,
#                          is.sig.par = is.sig.par,
#                          is.h.par = is.h.par,
#                          a_h = a_h, b_h = b_h)
# })
# 
# # theta, y, x_vec, der_mc_mat, H0,
# # a_h = 1, b_h = 1, sig2, 
# # is.h.par = TRUE
# theta_k <- c(1, 1)
# system.time(for(i in 1:100) {
#     res <- Rsolnp::solnp(pars = theta_k, fun = log_mar_lik_gp_der_mc_2_par,
#                          LB = lower[1:2], control = ctrl, y = y, x_vec = x, 
#                          H0 = H0, der_mc_mat = sample_t, 
#                          is.h.par = is.h.par, sig2 = 0.2^2,
#                          a_h = a_h, b_h = b_h)
# })
# 
# theta_k <- c(1, 1)
# system.time(for(i in 1:100) {
#     res <- nlminb(start = theta_k, log_mar_lik_gp_der_mc_2_par, 
#                   lower = lower[1:2], y = y, x_vec = x, 
#                   H0 = H0, der_mc_mat = sample_t, 
#                   is.h.par = is.h.par, sig2 = 1,
#                   a_h = a_h, b_h = b_h)
#     # res <- optim(par = theta_k, fn = log_mar_lik_gp_der_mc_2_par,
#     #              method = "L-BFGS-B",
#     #              lower = lower[1:2], control = list(trace = 0), 
#     #              y = y, x_vec = x, 
#     #              H0 = H0, der_mc_mat = sample_t, 
#     #              is.h.par = is.h.par, sig2 = 1,
#     #              a_h = a_h, b_h = b_h)
# })

