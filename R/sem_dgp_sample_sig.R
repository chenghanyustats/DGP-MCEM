sem_dgp_sample_sig <- function(y, x, H0, theta_init = c(1, 1), epsilon = 1e-4, 
                               D = 100, a = 0, b = 2, n_sample = 4000, 
                               max_iter = 100,
                               lower = c(0.001, 0.001),
                               upper = c(1/0.001, 1/0.001), init_t = 0.5, 
                               init_sig = 1,
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
    
    sigs <- rep(0, max_iter)
    taus <- rep(0, max_iter)
    n <- length(y)
    
    
    print("SEM with sampling t and sig Begins")
    
    eps <- epsilon + 1
    count <- 1
    
    n_mixture <- length(shape1)
    
    
    # Initialization of parameters
    thetas[1, ] <- theta_k <- theta_init
    theta_new <- theta_k

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
        
        #==============================================
        # ******** Stochastic E-step ************ #
        #==============================================
        

        print("1")
        if (count == 1) {
            sample_t <- matrix(stats::runif(D, min = a, max = b), 1, D)
            sample_sig <- rep(init_sig, D)
        } else {
            sample_t <- matrix(0L, 1, D)
            sample_sig <- rep(init_sig, D)
            print("2")
            
            opt_res_t <- Rsolnp::solnp(pars = c(theta_k, init_sig, init_t),
                                       fun = log_mar_lik_gp_der_t_sig,
                                       LB = c(lower, 0.001, a),
                                       UB = c(upper, 1/0.001, b), H0 = H0, 
                                       control = ctrl, y = y, x_vec = x)
            
            M_const <- max(exp(-opt_res_t$values[length(opt_res_t$values)]))
            
            s <- 1
            print("Sampling begins")
            while (s <= D) {
                # opt_res_t <- Rsolnp::solnp(pars = c(theta_k, init_t),
                #                            fun = log_mar_lik_gp_der_t_given_sig,
                #                            LB = c(lower, a),
                #                            UB = c(upper, b), H0 = H0, 
                #                            control = ctrl, y = y, x_vec = x, 
                #                            sig = sample_sig[s])

                if (s == 1) {
                    sig <- init_sig
                } else {
                    sig <- sample_sig[s - 1]
                }
                
                #-----------------------------------------------
                # Sampling t
                #-----------------------------------------------
                t1_star <- stats::runif(1, min = a, max = b)
                candi_den <- 1 / (b - a)
                

                log_den <- -log_mar_lik_gp_der_given_sig(theta = theta_k, y = y, 
                                                         x_vec = x, 
                                                         der_vec = t1_star,
                                                         H0 = H0,
                                                         sig = sig)   
                
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
                # cat("den", den)
                
                # stats::runif(1) * (M_const / tune)
                
                if (stats::runif(1) * (M_const / tune) < (den / candi_den)) {
                    sample_t[[s]] <- t1_star
                    
                    sample_sig[s] <- sqrt(sample_sig2(y = y, x_vec = x,
                                                      sample_t = sample_t[[s]], 
                                                      theta = theta_k, H0 = H0, 
                                                      ga_shape = ga_shape, 
                                                      ga_rate = ga_rate))
                    
                    cat("no. draw:", s, "\r")
                    s <- s + 1
                }
            }
        }
        
        print("3")
        #==============================================
        # ********   M-step for theta    ************ #
        #==============================================
        res <- Rsolnp::solnp(pars = theta_k, fun = log_mar_lik_gp_der_mc_sig2,
                             LB = lower, control = ctrl, y = y, x_vec = x, 
                             H0 = H0, der_mc_mat = sample_t,
                             is.h.par = is.h.par, sample_sig = sample_sig,
                             a_h = a_h, b_h = b_h)
        theta_k <- res$par
        mar_post_k <- res$values[length(res$values)]
        
        print("4")
        # ******** Update epsilon and parameters ************ #
        
        eps <- sum((theta_new - theta_k) ^ 2)
        theta_new <- theta_k
        count <- count + 1
        print(paste("epsilon =", round(eps, 4), " count", count))
        print(paste("(increasing) marginal posterior =", -mar_post_k))
        thetas[count, ] <- theta_new
    }
    
    print("5")
    
    # opt_res_t <- Rsolnp::solnp(pars = c(theta_k, init_sig, init_t),
    #                            fun = log_mar_lik_gp_der_t_sig,
    #                            LB = c(lower, 0.001, a),
    #                            UB = c(upper, 1/0.001, b), H0 = H0, 
    #                            control = ctrl, y = y, x_vec = x)
    
    
    opt_res_t <- Rsolnp::solnp(pars = c(theta_new, sum(sample_sig) / D, 
                                        sum(sample_t) / D),
                               fun = log_mar_lik_gp_der_t_sig,
                               LB = c(lower, 0.001, a),
                               UB = c(upper, 1/0.001, b), H0 = H0,
                               control = ctrl, y = y, x_vec = x)
    
    
    M_const <- max(exp(-opt_res_t$values[length(opt_res_t$values)]))
    
    
    # sample_t <- matrix(0L, 1, D)
    sample_t <- matrix(0L, n_sample, 1)
    sample_sig <- rep(mean(sample_sig), n_sample)
    
    s <- 1
    while (s <= n_sample) {
        t1_star <- stats::runif(1, min = a, max = b)
        candi_den <- 1 / (b - a)
        
        
        if (s == 1) {
            sig <- init_sig
        } else {
            sig <- sample_sig[s - 1]
        }
        
        # log_den <- -log_mar_lik_gp_der(theta = theta_new, y = y, 
        #                                x_vec = x, der_vec = t1_star,
        #                                H0 = H0,
        #                                ga_shape = ga_shape, 
        #                                ga_rate = ga_rate,
        #                                is.sig.par = is.sig.par)
        
        log_den <- -log_mar_lik_gp_der_given_sig(theta = theta_k, y = y, 
                                                 x_vec = x, 
                                                 der_vec = t1_star,
                                                 H0 = H0,
                                                 sig = sig)   
        
        
        
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
            
            sample_sig[s] <- sqrt(sample_sig2(y = y, 
                                              sample_t = sample_t[s, ], 
                                              theta = theta_k, H0 = H0, 
                                              ga_shape = ga_shape, 
                                              ga_rate = ga_rate))
            
            cat("final draw:", s, "\r")
            s <- s + 1
        }
    }
    cat("\n")
    print("Done!")
    colnames(thetas) <- c("tau", "h")
    
    return(list(sample_t = sample_t, sample_sig, 
                thetas = thetas[1:count, ]))
}

 

sem_dgp_sample_sig_mh <- function(y, x, H0, theta_init = c(1, 1), epsilon = 1e-4, 
                                  D = 100, a = 0, b = 2, n_sample = 4000, 
                                  max_iter = 100,
                                  lower = c(0.001, 0.001),
                                  upper = c(1/0.001, 1/0.001), init_t = 0.5, 
                                  init_sig = 1,
                                  shape1 = 1, shape2 = 1,
                                  ctrl = list(TOL = 1e-5, trace = 0),
                                  ga_shape = 1/2, ga_rate = 1/2,
                                  a_h = 1, b_h = 1, mixture_prob = c(0.5, 0.5),
                                  is.h.par = FALSE) {
    ################################
    # Arguments
    
    ################################
    thetas <- matrix(0L, max_iter, length(theta_init))
    
    print("SEM with sampling t and sig Begins")
    
    eps <- epsilon + 1
    count <- 1
    
    n_mixture <- length(shape1)
    
    # Initialization of parameters
    thetas[1, ] <- theta_k <- theta_init
    theta_new <- theta_k
    
    while(eps > epsilon) {
        
        if (count == 100) {
            print("Does not converge")
            break
        }
        
        sample_sig <- rep(init_sig, D)
        
        Kff <- se_ker(H0 = H0, tau = theta_k[1], h = theta_k[2])
        
        #==============================================
        # ******** Stochastic E-step ************ #
        #==============================================
        
        if (count == 1) {
            sample_t <- matrix(stats::runif(D, min = a, max = b), 1, D)
            
        } else {
            sample_t <- matrix(init_t, 1, D)
            
            s <- 1
            print("Sampling begins")
            while (s <= D) {
                
                if (s == 1) {
                    sig <- init_sig
                    t1_old <- runif(1, min = a, max = b)
                } else {
                    sig <- sample_sig[s - 1]
                    t1_old <- sample_t[[s - 1]]
                }
                
                
                #-----------------------------------------------
                # Sampling t
                #-----------------------------------------------
                t1_star <- stats::runif(1, min = a, max = b)
                
                # candi_den <- 1 / (b - a)
                
                log_post_den_res <- log_post_den_given_sig(t = t1_star, x = x, 
                                                           theta = theta_k, 
                                                           y = y, Kff = Kff, 
                                                           sig = sig,
                                                           shape1 = shape1, 
                                                           shape2 = shape2, 
                                                           a = a, b = b, 
                                                           mixture_prob = mixture_prob)
                
                log_post_den <- log_post_den_res$log_post_den
                
                log_post_den_old_res <- log_post_den_given_sig(t = t1_old, x = x, 
                                                               theta = theta_k, 
                                                               y = y, Kff = Kff, 
                                                               sig = sig,
                                                               shape1 = shape1, 
                                                               shape2 = shape2, 
                                                               a = a, b = b, 
                                                               mixture_prob = mixture_prob)
                log_post_den_old <- log_post_den_old_res$log_post_den
                
                ## MH-step
                if (log(runif(1)) < (log_post_den - log_post_den_old)) {
                    sample_t[[s]] <- t1_star
                    sample_sig[s] <- sqrt(sample_sig2(y = y, x_vec = x,
                                                      theta = theta_k, 
                                                      Kff = Kff, 
                                                      Kdf = log_post_den_res$Kdf,
                                                      Kdd = log_post_den_res$Kdd,
                                                      ga_shape = ga_shape, 
                                                      ga_rate = ga_rate))
                } else {
                    sample_t[[s]] <- t1_old
                    sample_sig[s] <- sqrt(sample_sig2(y = y, x_vec = x,
                                                      theta = theta_k, 
                                                      Kff = Kff, 
                                                      Kdf = log_post_den_old_res$Kdf,
                                                      Kdd = log_post_den_old_res$Kdd,
                                                      ga_shape = ga_shape, 
                                                      ga_rate = ga_rate))
                }
                cat("no. draw:", s, "\r")
                s <- s + 1
            }
        }
        
        
        #==============================================
        # ********   M-step for theta    ************ #
        #==============================================
        res <- Rsolnp::solnp(pars = theta_k, fun = log_mar_lik_gp_der_mc_sig2,
                             LB = lower, control = ctrl, y = y, x_vec = x, 
                             H0 = H0, der_mc_mat = sample_t,
                             is.h.par = is.h.par, sample_sig = sample_sig,
                             a_h = a_h, b_h = b_h)
        theta_k <- res$par
        mar_post_k <- res$values[length(res$values)]
        
        
        # ******** Update epsilon and parameters ************ #
        
        eps <- sum((theta_new - theta_k) ^ 2)
        theta_new <- theta_k
        count <- count + 1
        print(paste("epsilon =", round(eps, 4), " count", count))
        print(paste("(increasing) marginal posterior =", -mar_post_k))
        thetas[count, ] <- theta_new
    }
    
    ## Final samples
    sample_t <- matrix(0L, n_sample, 1)
    sample_sig <- rep(mean(sample_sig), n_sample)
    
    s <- 1
    while (s <= n_sample) {
        
        if (s == 1) {
            sig <- sample_sig[s]
            t1_old <- runif(1, min = a, max = b)
            
        } else {
            sig <- sample_sig[s - 1]
            t1_old <- sample_t[[s - 1]]
        }
        
        
        #-----------------------------------------------
        # Sampling t
        #-----------------------------------------------
        t1_star <- stats::runif(1, min = a, max = b)
        
        # candi_den <- 1 / (b - a)
        
        log_post_den_res <- log_post_den_given_sig(t = t1_star, x = x, 
                                                   theta = theta_k, 
                                                   y = y, Kff = Kff, 
                                                   sig = sig,
                                                   shape1 = shape1, 
                                                   shape2 = shape2, 
                                                   a = a, b = b, 
                                                   mixture_prob = mixture_prob)
        
        log_post_den <- log_post_den_res$log_post_den
        
        log_post_den_old_res <- log_post_den_given_sig(t = t1_old, x = x, 
                                                       theta = theta_k, 
                                                       y = y, Kff = Kff, 
                                                       sig = sig,
                                                       shape1 = shape1, 
                                                       shape2 = shape2, 
                                                       a = a, b = b, 
                                                       mixture_prob = mixture_prob)
        log_post_den_old <- log_post_den_old_res$log_post_den
        
        
        ## MH-step
        if (log(runif(1)) < (log_post_den - log_post_den_old)) {
            sample_t[[s]] <- t1_star
            sample_sig[s] <- sqrt(sample_sig2(y = y, x_vec = x,
                                              theta = theta_k, 
                                              Kff = Kff, 
                                              Kdf = log_post_den_res$Kdf,
                                              Kdd = log_post_den_res$Kdd,
                                              ga_shape = ga_shape, 
                                              ga_rate = ga_rate))
        } else {
            sample_t[[s]] <- t1_old
            sample_sig[s] <- sqrt(sample_sig2(y = y, x_vec = x,
                                              theta = theta_k, 
                                              Kff = Kff, 
                                              Kdf = log_post_den_old_res$Kdf,
                                              Kdd = log_post_den_old_res$Kdd,
                                              ga_shape = ga_shape, 
                                              ga_rate = ga_rate))
        }
        
        cat("final draw:", s, "\r")
        s <- s + 1
    }
    
    
    cat("\n")
    print("Done!")
    colnames(thetas) <- c("tau", "h")
    
    return(list(sample_t = sample_t, 
                sample_sig = sample_sig, 
                thetas = thetas[1:count, ]))
}

mcem_dgp <- function(y, x, H0, theta_init = c(1, 1), epsilon = 1e-4, 
                     D = 100, n_sample = 4000, a = 0, b = 2, 
                     max_iter = 100,
                     lower = c(0.001, 0.001),
                     upper = c(1/0.001, 1/0.001), init_t = 0.5, 
                     init_sig = 1,
                     shape1 = 1, shape2 = 1,
                     ctrl = list(TOL = 1e-5, trace = 0),
                     ga_shape = 1/2, ga_rate = 1/2,
                     a_h = 1, b_h = 1, mixture_prob = c(0.5, 0.5),
                     is.h.par = FALSE) {
    ################################
    # Arguments
    
    ################################
    thetas <- matrix(0L, max_iter, length(theta_init))
    
    print("MCEM with sampling t and sig Begins")
    
    eps <- epsilon + 1
    count <- 1
    
    n_mixture <- length(shape1)
    
    # Initialization of parameters
    thetas[1, ] <- theta_k <- theta_init
    theta_new <- theta_k
    
    while(eps > epsilon) {
        
        if (count == 100) {
            print("Does not converge")
            break
        }
        
        sample_sig <- rep(init_sig, n_sample)
        
        Kff <- se_ker(H0 = H0, tau = theta_k[1], h = theta_k[2])
        
        #==============================================
        # ********   Monte Carlo E-step  ************ #
        #==============================================
        
        if (count == 1) {
            sample_t <- matrix(stats::runif(n_sample, min = a, max = b), 
                               1, n_sample)
            
        } else {
            sample_t <- matrix(init_t, 1, n_sample)
            
            s <- 1
            print("Sampling begins")
            while (s <= n_sample) {
                
                if (s == 1) {
                    sig <- init_sig
                    t1_old <- runif(1, min = a, max = b)
                } else {
                    sig <- sample_sig[s - 1]
                    t1_old <- sample_t[[s - 1]]
                }
                
                
                #-----------------------------------------------
                # Sampling t
                #-----------------------------------------------
                t1_star <- stats::runif(1, min = a, max = b)
                
                # candi_den <- 1 / (b - a)
                
                log_post_den_res <- log_post_den_given_sig(t = t1_star, x = x, 
                                                           theta = theta_k, 
                                                           y = y, Kff = Kff, 
                                                           sig = sig,
                                                           shape1 = shape1, 
                                                           shape2 = shape2, 
                                                           a = a, b = b, 
                                                           mixture_prob = mixture_prob)
                
                log_post_den <- log_post_den_res$log_post_den
                
                log_post_den_old_res <- log_post_den_given_sig(t = t1_old, x = x, 
                                                               theta = theta_k, 
                                                               y = y, Kff = Kff, 
                                                               sig = sig,
                                                               shape1 = shape1, 
                                                               shape2 = shape2, 
                                                               a = a, b = b, 
                                                               mixture_prob = mixture_prob)
                log_post_den_old <- log_post_den_old_res$log_post_den
                
                ## MH-step
                if (log(runif(1)) < (log_post_den - log_post_den_old)) {
                    sample_t[[s]] <- t1_star
                    sample_sig[s] <- sqrt(sample_sig2(y = y, x_vec = x,
                                                      theta = theta_k, 
                                                      Kff = Kff, 
                                                      Kdf = log_post_den_res$Kdf,
                                                      Kdd = log_post_den_res$Kdd,
                                                      ga_shape = ga_shape, 
                                                      ga_rate = ga_rate))
                } else {
                    sample_t[[s]] <- t1_old
                    sample_sig[s] <- sqrt(sample_sig2(y = y, x_vec = x,
                                                      theta = theta_k, 
                                                      Kff = Kff, 
                                                      Kdf = log_post_den_old_res$Kdf,
                                                      Kdd = log_post_den_old_res$Kdd,
                                                      ga_shape = ga_shape, 
                                                      ga_rate = ga_rate))
                }
                cat("no. draw:", s, "\r")
                s <- s + 1
            }
        }
        
        
        #==============================================
        # ********   M-step for theta    ************ #
        #==============================================
        
        sample_t_mstep <- t(apply(sample_t, 1, sample, size = D))
        sample_sig_mstep <- sample(sample_sig, size = D)
        
        res <- Rsolnp::solnp(pars = theta_k, fun = log_mar_lik_gp_der_mc_sig2,
                             LB = lower, control = ctrl, y = y, x_vec = x, 
                             H0 = H0, der_mc_mat = sample_t_mstep,
                             is.h.par = is.h.par, sample_sig = sample_sig_mstep,
                             a_h = a_h, b_h = b_h)
        theta_k <- res$par
        mar_post_k <- res$values[length(res$values)]
        
        
        # ******** Update epsilon and parameters ************ #
        
        eps <- sum((theta_new - theta_k) ^ 2)
        theta_new <- theta_k
        count <- count + 1
        print(paste("epsilon =", round(eps, 4), " count", count))
        print(paste("(increasing) marginal posterior =", -mar_post_k))
        thetas[count, ] <- theta_new
    }
    
    ## Final samples
    # sample_t <- matrix(0L, n_sample, 1)
    # sample_sig <- rep(mean(sample_sig), n_sample)
    # 
    # s <- 1
    # while (s <= n_sample) {
    #     
    #     if (s == 1) {
    #         sig <- sample_sig[s]
    #         t1_old <- runif(1, min = a, max = b)
    #         
    #     } else {
    #         sig <- sample_sig[s - 1]
    #         t1_old <- sample_t[[s - 1]]
    #     }
    #     
    #     
    #     #-----------------------------------------------
    #     # Sampling t
    #     #-----------------------------------------------
    #     t1_star <- stats::runif(1, min = a, max = b)
    #     
    #     # candi_den <- 1 / (b - a)
    #     
    #     log_post_den_res <- log_post_den_given_sig(t = t1_star, x = x, 
    #                                                theta = theta_k, 
    #                                                y = y, Kff = Kff, 
    #                                                sig = sig,
    #                                                shape1 = shape1, 
    #                                                shape2 = shape2, 
    #                                                a = a, b = b, 
    #                                                mixture_prob = mixture_prob)
    #     
    #     log_post_den <- log_post_den_res$log_post_den
    #     
    #     log_post_den_old_res <- log_post_den_given_sig(t = t1_old, x = x, 
    #                                                    theta = theta_k, 
    #                                                    y = y, Kff = Kff, 
    #                                                    sig = sig,
    #                                                    shape1 = shape1, 
    #                                                    shape2 = shape2, 
    #                                                    a = a, b = b, 
    #                                                    mixture_prob = mixture_prob)
    #     log_post_den_old <- log_post_den_old_res$log_post_den
    #     
    #     
    #     ## MH-step
    #     if (log(runif(1)) < (log_post_den - log_post_den_old)) {
    #         sample_t[[s]] <- t1_star
    #         sample_sig[s] <- sqrt(sample_sig2(y = y, x_vec = x,
    #                                           theta = theta_k, 
    #                                           Kff = Kff, 
    #                                           Kdf = log_post_den_res$Kdf,
    #                                           Kdd = log_post_den_res$Kdd,
    #                                           ga_shape = ga_shape, 
    #                                           ga_rate = ga_rate))
    #     } else {
    #         sample_t[[s]] <- t1_old
    #         sample_sig[s] <- sqrt(sample_sig2(y = y, x_vec = x,
    #                                           theta = theta_k, 
    #                                           Kff = Kff, 
    #                                           Kdf = log_post_den_old_res$Kdf,
    #                                           Kdd = log_post_den_old_res$Kdd,
    #                                           ga_shape = ga_shape, 
    #                                           ga_rate = ga_rate))
    #     }
    #     
    #     cat("final draw:", s, "\r")
    #     s <- s + 1
    # }
    
    
    cat("\n")
    print("Done!")
    colnames(thetas) <- c("tau", "h")
    
    return(list(sample_t = t(sample_t), 
                sample_sig = sample_sig, 
                thetas = thetas[1:count, ]))
}


mcem_dgp_matern <- function(y, x, H0, theta_init = c(1, 1), epsilon = 1e-4, 
                     D = 100, n_sample = 4000, a = 0, b = 2, 
                     max_iter = 100,
                     lower = c(0.001, 0.001),
                     upper = c(1/0.001, 1/0.001), init_t = 0.5, 
                     init_sig = 1,
                     shape1 = 1, shape2 = 1,
                     ctrl = list(TOL = 1e-5, trace = 0),
                     ga_shape = 1/2, ga_rate = 1/2,
                     a_h = 1, b_h = 1, nu = 1.5,
                     mixture_prob = c(0.5, 0.5),
                     is.h.par = FALSE) {
    ################################
    # Arguments
    
    ################################
    thetas <- matrix(0L, max_iter, length(theta_init))
    
    print("MCEM with sampling t and sig Begins")
    
    eps <- epsilon + 1
    count <- 1
    
    n_mixture <- length(shape1)
    
    # Initialization of parameters
    thetas[1, ] <- theta_k <- theta_init
    theta_new <- theta_k
    
    while(eps > epsilon) {
        
        if (count == 100) {
            print("Does not converge")
            break
        }
        
        sample_sig <- rep(init_sig, n_sample)
        
        Kff <- matern_ker(H0 = H0, nu = nu, tau = theta_k[1], l = theta_k[2])
        
        #==============================================
        # ********   Monte Carlo E-step  ************ #
        #==============================================
        
        if (count == 1) {
            sample_t <- matrix(stats::runif(n_sample, min = a, max = b), 
                               1, n_sample)
            
        } else {
            sample_t <- matrix(init_t, 1, n_sample)
            
            s <- 1
            print("Sampling begins")
            while (s <= n_sample) {
                
                if (s == 1) {
                    sig <- init_sig
                    t1_old <- runif(1, min = a, max = b)
                } else {
                    sig <- sample_sig[s - 1]
                    t1_old <- sample_t[[s - 1]]
                }
                
                
                #-----------------------------------------------
                # Sampling t
                #-----------------------------------------------
                t1_star <- stats::runif(1, min = a, max = b)
                
                # candi_den <- 1 / (b - a)
                
                log_post_den_res <- log_post_den_given_sig_matern(
                    t = t1_star, x = x, 
                    theta = theta_k, 
                    y = y, Kff = Kff, 
                    sig = sig,
                    nu = nu,
                    shape1 = shape1, 
                    shape2 = shape2, 
                    a = a, b = b, 
                    mixture_prob = mixture_prob
                    )
                
                log_post_den <- log_post_den_res$log_post_den
                
                log_post_den_old_res <- log_post_den_given_sig_matern(
                    t = t1_old, x = x, 
                    theta = theta_k, 
                    y = y, Kff = Kff, 
                    sig = sig,
                    nu = nu,
                    shape1 = shape1, 
                    shape2 = shape2, 
                    a = a, b = b, 
                    mixture_prob = mixture_prob
                    )
                log_post_den_old <- log_post_den_old_res$log_post_den
                
                ## MH-step
                if (log(runif(1)) < (log_post_den - log_post_den_old)) {
                    sample_t[[s]] <- t1_star
                    sample_sig[s] <- sqrt(sample_sig2(y = y, x_vec = x,
                                                      theta = theta_k, 
                                                      Kff = Kff, 
                                                      Kdf = log_post_den_res$Kdf,
                                                      Kdd = log_post_den_res$Kdd,
                                                      ga_shape = ga_shape, 
                                                      ga_rate = ga_rate))
                } else {
                    sample_t[[s]] <- t1_old
                    sample_sig[s] <- sqrt(sample_sig2(y = y, x_vec = x,
                                                      theta = theta_k, 
                                                      Kff = Kff, 
                                                      Kdf = log_post_den_old_res$Kdf,
                                                      Kdd = log_post_den_old_res$Kdd,
                                                      ga_shape = ga_shape, 
                                                      ga_rate = ga_rate))
                }
                cat("no. draw:", s, "\r")
                s <- s + 1
            }
        }
        
        
        #==============================================
        # ********   M-step for theta    ************ #
        #==============================================
        
        sample_t_mstep <- t(apply(sample_t, 1, sample, size = D))
        sample_sig_mstep <- sample(sample_sig, size = D)
        
        res <- Rsolnp::solnp(pars = theta_k, fun = log_mar_lik_gp_der_mc_sig2_matern,
                             LB = lower, control = ctrl, y = y, x_vec = x, 
                             H0 = H0, der_mc_mat = sample_t_mstep,
                             is.h.par = is.h.par, sample_sig = sample_sig_mstep,
                             a_h = a_h, b_h = b_h, nu = nu)
        theta_k <- res$par
        mar_post_k <- res$values[length(res$values)]
        
        
        # ******** Update epsilon and parameters ************ #
        
        eps <- sum((theta_new - theta_k) ^ 2)
        theta_new <- theta_k
        count <- count + 1
        print(paste("epsilon =", round(eps, 4), " count", count))
        print(paste("(increasing) marginal posterior =", -mar_post_k))
        thetas[count, ] <- theta_new
    }
    
    ## Final samples
    # sample_t <- matrix(0L, n_sample, 1)
    # sample_sig <- rep(mean(sample_sig), n_sample)
    # 
    # s <- 1
    # while (s <= n_sample) {
    #     
    #     if (s == 1) {
    #         sig <- sample_sig[s]
    #         t1_old <- runif(1, min = a, max = b)
    #         
    #     } else {
    #         sig <- sample_sig[s - 1]
    #         t1_old <- sample_t[[s - 1]]
    #     }
    #     
    #     
    #     #-----------------------------------------------
    #     # Sampling t
    #     #-----------------------------------------------
    #     t1_star <- stats::runif(1, min = a, max = b)
    #     
    #     # candi_den <- 1 / (b - a)
    #     
    #     log_post_den_res <- log_post_den_given_sig(t = t1_star, x = x, 
    #                                                theta = theta_k, 
    #                                                y = y, Kff = Kff, 
    #                                                sig = sig,
    #                                                shape1 = shape1, 
    #                                                shape2 = shape2, 
    #                                                a = a, b = b, 
    #                                                mixture_prob = mixture_prob)
    #     
    #     log_post_den <- log_post_den_res$log_post_den
    #     
    #     log_post_den_old_res <- log_post_den_given_sig(t = t1_old, x = x, 
    #                                                    theta = theta_k, 
    #                                                    y = y, Kff = Kff, 
    #                                                    sig = sig,
    #                                                    shape1 = shape1, 
    #                                                    shape2 = shape2, 
    #                                                    a = a, b = b, 
    #                                                    mixture_prob = mixture_prob)
    #     log_post_den_old <- log_post_den_old_res$log_post_den
    #     
    #     
    #     ## MH-step
    #     if (log(runif(1)) < (log_post_den - log_post_den_old)) {
    #         sample_t[[s]] <- t1_star
    #         sample_sig[s] <- sqrt(sample_sig2(y = y, x_vec = x,
    #                                           theta = theta_k, 
    #                                           Kff = Kff, 
    #                                           Kdf = log_post_den_res$Kdf,
    #                                           Kdd = log_post_den_res$Kdd,
    #                                           ga_shape = ga_shape, 
    #                                           ga_rate = ga_rate))
    #     } else {
    #         sample_t[[s]] <- t1_old
    #         sample_sig[s] <- sqrt(sample_sig2(y = y, x_vec = x,
    #                                           theta = theta_k, 
    #                                           Kff = Kff, 
    #                                           Kdf = log_post_den_old_res$Kdf,
    #                                           Kdd = log_post_den_old_res$Kdd,
    #                                           ga_shape = ga_shape, 
    #                                           ga_rate = ga_rate))
    #     }
    #     
    #     cat("final draw:", s, "\r")
    #     s <- s + 1
    # }
    
    
    cat("\n")
    print("Done!")
    colnames(thetas) <- c("tau", "h")
    
    return(list(sample_t = t(sample_t), 
                sample_sig = sample_sig, 
                thetas = thetas[1:count, ]))
}



sem_dgp_multi_sample_sig_mh <- function(y, x, H0, theta_init = c(1, 1), epsilon = 1e-4, 
                                  D = 100, a_vec = c(0, 1), b_vec = c(1, 2), 
                                  n_sample = 4000, 
                                  max_iter = 100,
                                  lower = c(0.001, 0.001),
                                  upper = c(1/0.001, 1/0.001), init_t = c(0.5, 0.5), 
                                  init_sig = 1,
                                  shape1 = 1, shape2 = 1,
                                  ctrl = list(TOL = 1e-5, trace = 0),
                                  ga_shape = 1/2, ga_rate = 1/2,
                                  mixture_prob = c(0.5, 0.5),
                                  a_h = 1, b_h = 1, is.h.par = FALSE) {
    ################################
    # Arguments
    ################################
    thetas <- matrix(0L, max_iter, length(theta_init))
    
    print("SEM (multiple t) with sampling t and sig Begins")
    
    eps <- epsilon + 1
    count <- 1
    
    nn <- length(init_t)  ## number of critical points
    
    # Initialization of parameters
    thetas[1, ] <- theta_k <- theta_init
    theta_new <- theta_k
    # log_prior_den <-  log(prod(1 / (b_vec - a_vec)))
    
    while(eps > epsilon) {
        
        if (count == 100) {
            print("Does not converge")
            break
        }
        
        sample_sig <- rep(init_sig, D)
        
        Kff <- se_ker(H0 = H0, tau = theta_k[1], h = theta_k[2])
        
        #==============================================
        # ******** Stochastic E-step ************ #
        #==============================================
        
        if (count == 1) {
            sample_t_mat <- matrix(stats::runif(D * nn, min = a_vec,
                                                max = b_vec), nn, D)
            
        } else {
            sample_t_mat <- matrix(init_t, nn, D)
            
            s <- 1
            print("Sampling begins")
            while (s <= D) {
                
                if (s == 1) {
                    sig <- init_sig
                    t_old <- runif(nn, min = a_vec, max = b_vec)
                    # print(t_old)
                } else {
                    sig <- sample_sig[s - 1]
                    t_old <- sample_t_mat[, s - 1]
                }
                
                
                #-----------------------------------------------
                # Sampling t
                #-----------------------------------------------
                t_star <- stats::runif(nn, min = a_vec, max = b_vec)
                
                for (i in 1:nn) {
                    t_new <- t_old
                    t_new[i] <- t_star[i]
                    log_post_den_res <- log_post_den_given_sig(t = t_new,
                                                               x = x, 
                                                               theta = theta_k, 
                                                               y = y, Kff = Kff, 
                                                               sig = sig,
                                                               shape1 = shape1, 
                                                               shape2 = shape2, 
                                                               a = a_vec, b = b_vec, 
                                                               mixture_prob = mixture_prob)
                    log_post_den <- log_post_den_res$log_post_den
                    
                    log_post_den_old_res <- log_post_den_given_sig(t = t_old, x = x, 
                                                                   theta = theta_k, 
                                                                   y = y, Kff = Kff, 
                                                                   sig = sig,
                                                                   shape1 = shape1, 
                                                                   shape2 = shape2, 
                                                                   a = a_vec, b = b_vec, 
                                                                   mixture_prob = mixture_prob)
                    log_post_den_old <- log_post_den_old_res$log_post_den
                    
                    
                    # print(t_old)
                    ## MH-step
                    if (log(runif(1)) < (log_post_den - log_post_den_old)) {
                        sample_t_mat[i, s] <- t_star[i]
                        t_old[i] <- t_star[i]
                    } else {
                        sample_t_mat[i, s] <- t_old[i]
                    }
                }
                
                #-----------------------------------------------
                # Sampling sig
                #-----------------------------------------------
                log_post_den_res <- log_post_den_given_sig(t = sample_t_mat[, s],
                                                           x = x, 
                                                           theta = theta_k, 
                                                           y = y, Kff = Kff, 
                                                           sig = sig,
                                                           shape1 = shape1, 
                                                           shape2 = shape2, 
                                                           a = a_vec, b = b_vec, 
                                                           mixture_prob = mixture_prob)
                
                sample_sig[s] <- sqrt(sample_sig2(y = y, x_vec = x,
                                                  theta = theta_k, 
                                                  Kff = Kff, 
                                                  Kdf = log_post_den_res$Kdf,
                                                  Kdd = log_post_den_res$Kdd,
                                                  ga_shape = ga_shape, 
                                                  ga_rate = ga_rate))

                cat("no. draw:", s, "\r")
                s <- s + 1
            }
        }
        
        
        #==============================================
        # ********   M-step for theta    ************ #
        #==============================================
        res <- Rsolnp::solnp(pars = theta_k, fun = log_mar_lik_gp_der_mc_sig2,
                             LB = lower, control = ctrl, y = y, x_vec = x, 
                             H0 = H0, der_mc_mat = sample_t_mat,
                             is.h.par = is.h.par, sample_sig = sample_sig,
                             a_h = a_h, b_h = b_h, is.multi = TRUE)
        theta_k <- res$par
        # mar_post_k <- res$values[length(res$values)]
        
        
        # ******** Update epsilon and parameters ************ #
        
        eps <- sum((theta_new - theta_k) ^ 2)
        theta_new <- theta_k
        count <- count + 1
        print(paste("epsilon =", round(eps, 4), " count", count))
        # print(paste("(increasing) marginal posterior =", -mar_post_k))
        thetas[count, ] <- theta_new
    }
    
    ## Final samples
    sample_t_mat <- matrix(0L, nrow = n_sample, ncol = nn)
    sample_sig <- rep(mean(sample_sig), n_sample)
    
    s <- 1
    while (s <= n_sample) {
        
        if (s == 1) {
            sig <- init_sig
            t1_old <- runif(nn, min = a_vec, max = b_vec)
        } else {
            sig <- sample_sig[s - 1]
            t_old <- sample_t_mat[s - 1, ]
        }
        
        
        #-----------------------------------------------
        # Sampling t
        #-----------------------------------------------
        t_star <- stats::runif(nn, min = a_vec, max = b_vec)
        
        for (i in 1:nn) {
            t_new <- t_old
            t_new[i] <- t_star[i]
            log_post_den_res <- log_post_den_given_sig(t = t_new,
                                                       x = x, 
                                                       theta = theta_k, 
                                                       y = y, Kff = Kff, 
                                                       sig = sig,
                                                       shape1 = shape1, 
                                                       shape2 = shape2, 
                                                       a = a_vec, b = b_vec, 
                                                       mixture_prob = mixture_prob)
            log_post_den <- log_post_den_res$log_post_den
            
            log_post_den_old_res <- log_post_den_given_sig(t = t_old, x = x, 
                                                           theta = theta_k, 
                                                           y = y, Kff = Kff, 
                                                           sig = sig,
                                                           shape1 = shape1, 
                                                           shape2 = shape2, 
                                                           a = a_vec, b = b_vec, 
                                                           mixture_prob = mixture_prob)
            log_post_den_old <- log_post_den_old_res$log_post_den
            
            
            
            ## MH-step
            if (log(runif(1)) < (log_post_den - log_post_den_old)) {
                sample_t_mat[s, i] <- t_star[i]
                t_old[i] <- t_star[i]
            } else {
                sample_t_mat[s, i] <- t_old[i]
            }
        }
        
        #-----------------------------------------------
        # Sampling sig
        #-----------------------------------------------
        log_post_den_res <- log_post_den_given_sig(t = sample_t_mat[s, ],
                                                   x = x, 
                                                   theta = theta_k, 
                                                   y = y, Kff = Kff, 
                                                   sig = sig,
                                                   shape1 = shape1, 
                                                   shape2 = shape2, 
                                                   a = a_vec, b = b_vec, 
                                                   mixture_prob = mixture_prob)
        
        sample_sig[s] <- sqrt(sample_sig2(y = y, x_vec = x,
                                          theta = theta_k, 
                                          Kff = Kff, 
                                          Kdf = log_post_den_res$Kdf,
                                          Kdd = log_post_den_res$Kdd,
                                          ga_shape = ga_shape, 
                                          ga_rate = ga_rate))
        
        cat("final draw:", s, "\r")
        s <- s + 1
    }
    
    
    cat("\n")
    print("Done!")
    colnames(thetas) <- c("tau", "h")
    
    return(list(sample_t_mat = sample_t_mat, 
                sample_sig = sample_sig, 
                thetas = thetas[1:count, ]))
}


mcem_dgp_multi <- function(y, x, H0, theta_init = c(1, 1), epsilon = 1e-4, 
                           D = 100, a_vec = c(0, 1), b_vec = c(1, 2), 
                           n_sample = 4000, 
                           max_iter = 100,
                           lower = c(0.001, 0.001),
                           upper = c(1/0.001, 1/0.001), init_t = c(0.5, 0.5), 
                           init_sig = 1,
                           shape1 = 1, shape2 = 1,
                           ctrl = list(TOL = 1e-5, trace = 0),
                           ga_shape = 1/2, ga_rate = 1/2, 
                           mixture_prob = c(0.5, 0.5),
                           a_h = 1, b_h = 1, is.h.par = FALSE) {
    ################################
    # Arguments
    ################################
    thetas <- matrix(0L, max_iter, length(theta_init))
    
    print("MCEM (multiple t) with sampling t and sig Begins")
    
    eps <- epsilon + 1
    count <- 1
    
    nn <- length(init_t)  ## number of critical points
    
    # Initialization of parameters
    thetas[1, ] <- theta_k <- theta_init
    theta_new <- theta_k
    # log_prior_den <-  log(prod(1 / (b_vec - a_vec)))
    
    while(eps > epsilon) {
        
        if (count == 100) {
            print("Does not converge")
            break
        }
        
        sample_sig <- rep(init_sig, n_sample)
        
        Kff <- se_ker(H0 = H0, tau = theta_k[1], h = theta_k[2])
        
        #==============================================
        # ******** Stochastic E-step ************ #
        #==============================================
        
        if (count == 1) {
            sample_t_mat <- matrix(stats::runif(D * nn, min = a_vec,
                                                max = b_vec), nn, n_sample)
            
        } else {
            sample_t_mat <- matrix(init_t, nn, n_sample)
            
            s <- 1
            print("Sampling begins")
            while (s <= n_sample) {
                
                if (s == 1) {
                    sig <- init_sig
                    t_old <- runif(nn, min = a_vec, max = b_vec)
                    # print(t_old)
                } else {
                    sig <- sample_sig[s - 1]
                    t_old <- sample_t_mat[, s - 1]
                }
                
                
                #-----------------------------------------------
                # Sampling t
                #-----------------------------------------------
                t_star <- stats::runif(nn, min = a_vec, max = b_vec)
                
                for (i in 1:nn) {
                    t_new <- t_old
                    t_new[i] <- t_star[i]
                    log_post_den_res <- log_post_den_given_sig(t = t_new,
                                                               x = x,
                                                               theta = theta_k,
                                                               y = y, Kff = Kff,
                                                               sig = sig,
                                                               shape1 = shape1,
                                                               shape2 = shape2,
                                                               a = a_vec, b = b_vec,
                                                               mixture_prob = mixture_prob)
                    log_post_den <- log_post_den_res$log_post_den

                    log_post_den_old_res <- log_post_den_given_sig(t = t_old, x = x,
                                                                   theta = theta_k,
                                                                   y = y, Kff = Kff,
                                                                   sig = sig,
                                                                   shape1 = shape1,
                                                                   shape2 = shape2,
                                                                   a = a_vec, b = b_vec,
                                                                   mixture_prob = mixture_prob)
                    log_post_den_old <- log_post_den_old_res$log_post_den


                    # print(t_old)
                    ## MH-step
                    if (log(runif(1)) < (log_post_den - log_post_den_old)) {
                        sample_t_mat[i, s] <- t_star[i]
                        t_old[i] <- t_star[i]
                    } else {
                        sample_t_mat[i, s] <- t_old[i]
                    }
                }
                
                
                
                # t_new <- t_old
                # t_star <- stats::runif(nn, min = a_vec, max = b_vec)
                # # t_new[i] <- t_star[i]
                # log_post_den_res <- log_post_den_given_sig(t = t_star,
                #                                            x = x,
                #                                            theta = theta_k,
                #                                            y = y, Kff = Kff,
                #                                            sig = sig,
                #                                            shape1 = shape1,
                #                                            shape2 = shape2,
                #                                            a = a_vec, b = b_vec,
                #                                            mixture_prob = mixture_prob)
                # log_post_den <- log_post_den_res$log_post_den
                # 
                # log_post_den_old_res <- log_post_den_given_sig(t = t_old, x = x,
                #                                                theta = theta_k,
                #                                                y = y, Kff = Kff,
                #                                                sig = sig,
                #                                                shape1 = shape1,
                #                                                shape2 = shape2,
                #                                                a = a_vec, b = b_vec,
                #                                                mixture_prob = mixture_prob)
                # log_post_den_old <- log_post_den_old_res$log_post_den
                # 
                # 
                # # print(t_old)
                # ## MH-step
                # if (log(runif(1)) < (log_post_den - log_post_den_old)) {
                #     sample_t_mat[, s] <- t_star
                #     # t_old <- t_star
                # } else {
                #     sample_t_mat[, s] <- t_old
                # }
                
                #-----------------------------------------------
                # Sampling sig
                #-----------------------------------------------
                log_post_den_res <- log_post_den_given_sig(t = sample_t_mat[, s],
                                                           x = x, 
                                                           theta = theta_k, 
                                                           y = y, Kff = Kff, 
                                                           sig = sig,
                                                           shape1 = shape1, 
                                                           shape2 = shape2, 
                                                           a = a_vec, b = b_vec, 
                                                           mixture_prob = mixture_prob)
                
                sample_sig[s] <- sqrt(sample_sig2(y = y, x_vec = x,
                                                  theta = theta_k, 
                                                  Kff = Kff, 
                                                  Kdf = log_post_den_res$Kdf,
                                                  Kdd = log_post_den_res$Kdd,
                                                  ga_shape = ga_shape, 
                                                  ga_rate = ga_rate))
                
                cat("no. draw:", s, "\r")
                s <- s + 1
            }
        }
        
        
        #==============================================
        # ********   M-step for theta    ************ #
        #==============================================
        
        sample_t_mstep <- t(apply(sample_t_mat, 1, sample, size = D))
        sample_sig_mstep <- sample(sample_sig, size = D)
        
        
        res <- Rsolnp::solnp(pars = theta_k, fun = log_mar_lik_gp_der_mc_sig2,
                             LB = lower, control = ctrl, y = y, x_vec = x, 
                             H0 = H0, der_mc_mat = sample_t_mstep,
                             is.h.par = is.h.par, sample_sig = sample_sig_mstep,
                             a_h = a_h, b_h = b_h, is.multi = TRUE)
        theta_k <- res$par
        # mar_post_k <- res$values[length(res$values)]
        
        
        # ******** Update epsilon and parameters ************ #
        
        eps <- sum((theta_new - theta_k) ^ 2)
        theta_new <- theta_k
        count <- count + 1
        print(paste("epsilon =", round(eps, 4), " count", count))
        # print(paste("(increasing) marginal posterior =", -mar_post_k))
        thetas[count, ] <- theta_new
    }
    
    # ## Final samples
    # sample_t_mat <- matrix(0L, nrow = n_sample, ncol = nn)
    # sample_sig <- rep(mean(sample_sig), n_sample)
    # 
    # s <- 1
    # while (s <= n_sample) {
    #     
    #     if (s == 1) {
    #         sig <- init_sig
    #         t1_old <- runif(nn, min = a_vec, max = b_vec)
    #     } else {
    #         sig <- sample_sig[s - 1]
    #         t_old <- sample_t_mat[s - 1, ]
    #     }
    #     
    #     
    #     #-----------------------------------------------
    #     # Sampling t
    #     #-----------------------------------------------
    #     t_star <- stats::runif(nn, min = a_vec, max = b_vec)
    #     
    #     for (i in 1:nn) {
    #         t_new <- t_old
    #         t_new[i] <- t_star[i]
    #         log_post_den_res <- log_post_den_given_sig(t = t_new,
    #                                                    x = x, 
    #                                                    theta = theta_k, 
    #                                                    y = y, Kff = Kff, 
    #                                                    sig = sig,
    #                                                    shape1 = shape1, 
    #                                                    shape2 = shape2, 
    #                                                    a = a_vec, b = b_vec, 
    #                                                    mixture_prob = mixture_prob)
    #         log_post_den <- log_post_den_res$log_post_den
    #         
    #         log_post_den_old_res <- log_post_den_given_sig(t = t_old, x = x, 
    #                                                        theta = theta_k, 
    #                                                        y = y, Kff = Kff, 
    #                                                        sig = sig,
    #                                                        shape1 = shape1, 
    #                                                        shape2 = shape2, 
    #                                                        a = a_vec, b = b_vec, 
    #                                                        mixture_prob = mixture_prob)
    #         log_post_den_old <- log_post_den_old_res$log_post_den
    #         
    #         
    #         
    #         ## MH-step
    #         if (log(runif(1)) < (log_post_den - log_post_den_old)) {
    #             sample_t_mat[s, i] <- t_star[i]
    #             t_old[i] <- t_star[i]
    #         } else {
    #             sample_t_mat[s, i] <- t_old[i]
    #         }
    #     }
    #     
    #     #-----------------------------------------------
    #     # Sampling sig
    #     #-----------------------------------------------
    #     log_post_den_res <- log_post_den_given_sig(t = sample_t_mat[s, ],
    #                                                x = x, 
    #                                                theta = theta_k, 
    #                                                y = y, Kff = Kff, 
    #                                                sig = sig,
    #                                                shape1 = shape1, 
    #                                                shape2 = shape2, 
    #                                                a = a_vec, b = b_vec, 
    #                                                mixture_prob = mixture_prob)
    #     
    #     sample_sig[s] <- sqrt(sample_sig2(y = y, x_vec = x,
    #                                       theta = theta_k, 
    #                                       Kff = Kff, 
    #                                       Kdf = log_post_den_res$Kdf,
    #                                       Kdd = log_post_den_res$Kdd,
    #                                       ga_shape = ga_shape, 
    #                                       ga_rate = ga_rate))
    #     
    #     cat("final draw:", s, "\r")
    #     s <- s + 1
    # }
    
    
    cat("\n")
    print("Done!")
    colnames(thetas) <- c("tau", "h")
    
    return(list(sample_t_mat = sample_t_mat, 
                sample_sig = sample_sig, 
                thetas = thetas[1:count, ]))
}


mcem_dgp_multi_matern <- function(y, x, H0, theta_init = c(1, 1), epsilon = 1e-4, 
                           D = 100, a_vec = c(0, 1), b_vec = c(1, 2), 
                           n_sample = 4000, 
                           max_iter = 100,
                           lower = c(0.001, 0.001),
                           upper = c(1/0.001, 1/0.001), init_t = c(0.5, 0.5), 
                           init_sig = 1,
                           shape1 = 1, shape2 = 1,
                           ctrl = list(TOL = 1e-5, trace = 0),
                           ga_shape = 1/2, ga_rate = 1/2, 
                           mixture_prob = c(0.5, 0.5),
                           a_h = 1, b_h = 1, nu = 1.5,
                           is.h.par = FALSE) {
    ################################
    # Arguments
    ################################
    thetas <- matrix(0L, max_iter, length(theta_init))
    
    print("MCEM (multiple t) with sampling t and sig Begins")
    
    eps <- epsilon + 1
    count <- 1
    
    nn <- length(init_t)  ## number of critical points
    
    # Initialization of parameters
    thetas[1, ] <- theta_k <- theta_init
    theta_new <- theta_k
    # log_prior_den <-  log(prod(1 / (b_vec - a_vec)))
    
    while(eps > epsilon) {
        
        if (count == 100) {
            print("Does not converge")
            break
        }
        
        sample_sig <- rep(init_sig, n_sample)
        
        # Kff <- se_ker(H0 = H0, tau = theta_k[1], h = theta_k[2])
        Kff <- matern_ker(H0 = H0, nu = nu, tau = theta_k[1], l = theta_k[2])
        #==============================================
        # ******** Stochastic E-step ************ #
        #==============================================
        
        if (count == 1) {
            sample_t_mat <- matrix(stats::runif(D * nn, min = a_vec,
                                                max = b_vec), nn, n_sample)
            
        } else {
            sample_t_mat <- matrix(init_t, nn, n_sample)
            
            s <- 1
            print("Sampling begins")
            while (s <= n_sample) {
                
                if (s == 1) {
                    sig <- init_sig
                    t_old <- runif(nn, min = a_vec, max = b_vec)
                    # print(t_old)
                } else {
                    sig <- sample_sig[s - 1]
                    t_old <- sample_t_mat[, s - 1]
                }
                
                
                #-----------------------------------------------
                # Sampling t
                #-----------------------------------------------
                t_star <- stats::runif(nn, min = a_vec, max = b_vec)
                
                for (i in 1:nn) {
                    t_new <- t_old
                    t_new[i] <- t_star[i]
                    log_post_den_res <- log_post_den_given_sig_matern(
                        t = t_new,
                        x = x,
                        theta = theta_k,
                        y = y, Kff = Kff,
                        sig = sig,
                        nu = nu,
                        shape1 = shape1,
                        shape2 = shape2,
                        a = a_vec, b = b_vec,
                        mixture_prob = mixture_prob
                        )
                    log_post_den <- log_post_den_res$log_post_den
                    
                    log_post_den_old_res <- log_post_den_given_sig_matern(
                        t = t_old, x = x,
                        theta = theta_k,
                        y = y, Kff = Kff,
                        sig = sig,
                        nu = nu,
                        shape1 = shape1,
                        shape2 = shape2,
                        a = a_vec, b = b_vec,
                        mixture_prob = mixture_prob
                        )
                    log_post_den_old <- log_post_den_old_res$log_post_den
                    
                    
                    # print(t_old)
                    ## MH-step
                    if (log(runif(1)) < (log_post_den - log_post_den_old)) {
                        sample_t_mat[i, s] <- t_star[i]
                        t_old[i] <- t_star[i]
                    } else {
                        sample_t_mat[i, s] <- t_old[i]
                    }
                }
                
                
                
                # t_new <- t_old
                # t_star <- stats::runif(nn, min = a_vec, max = b_vec)
                # # t_new[i] <- t_star[i]
                # log_post_den_res <- log_post_den_given_sig(t = t_star,
                #                                            x = x,
                #                                            theta = theta_k,
                #                                            y = y, Kff = Kff,
                #                                            sig = sig,
                #                                            shape1 = shape1,
                #                                            shape2 = shape2,
                #                                            a = a_vec, b = b_vec,
                #                                            mixture_prob = mixture_prob)
                # log_post_den <- log_post_den_res$log_post_den
                # 
                # log_post_den_old_res <- log_post_den_given_sig(t = t_old, x = x,
                #                                                theta = theta_k,
                #                                                y = y, Kff = Kff,
                #                                                sig = sig,
                #                                                shape1 = shape1,
                #                                                shape2 = shape2,
                #                                                a = a_vec, b = b_vec,
                #                                                mixture_prob = mixture_prob)
                # log_post_den_old <- log_post_den_old_res$log_post_den
                # 
                # 
                # # print(t_old)
                # ## MH-step
                # if (log(runif(1)) < (log_post_den - log_post_den_old)) {
                #     sample_t_mat[, s] <- t_star
                #     # t_old <- t_star
                # } else {
                #     sample_t_mat[, s] <- t_old
                # }
                
                #-----------------------------------------------
                # Sampling sig
                #-----------------------------------------------
                log_post_den_res <- log_post_den_given_sig_matern(
                    t = sample_t_mat[, s],
                    x = x, 
                    theta = theta_k, 
                    y = y, Kff = Kff, 
                    sig = sig,
                    nu = nu,
                    shape1 = shape1, 
                    shape2 = shape2, 
                    a = a_vec, b = b_vec, 
                    mixture_prob = mixture_prob)
                
                sample_sig[s] <- sqrt(sample_sig2(y = y, x_vec = x,
                                                  theta = theta_k, 
                                                  Kff = Kff, 
                                                  Kdf = log_post_den_res$Kdf,
                                                  Kdd = log_post_den_res$Kdd,
                                                  ga_shape = ga_shape, 
                                                  ga_rate = ga_rate))
                
                cat("no. draw:", s, "\r")
                s <- s + 1
            }
        }
        
        
        #==============================================
        # ********   M-step for theta    ************ #
        #==============================================
        
        sample_t_mstep <- t(apply(sample_t_mat, 1, sample, size = D))
        sample_sig_mstep <- sample(sample_sig, size = D)
        
        
        res <- Rsolnp::solnp(pars = theta_k, fun = log_mar_lik_gp_der_mc_sig2_matern,
                             LB = lower, control = ctrl, y = y, x_vec = x, 
                             H0 = H0, der_mc_mat = sample_t_mstep,
                             is.h.par = is.h.par, sample_sig = sample_sig_mstep,
                             a_h = a_h, b_h = b_h, is.multi = TRUE, nu = nu)
        theta_k <- res$par
        # mar_post_k <- res$values[length(res$values)]
        
        
        # ******** Update epsilon and parameters ************ #
        
        eps <- sum((theta_new - theta_k) ^ 2)
        theta_new <- theta_k
        count <- count + 1
        print(paste("epsilon =", round(eps, 4), " count", count))
        # print(paste("(increasing) marginal posterior =", -mar_post_k))
        thetas[count, ] <- theta_new
    }
    
    # ## Final samples
    # sample_t_mat <- matrix(0L, nrow = n_sample, ncol = nn)
    # sample_sig <- rep(mean(sample_sig), n_sample)
    # 
    # s <- 1
    # while (s <= n_sample) {
    #     
    #     if (s == 1) {
    #         sig <- init_sig
    #         t1_old <- runif(nn, min = a_vec, max = b_vec)
    #     } else {
    #         sig <- sample_sig[s - 1]
    #         t_old <- sample_t_mat[s - 1, ]
    #     }
    #     
    #     
    #     #-----------------------------------------------
    #     # Sampling t
    #     #-----------------------------------------------
    #     t_star <- stats::runif(nn, min = a_vec, max = b_vec)
    #     
    #     for (i in 1:nn) {
    #         t_new <- t_old
    #         t_new[i] <- t_star[i]
    #         log_post_den_res <- log_post_den_given_sig(t = t_new,
    #                                                    x = x, 
    #                                                    theta = theta_k, 
    #                                                    y = y, Kff = Kff, 
    #                                                    sig = sig,
    #                                                    shape1 = shape1, 
    #                                                    shape2 = shape2, 
    #                                                    a = a_vec, b = b_vec, 
    #                                                    mixture_prob = mixture_prob)
    #         log_post_den <- log_post_den_res$log_post_den
    #         
    #         log_post_den_old_res <- log_post_den_given_sig(t = t_old, x = x, 
    #                                                        theta = theta_k, 
    #                                                        y = y, Kff = Kff, 
    #                                                        sig = sig,
    #                                                        shape1 = shape1, 
    #                                                        shape2 = shape2, 
    #                                                        a = a_vec, b = b_vec, 
    #                                                        mixture_prob = mixture_prob)
    #         log_post_den_old <- log_post_den_old_res$log_post_den
    #         
    #         
    #         
    #         ## MH-step
    #         if (log(runif(1)) < (log_post_den - log_post_den_old)) {
    #             sample_t_mat[s, i] <- t_star[i]
    #             t_old[i] <- t_star[i]
    #         } else {
    #             sample_t_mat[s, i] <- t_old[i]
    #         }
    #     }
    #     
    #     #-----------------------------------------------
    #     # Sampling sig
    #     #-----------------------------------------------
    #     log_post_den_res <- log_post_den_given_sig(t = sample_t_mat[s, ],
    #                                                x = x, 
    #                                                theta = theta_k, 
    #                                                y = y, Kff = Kff, 
    #                                                sig = sig,
    #                                                shape1 = shape1, 
    #                                                shape2 = shape2, 
    #                                                a = a_vec, b = b_vec, 
    #                                                mixture_prob = mixture_prob)
    #     
    #     sample_sig[s] <- sqrt(sample_sig2(y = y, x_vec = x,
    #                                       theta = theta_k, 
    #                                       Kff = Kff, 
    #                                       Kdf = log_post_den_res$Kdf,
    #                                       Kdd = log_post_den_res$Kdd,
    #                                       ga_shape = ga_shape, 
    #                                       ga_rate = ga_rate))
    #     
    #     cat("final draw:", s, "\r")
    #     s <- s + 1
    # }
    
    
    cat("\n")
    print("Done!")
    colnames(thetas) <- c("tau", "h")
    
    return(list(sample_t_mat = sample_t_mat, 
                sample_sig = sample_sig, 
                thetas = thetas[1:count, ]))
}






metropolis_log_ratio_ts_subj <- function(ys, ts, ws = NULL, 
                                         ts_new, ws_new = NULL, 
                                         Sigma_s, Sigma_s_new,
                                         sig, s, shape1, shape2, a, b,
                                         proposal_type = NULL) {
    # conditional of ts (Metropolis ratio)
    # ============================================
    n <- length(ys)
    # n_subj <- length(t_vec)
    # sig_diag <- diag(theta[1] ^ 2, n)
    sig_diag <- diag(sig ^ 2, n)
    mu <- rep(0L, n)
    
    log_normal <- dmvn(ys, mu = mu, sigma = sig_diag + Sigma_s, log = TRUE)
    
    log_normal_new <- dmvn(ys, mu = mu, sigma = sig_diag + Sigma_s_new, 
                           log = TRUE)
    
    log_beta_new <- log_general_dbeta(ts_new, shape1 = shape1,
                                      shape2 = shape2,
                                      a = a, b = b)
    log_beta <- log_general_dbeta(ts, shape1 = shape1,
                                  shape2 = shape2,
                                  a = a, b = b)
    if (proposal_type == "transform_normal") {
        log_normal <- log_normal + ws + log(b - a) - 
            2 * log(1 + exp(ws))
        log_normal_new <- log_normal_new + ws_new + log(b - a) - 
            2 * log(1 + exp(ws_new))
    } 
    
    return(log_normal_new + log_beta_new - log_normal - log_beta)
}



mcmc_Estep_subj <- function(multi_y, x_vec, n_subj, theta_k, H0,
                            start.lst, n_mcmc = 1000, burn = 100, thin = 1, 
                            name.par, adapt = TRUE, tune.lst, keep.lst,
                            tune.len = 50, 
                            target.accept.rate = 0.35, 
                            proposal_type = "uniform", a = NULL, b = NULL) {
    
    
    # proposal_type = "normal", "uniform", "transform_normal"
    
    #-------------------------------
    # Functions used in the algorithm
    #-------------------------------
    get.tune <- function(tune, keep, k, target = target.accept.rate){  
        # adaptive tuning
        a <- min(0.025, 1 / sqrt(k))
        exp(ifelse(keep < target, log(tune) - a, log(tune) + a))
    }
    
    
    #-------------------------------
    # Adaptive tuning
    #------------------------------- 
    Tb <- tune.len  # frequency of adaptive tuning
    # for t
    keep <- keep.lst
    keep.tmp <- keep  # track MH accpetance rate for adaptive tuning
    
    #-------------------------------
    # Here create some objects that are used later.
    #-------------------------------
    # sig <- theta_k[1]
    tau <- theta_k[1]
    h <- theta_k[2]
    
    Kff <- se_ker(H0 = H0, tau = tau, h = h)
    
    
    #-------------------------------
    # Storage
    #-------------------------------
    # name.par should include t and sig
    sampleidx <- seq(from = (burn + thin), to = n_mcmc, by = thin)
    draws <- matrix(NA, nrow = length(sampleidx), ncol = length(name.par))
    colnames(draws) <- name.par
    
    
    #-------------------------------
    # Starting values
    #-------------------------------
    t_vec <- start.lst$t
    sig <- start.lst$sig
    #-------------------------------
    # M-H algorithm
    #------------------------------- 
    for (i in 1:n_mcmc) {
        if (i %% 5 == 0) cat("iter:", i, "\r")
        flush.console()
        
        #------------------------------------------------------------------
        # Sampling posterior distribution of t given theta^(i)
        #------------------------------------------------------------------
        Sigma_lst <- lapply(1:n_subj, function(s) {
            Kdf <- computeCovDer1(idx1 = t_vec[s], idx2 = x_vec, tau = tau, 
                                  h = h) 
            
            Kdd <- computeCovDer2(idx1 = t_vec[s], tau = tau, h = h)
            Sigma <- Kff - quad.form.inv(Kdd, Kdf)
            Sigma <- (Sigma + t(Sigma))/2
            # if(!is.positive.definite(Sigma)) {
            #     Sigma <- as.matrix(nearPD(Sigma)$mat)
            # }
            return(Sigma)
        })
        Sigma_lst_update <- Sigma_lst
        
        for (s in 1:n_subj) {
            # Update tuning parameter
            #------------------------
            if(adapt == TRUE & i %% Tb == 0) {
                # Adaptive tuning
                keep.tmp[[s]] <- keep.tmp[[s]] / Tb
                tune.lst[[s]] <- get.tune(tune.lst[[s]], keep.tmp[[s]], i)
                keep.tmp[[s]] <- 0
            }
            
            ts_star <- stats::runif(1, min = a, max = b)
            
            Kdf_new <- computeCovDer1(idx1 = ts_star, idx2 = x_vec, tau = tau, 
                                      h = h) 
            
            Kdd_new <- computeCovDer2(idx1 = ts_star, tau = tau, h = h)
            
            Sigma_new <- Kff - quad.form.inv(Kdd_new, Kdf_new)
            
            log_rho <- metropolis_log_ratio_ts_subj(ys = multi_y[, s],
                                                    ts = t_vec[s], 
                                                    ts_new = ts_star,
                                                    Sigma_s = Sigma_lst_update[[s]], 
                                                    Sigma_s_new = Sigma_new,
                                                    sig = sig, s = s,
                                                    shape1 = shape1, 
                                                    shape2 = shape2, 
                                                    a = a, b = b,
                                                    proposal_type = proposal_type)
            if (log(runif(1)) < log_rho) {
                t_vec[s] <- ts_star
                keep[[s]] <- keep[[s]] + 1
                keep.tmp[[s]] <- keep.tmp[[s]] + 1
                Sigma_lst_update[[s]] <- Sigma_new
            }
            
            
        }
        
        sig <- sqrt(sample_sig2_subj(multi_y = multi_y, x_vec = x,
                                     theta = theta_k, 
                                     Sigma = Sigma_lst_update,
                                     ga_shape = ga_shape, 
                                     ga_rate = ga_rate))
        
        #  Save samples of t and sig
        # -----------------------------------------------------
        if (i > burn) {
            if (i %% thin == 0) {
                draws[(i - burn) %/% thin, ] <- c(t_vec, sig)
            } 
        }
        

    }
    
    # Acceptance Probability
    #----------------------------
    keep <- lapply(keep, function(x) x / n_mcmc)
    # Write output
    #--------------
    return(list(draws = draws, 
                accept = keep,
                start = start.lst, 
                tune = tune.lst, 
                burn = burn, 
                thin = thin, 
                n_mcmc = n_mcmc, 
                sampleidx = sampleidx))
}


sem_dgp_subj_sample_sig_mh <- function(multi_y, x, H0, theta_init = c(1, 1), 
                                       epsilon = 1e-4, D = 100, a = 0, b = 2, 
                                       max_iter = 100,
                                       lower = c(0.001, 0.001, 0.001),
                                       upper = c(1/0.001, 1/0.001, 1/0.001), 
                                       shape1 = 1, shape2 = 1,
                                       ctrl = list(TOL = 1e-5, trace = 0),
                                       ga_shape = 1/2, ga_rate = 1/2,
                                       n_mcmc = 12000, start.lst, 
                                       burn_e = 2000, burn_final = 2000, 
                                       thin_e = 100, thin_final = 1, 
                                       name.par, adapt = TRUE, 
                                       tune.lst, keep.lst, tune.len = 20, 
                                       target.accept.rate = 0.35,
                                       proposal_type = "uniform",
                                       is.h.par = FALSE, 
                                       a_h = 10, b_h = 10) {
    ################################
    # Arguments
    
    ################################
    eps <- epsilon + 1
    count <- 1
    n_subj <- ncol(multi_y)
    # print(paste("epsilon =", round(eps, decimalplaces(epsilon) * 2), 
    #             " count", count))
    
    # Initialization of parameters
    n_theta <- length(theta_init)
    theta_mat <- matrix(0L, max_iter, n_theta)
    theta_mat[1, ] <- theta_k <- theta_init
    theta_new <- theta_k
    mar_post_new <- 0
    
    # n_mixture <- length(shape1)
    
    if (shape1 == 1 && shape2 == 1) {
        log_prior_den <- log_general_dbeta((a + b) / 2, shape1 = shape1,
                                           shape2 = shape2, 
                                           a = a, b = b)
    }
    
    print("MCEM w/ sampling t and sig (multi-subj mcmc) Begins")
    while(eps > epsilon) {
        
        if (count == 100) {
            print("Does not converge")
            break
        }
        
        # *****************************************
        # ******** Stochastic E-step ************ #
        #  ****************************************
        ## MCMC for the Stochastic E-step
        
        #-----------------------------------------------
        # Sampling t given sig, theta = (tau, h), y
        #-----------------------------------------------
        print("E-step MCMC sampling for t and sig begins")
        if (count == 1) {
            start.lst <- start.lst
        } else{
            post_mean <- apply(mcmc_out$draws, 2, mean)
            start.lst <- list(t = post_mean[1:n_subj],
                              sig = post_mean[n_subj + 1])
        }
        mcmc_out <- mcmc_Estep_subj(multi_y = multi_y, x_vec = x, 
                                    n_subj = n_subj, 
                                    theta_k = theta_k, H0 = H0, 
                                    start.lst = start.lst,
                                    n_mcmc = n_mcmc, 
                                    burn = burn_e, 
                                    thin = thin_e, 
                                    name.par = name.par, adapt = adapt, 
                                    tune.lst = tune.lst, 
                                    keep.lst = keep.lst,
                                    tune.len = tune.len, 
                                    target.accept.rate = target.accept.rate,
                                    proposal_type = proposal_type, 
                                    a = a, b = b)
        
        len <- length(mcmc_out$sampleidx)
        sample_t <- mcmc_out$draws[, 1:n_subj]
        sample_sig <- mcmc_out$draws[, (n_subj + 1)]
        
        sample_t_mstep <- apply(sample_t, 2, sample, size = D)
        sample_sig_mstep <- sample(sample_sig, size = D)
        
        
        # *****************************************
        # ******** M-step for theta ************* #
        # *****************************************
        print("M-step updating kernel parameters")
        
        res <- Rsolnp::solnp(pars = theta_k, fun = marg_lik_gp_der_mc_sig_subj,
                             LB = lower, control = ctrl, multi_y = multi_y, 
                             x_vec = x,
                             H0 = H0, der_mc_mat = sample_t_mstep,
                             sample_sig = sample_sig_mstep,
                             a_h = a_h, b_h = b_h, is.h.par = is.h.par)
        theta_k <- res$par
        mar_post_k <- res$values[length(res$values)]
        print(paste("theta =", theta_k))
        print(paste("(increasing) marginal posterior =", -mar_post_k))
        
        # *****************************************************
        # ******** Update epsilon and parameters ************ #
        # *****************************************************
        eps <- sum((theta_new - theta_k) ^ 2)
        theta_new <- theta_k
        mar_post_new <- mar_post_k
        count <- count + 1
        print(paste("epsilon =", round(eps, 4), " count", count))
        theta_mat[count, ] <- theta_new
    }
    
    ## =========================================================================
    ## Final MCMC samples for t
    ## =========================================================================
    # print("Final MCMC sampling for t and sigma")
    # post_mean <- apply(mcmc_out$draws, 2, mean)
    # start.lst <- list(t = post_mean[1:n_subj],
    #                   sig = post_mean[n_subj + 1])
    # mcmc_final <- mcmc_Estep_subj(multi_y = multi_y, x_vec = x, 
    #                               n_subj = n_subj,
    #                               theta_k = theta_k, H0 = H0, 
    #                               start.lst = start.lst,
    #                               n_mcmc = n_mcmc, 
    #                               burn = burn_e, 
    #                               thin = thin_e, 
    #                               name.par = name.par, adapt = adapt, 
    #                               tune.lst = tune.lst, 
    #                               keep.lst = keep.lst,
    #                               tune.len = tune.len, 
    #                               target.accept.rate = target.accept.rate,
    #                               proposal_type = proposal_type, 
    #                               a = a, b = b)
    cat("\n")
    print("Done!")
    colnames(theta_mat) <- c("tau", "h")
    
    return(list(mcmc_output = mcmc_out, 
                theta_mat = theta_mat[1:count, ]))
}



marg_lik_gp_der_mc_sig_subj <- function(theta, multi_y, x_vec, der_mc_mat, H0,
                                          a_h = 1, b_h = 1, sample_sig,
                                          is.h.par = TRUE) {
    n_subj <- ncol(der_mc_mat)  ## one t 
    D <- nrow(der_mc_mat) 
    n <- nrow(multi_y)
    mu <- rep(0L, n)

    Kff <- se_ker(H0 = H0, tau = theta[1], h = theta[2])
    nn <- 1
    log_all <- sapply(1:n_subj, function(s) {
        Kdf_lst <- as.list(data.frame(apply(matrix(der_mc_mat[, s], nn, D), 2, 
                                            computeCovDer1, 
                                            idx2 = x_vec, 
                                            tau = theta[1], h = theta[2])))
        Kdf_lst_new <- lapply(Kdf_lst, function(x) {
            matrix(x, nrow = nn, ncol = n)
        })
        
        den_value_lst_new <- lapply(1:D, function(m) {
            Sigma <- Kff - crossprod(Kdf_lst_new[[m]]) / ((theta[1] / theta[2]) ^ 2)
            Psi <- Sigma + diag(sample_sig[m] ^ 2, n)
            mvnfast::dmvn(multi_y[, s], mu = mu, sigma = Psi)
        })
        return(log(sum(unlist(den_value_lst_new)) / D))
    })
    
    log_summ <- sum(log_all)
    
    if(is.h.par) {
        H <- stats::dgamma(theta[2], shape = a_h, rate = b_h, log = TRUE)
        return(-log_summ - H)
    } else {
        return(-log_summ)
    }
}


