## compute conditional mean and variance of (y1 | y2) (x_test | x) from multivariate normal
# depends on matrixcalc, Matrix
cond_norm_mean_var <- function(mean_vec = NULL, cov_mat = NULL, 
                               idx1 = NULL, idx2 = NULL,
                               mean_vec_1, mean_vec_2, obs_2,
                               cov_mat_1, cov_mat_2, cov_mat_12) {
    # if (missing(cov_mat_1)) {
    #     cov_mat_1 <- cov_mat[idx1, idx1]
    #     cov_mat_12 <- cov_mat[idx2, idx1, drop = FALSE]
    #     cov_mat_2 <- cov_mat[idx2, idx2]
    #     mean_vec_1 <- mean_vec[idx1]
    #     mean_vec_2 <- mean_vec[idx2]
    # }
    
    if (missing(cov_mat_1)) {
        if (is.null(cov_mat) || is.null(idx1)) {
            stop("must provide cov_mat and idx1")
        } else {
            cov_mat_1 <- cov_mat[idx1, idx1]
        }
    }
    
    if (missing(cov_mat_12)) {
        if (is.null(cov_mat) || is.null(idx1) || is.null(idx2)) {
            stop("must provide cov_mat, idx1 and idx2")
        } else {
            cov_mat_12 <- cov_mat[idx1, idx2, drop = FALSE]
        }
    }
    
    if (missing(cov_mat_2)) {
        if (is.null(cov_mat) || is.null(idx2)) {
            stop("must provide cov_mat and idx2")
        } else {
            cov_mat_2 <- cov_mat[idx2, idx2]
        }
    }
    
    if (missing(mean_vec_1)) {
        if (is.null(mean_vec) || is.null(idx1)) {
            stop("must provide mean_vec and idx1")
        } else {
            mean_vec_1 <- mean_vec[idx1]
        }
    }
    
    if (missing(mean_vec_2)) {
        if (is.null(mean_vec) || is.null(idx2)) {
            stop("must provide mean_vec and idx2")
        } else {
            mean_vec_2 <- mean_vec[idx2]
        }
    }
    
    B <- cov_mat_12 %*% chol2inv(chol(cov_mat_2))
    
    if (mean_vec_1 == 0 && mean_vec_2 == 0) {
        cMu <- B %*% obs_2
    } else {
        cMu <- mean_vec_1 + B %*% (obs_2 - mean_vec_2)
    }
    # cVar <- cov_mat_1 - B %*% t(cov_mat_12)
    cVar <- cov_mat_1 - tcrossprod(B, cov_mat_12)
    # if (!matrixcalc::is.positive.definite(cVar)) {
    #     cVar <- as.matrix(Matrix::nearPD(cVar)$mat)
    # }
    cVar <- as.matrix(Matrix::nearPD(cVar)$mat)
    return(list(condMean = cMu, condVar = cVar))
}


# ## speed test
# mat <- matrix(rnorm(10000), 100, 100)
# mat2 <- matrix(rnorm(1000), 100, 10)
# mat <- t(mat) %*% mat
# 
# system.time(
#     for (i in 1:1000) {
#         # xx <- chol2inv(chol(mat))
#         quad.form.inv(mat, mat2)
#         t(mat2) %*% chol2inv(chol(mat))
#     }
# )
# 
# faster
# system.time(
#     for (i in 1:1000) {
#         xx2 <- t(mat2) %*% chol2inv(chol(mat))
#         xx2 %*% mat2
#     }
# )
# 
# system.time(
#     for (i in 1:1000) {
#         xx2 <- t(mat2) %*% solve(mat)
#         xx2 %*% mat2
#     }
# )

# system.time(
#     for (i in 1:5000) {
#         diag(mat) <- diag(mat) + 1
#     }
# )
# 
# # faster
# system.time(
#     for (i in 1:5000) {
#         mat <- mat + diag(1, 100)
#     }
# )
# 
# system.time(
#     for (i in 1:5000) {
#         mat <- mat + 1 * diag(100)
#     }
# )


# x_mat <- matrix(rnorm(10000), 100, 100)
# y_mat <- matrix(rnorm(10000), 100, 100)
# 
# microbenchmark(tcrossprod(x_mat, y_mat), x_mat %*% t(y_mat), times = 5000)
# microbenchmark(crossprod(x_mat, y_mat), t(x_mat) %*% y_mat, times = 5000)
