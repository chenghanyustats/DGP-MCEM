deriv_2nd <- function(x, y) {
    deriv_1st(middle_pts(x), deriv_1st(x, y))
}