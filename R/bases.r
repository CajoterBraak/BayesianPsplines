tpower <- function(x, t, p)
# Truncated p-th power function
    (x - t) ^ p * (x > t)


bbase <- function(x, xl = min(x), xr = max(x), nseg = 10, deg = 3, eps = 1e-5){
# Construct B-spline basis
    dx <- (xr - xl) / nseg
    knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
    P <- outer(x, knots, tpower, deg)
    n <- dim(P)[2]
    D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
    B <- (-1) ^ (deg + 1) * P %*% t(D)
    B[abs(B)<eps] = 0
    Matrix(B)
}