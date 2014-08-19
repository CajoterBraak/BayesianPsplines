#' @title Utility function lntransform
#'
#' @description
#' \code{lntransform} takes the natural logarithm of the argument vector y which can contains 
#' zeroes or be negative  
#'
#' @details  If y is nonnegative and contains zero(es), it adds the minimum non-zero value, say eps, to y and
#' the result is ln(y+eps) with attribute eps.
#' If y is negative, it subtracts the minimum value and then applies the procedure for nonnegative y. 
#
#' @param x a numeric vector.
#' @param eps value. If set, used for the result ln(y+eps),
#' @return ln(y+eps) with attribute eps.
#' @example /demo/example_lntransform.r

lntransform <- function(x,eps){
  if (!missing(eps)){ 
    eps1 = eps
  } else {
    eps1 = 0
    minx = min(x, na.rm = TRUE)
    if (minx<=0) {
      if (minx==0){
        eps1 = min(x[x > 0], na.rm = TRUE)
      } else{  # <0
        minxx = min(x[x>minx], na.rm = TRUE)
        eps = minxx - minx
        eps1 = -(minx-eps)
      }
    }
  }
  z= log(x+eps1); attr(z, "eps")= eps1;
  z
} 