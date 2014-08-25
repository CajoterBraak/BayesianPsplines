#' @title extend_range
#'
#' @description
#' \code{extend_range} extends the range of the predictors. Extension and bounds can be specified
#'   
#'
#' @param x  A numeric vector 
#' @param extend A positive scalar (default = 0.1, meaning 10\% extension at both sides)
#' @param bounds If set, extension is truncated. A non-negative restriction is set by \code{bounds = c(0,NA)}. 
#' 
#' @details \code{extend_range} multiplies the range in the data by
#' the \code{extend} paramater  and subtracts/adds this to the upper and lower limit.
#' If \code{bounds} parameter is given, the intersection with the bounds is taken.

extend_range <- function(x, extend = 0.1, bounds){
  rr =range(x)
  dx = (rr[2]-rr[1])*extend
  xrange = c(rr[1]-dx, rr[2] +dx)
  if (!missing(bounds)){
    if(!is.na(bounds[1])) if (xrange[1]<bounds[1]) xrange[1]= bounds[1]
    if(!is.na(bounds[2])) if (xrange[2]>bounds[2]) xrange[2]= bounds[2]
  }
  xrange
}