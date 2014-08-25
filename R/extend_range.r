extend_range <- function(x, extend = 0.2, bounds){
  rr =range(x)
  dx = (rr[2]-rr[1])*extend
  xrange = c(rr[1]-dx, rr[2] +dx)
  if (!missing(bounds)){
    if(!is.na(bounds[1])) if (xrange[1]<bounds[1]) xrange[1]= bounds[1]
    if(!is.na(bounds[2])) if (xrange[2]>bounds[2]) xrange[2]= bounds[2]
  }
  xrange
}