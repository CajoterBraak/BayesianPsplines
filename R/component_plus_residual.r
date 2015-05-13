#' @title Component plus residual plots for smooth_inla fits
#'
#' @description
#' \code{component_plus_residual} calculates and plots the smooth produced by \code{\link{smooth_inla}}. 
#' 
#'
#' @details   \code{component_plus_residual} postprocesses the output of \code{\link{smooth_inla}} that also
#' contains the inla object.  See also the example in \code{\link{smooth_inla}}.
#'@param  rs is an object created by smooth_inla
#'@param  k is the number of predictor to fit (default 1)
#'@param  transform a function to transform the observed scale to the linear predictor (default identity) 
#'@param  backtransform a function to transform from the linear predictor to the observed scale (default identity) 
#'@param  mplot what to plot.  2  (default): plot component plus residual, 1: plot raw data or proportion (if Ntrials set)
#'0: no plotting
#'@param yscale multiplier for plotting observed data when mplot = 1; when Ntrials is set, yscale = 100 with mplot= 1 then plots row percentages 
#'@examples 
#'# See example in smooth_inla for mydata2
#'
#'  rs2 = smooth_inla(mydata2)
#'  predz = component_plus_residual(rs2)
#'  predx = component_plus_residual(rs2, k=2)
#'# See demo/BayesPspline.r for use of transform/backtransform
#' @export
component_plus_residual = function(rs, k=1, transform = function(x){x}, backtransform = function(x){x}, mplot = 2, yscale = 1, ylim = NULL, ...){
  x_grid = rs$x_grid[,k] ; irange = rs$indices[[k]]$grid#
  y_pred_mean = rs$pred[irange,"mean"] 
  y_pred_sd =  rs$pred[irange,"sd"]
  y_pred = rs$pred[irange,5] # median
  y_lo =  rs$pred[irange,4] # lower quantile
  y_hi =  rs$pred[irange,6] # upper quantile
  y_pred_mean = rs$pred[irange,"mean"] # mean on transformed scale
  y_pred_sd = rs$pred[irange,"sd"] # sd on transformed scale
  # points(x, 100*y/Ntrials) # raw data
  # form residuals on the transformed scale
  yhat = rs$fitted[,"0.5quant"] # transformed scale
  if (is.null(rs$Ntrials)) yy = rs$y  else yy = rs$y/rs$Ntrials
  observed.tr = transform(yy)  # on transformed scale
  residual = observed.tr - yhat 
  if (!is.null(rs$indices_B_grid[[k]]$x) && (ncol(rs$x_grid)>1||!is.null(rs$group)|| nrow(rs$model$summary.random$idx0)>1))  {
    irangex = rs$indices_B_grid[[k]]$x 
    componentplus = residual + rs$pred[irangex,"0.5quant"] 
  }else componentplus = observed.tr  # or find closest x_grid to x and use that fit!
  if (!is.null(rs$offset)) { 
    componentplus = componentplus - rs$offset 
   # if (rs$model$".args"$family =="gaussian") yy = yy - rs$offset else yy = yy/exp(rs$offset)
  }
  
  if (mplot){
    # yy for simple plot
    if (mplot == 2) yyy = backtransform(componentplus) else yyy = yscale*yy
    if (is.null(ylim)) {pmax = max( c(backtransform(y_hi),yyy), na.rm = TRUE)
                        pmin = min( c(backtransform(y_lo),yyy), na.rm = TRUE)
                        plim = c(pmin, pmax)} else plim = ylim
    plot(x_grid, backtransform(y_pred), type = "l", ylim = plim,  ...)
    lines(x_grid, backtransform(y_hi), lty="dashed")
    lines(x_grid, backtransform(y_lo), lty="dashed")
    points(rs$x[,k],yyy,... ) # component-plus-residual plot
  } # mplot
  pred.out = data.frame(x_grid,mean=rs$pred[irange,"mean"] , sd=rs$pred[irange,"sd"],y_lo, y_pred, y_hi) # on transformed scale
  names(pred.out)[-1]= names(rs$pred)[-c(1,7,8)] 
  data.out = data.frame(x=rs$x[,k],y = rs$y, observed.tr, yhat,residual,componentplus,rs$fitted) # on transformed scale
  data.out$Ntrials = rs$Ntrials
  data.out$group = rs$group
  list(pred.grid = pred.out, pred.data = data.out)
}