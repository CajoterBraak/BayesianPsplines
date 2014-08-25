#' @export
component_plus_residual = function(rs, k=1, transform = function(x){x}, backtransform = function(x){x}, mplot = 2, yscale = 1, ylim = NULL, ...){
  # rs is an object created using a call to smooth_inla
  # k is the number of predictor to fit
  # specialized at present for (beta)binomial data; see transform and backtransform
  # mplot ==1 : plot raw data or proportion (if Ntrials set), 
  # mplot ==2 : plot component _ plus_ residual
  # yscale = multiplier of observed data when mplot = 1 (and backtransform is not used)
  #          use yscale = 100 to obtain percentages
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
  observed.tr = transform(yy) # on transformed scale
  residual = observed.tr - yhat
  if (!is.null(rs$indices_B_grid[[k]]$x) && (ncol(rs$x_grid)>1||!is.null(rs$group)))  {
    irangex = rs$indices_B_grid[[k]]$x 
    componentplus = residual + rs$pred[irangex,"0.5quant"] 
  }else componentplus = observed.tr  # or find closest x_grid to x and use that fit!

  if (mplot){
    if (is.null(ylim)) {pmax = max( backtransform(c(y_hi,componentplus)))
                        pmin = min(backtransform(c(y_lo,componentplus)))
                        plim = c(pmin, pmax)} else plim = ylim
    plot(x_grid, backtransform(y_pred), type = "l", ylim = plim,  ...)
    lines(x_grid, backtransform(y_hi), lty="dashed")
    lines(x_grid, backtransform(y_lo), lty="dashed")
    if (mplot == 2){
      points(rs$x[,k], backtransform(componentplus),... ) # component-plus-residual plot
      } else 
      points(rs$x[,k],yscale*yy,... ) # simple plot if mplot = 1 (useful if y~x)
  } # mplot
  pred.out = data.frame(x_grid,mean=rs$pred[irange,"mean"] , sd=rs$pred[irange,"sd"],y_lo, y_pred, y_hi) # on transformed scale
  names(pred.out)[-1]= names(rs$pred)[-c(1,7,8)] 
  data.out = data.frame(x=rs$x[,k],y = rs$y, Ntrials = rs$Ntrials, group = rs$group, observed.tr, yhat,residual,componentplus,rs$fitted) # on transformed scale
  list(pred.grid = pred.out, pred.data = data.out)
}