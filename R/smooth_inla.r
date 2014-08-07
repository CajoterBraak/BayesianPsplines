#' @title Bayesian P-splines and GAMs using INLA
#'
#' @description
#' \code{smooth_inla} Smooths data y in a data frame with respect to the
#' zero, one or more predictors with or without an additional random group factor
#'   
#'
#' @details   aha
#' @param data A data frame with columns y [group] x1 x2 ... , all numeric except
#' group which is a factor and is optional. The response y should always be the first 
#' column, and, if present, group second. All further vectors are used a predictors
#' in a GAM. 
#' @param Ntrials Number of binomial trials with same length as y, if family is binomial
#' @param family The family argument passed to INLA, see \code{\link[INLA]{inla}} 
#' @param hyperB A list of hyperparameter for each of the predictors
#' @param xrange A matrix with in rows the min and max of the range of 
#' the B-spline basis for each predictor. Default: the data range of each predictor
#' @param ngrid Number of grid points  for prediction. Default 100
 
#' @examples Example program in BayesPspline.r 
#' @export


smooth_inla <- function(data,  Ntrials =NULL, family = "gaussian", hyperB, 
     xrange= NULL, ngrid = 100, nseg = 20 , degree = 3, 
    diff.order = rep(2,ncol(data)), grid_with_x = TRUE, q = c(0.05, 0.5, 0.95), verbose = FALSE ){
# smooth_inla assumes that data consists of the columns (in this order)
#  y and optionally a group factor
#  zero, one or more predictor variables to be used
#   and no other variables!!
# author: Cajo.terBraak@wur.nl 
# copyright: CC-BY (Creative Commons Attribution 3.0 Netherlands)
  basisP = list()
  if (is.null(data$group)) dataX = data[,-1, drop = FALSE] else dataX = data[,-c(1,2), drop = FALSE]
  if (is.null(data$group)) Group = NULL else Group = model.matrix(~-1+group, data)
  Amlist = list(intercept = matrix(1, nrow = length(data$y), ncol = 1))
  for (k in seq_len(ncol(dataX))) {
    if (is.null(xrange)) xrange.k = c(min(dataX[[k]]),max(dataX[[k]])) else 
      if (is.matrix(xrange)) xrange.k = xrange[k,] else xrange.k = xrange
    basisP[[k]] = prepare_basis_P(dataX[[k]], xrange.k,diff.order[k], ngrid, nseg, degree, grid_with_x=grid_with_x)     
    Amlist[[k+1]] = basisP[[k]]$B
  }
  Amlist$Group = Group
  pA = sapply(Amlist,ncol)
  names(pA)= names(Amlist)
  npar = sum(pA); n.pA = length(pA)
  ipar = c(0,cumsum(pA[-n.pA]))
  datalist = list()  #matrix(NA,nrow = npar,ncol = n.pA)
  for (i in seq_len(n.pA)){
    idxi = rep(NA, npar); idxi[ipar[i]+seq_len(pA[i])]=seq_len(pA[i]); datalist[[i]] =idxi
  }
  names(datalist)= paste("idx",seq_len(n.pA)-1, sep="") # idx0, idx1, idx2,....
  A.matrix = Amlist[[1]]
  for (i in seq_len(length(Amlist))[-1]) A.matrix= Matrix::cbind2(A.matrix,Amlist[[i]])
  nB = ncol(A.matrix) # before the Aplus trick
  A.matrix = rBind(Diagonal(nB), A.matrix) 
  # this is Aplus: to allow predictions with credible regions via B_grid to work
  yplus = c(rep(NA, nB), data$y)
  if (!is.null(Ntrials)) Ntrials_plus = c(rep(1, nB), Ntrials)else Ntrials_plus = NULL
  datalist$y =yplus


B_grid = matrix(1, nrow =1, ncol = 1)
indices = list()
if (grid_with_x) ngridn = ngrid + nrow(dataX) else ngridn = ngrid
for (k in seq_len(ncol(dataX))) {
     nb =  basisP[[k]]$nb          
     B_grid = cBind(B_grid,matrix(0,nrow = nrow(B_grid), ncol = nb ))
     B_grid = rBind(B_grid, cBind(matrix(0,nrow = ngridn, ncol = ncol(B_grid)-nb ), basisP[[k]]$B_grid) )
     indices[[k]] = sapply(basisP[[k]]$indices_B_grid, function(x){x+ngridn*(k-1) + 1 })
}
B_grid[,1]= 1 # intercept added to all predictions
#
# B_grid is a block diagonal matrix with first column = 1
# not counting the 1.1 entry, the k-th block contains the basis for the k-th predictor and its grid
# the indices are indices
# indices_B_grid[[k]]$xgrid and $x give the rows in B_grid of the grid and the data of the k-th predictor
# 
lc.grid = inla.make.lincombs(Predictor = B_grid)

hyper.group = list(prec = list(prior = "loggamma", param = c(1, 0.01)))
hyper.fixed = list(prec = list(initial = log(0.001), fixed=TRUE))
prior.observation.precision = list(prior = "flat", param = numeric(0))

# make formula
formula.txt = paste(" formula.P = y ~ -1 + f(idx0,  model=\"iid\", hyper = hyper.fixed, constr = FALSE)") # intercept
for (k in seq_len(ncol(dataX))) {
  txt = paste("+ f(idx", k,", model=\"generic\", Cmatrix = basisP[[",k,"]]$P, constr = TRUE,  hyper = hyperB)", sep = "")
  formula.txt= paste(formula.txt, txt,sep = "")
}
if (!is.null(data$group)) {
  txt = paste("+ f(idx", ncol(dataX)+1,", model=\"iid\", constr = TRUE,  hyper = hyper.group)", sep = "")
  formula.txt= paste(formula.txt, txt,sep = "")
}

formula.P = as.formula(formula.txt)

contr.comp = inla.set.control.compute.default()
contr.comp$dic = TRUE
mod.P = NULL
if (family %in% c("binomial","betabinomial")){
mod.P = inla(formula.P, family = family, Ntrials = Ntrials_plus, data = datalist, 
      control.predictor = list(A = A.matrix, compute = TRUE, quantiles = q, link = NULL), 
              control.compute = contr.comp,
              lincomb = lc.grid)
} else if (family %in% c("gaussian","t","lognormal")){
  mod.P = inla(formula.P, family = family, data = datalist, 
               control.predictor = list(A = A.matrix, compute = TRUE, quantiles = q, link = NULL), 
               control.family = list(hyper = list(prec = prior.observation.precision)),
               control.compute = contr.comp,
               lincomb = lc.grid)
} else { 
mod.P = inla(formula.P, family = family, data = datalist, 
             control.predictor = list(A = A.matrix, compute = TRUE, quantiles = q, link = NULL), 
             control.compute = contr.comp,
             lincomb = lc.grid) 
}

if (verbose) summary(mod.P)
# get the default quantiles
intercept = mod.P$summary.lincomb.derived[1,]
Pred.rw = mod.P$summary.lincomb.derived

x_grid = NULL
for (k in seq_along(basisP))x_grid = cbind(x_grid, basisP[[k]]$x_grid)
fitted = mod.P$summary.fitted.values[nB + (seq_len(nrow(dataX))),]

list(model_inla = mod.P, x = dataX, y=data$y, Ntrials=Ntrials, intercept = intercept, fitted = fitted, pred = Pred.rw, x_grid = x_grid, B_grid=B_grid, indices_B_grid = indices)
}




