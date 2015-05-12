#' @title Bayesian P-splines and GAMs using INLA
#'
#' @description
#' \code{smooth_inla} smooths data y in a data frame with respect to the
#' zero, one or more predictors with or without an additional random group factor.
#'   
#'
#' @details   \code{smooth_inla} is a wrapper that uses INLA for the calculations. see \code{\link[INLA]{inla}}
#' The predictors are transformed into B-splines and fitted using a penalty matrix tau*P 
#' with P =D\'D with D differences of \code{diff.order}.
#' @param data A data frame with columns y [group] x1 x2 ... , all numeric except
#' group which is a factor and is optional. The response y should always be the first 
#' column, and, if present, group second. All further vectors are used a predictors
#' in a GAM. Response and predictors can have user-defined names, but the grouping factor must have name 'group'
#' @param Ntrials Number of binomial trials with same length as y, if family is binomial
#' @param family The family argument passed to INLA, see \code{\link[INLA]{inla}} 
#' @param hyperB A list of hyperparameter for each of the predictors. If missing,
#' the PC-prior for precision, i.e. a type-2 Gumbel with (default) lambda = 1.427 (sd(u) = 1, u=3.22, alpha = 0.01). 
#' See inla.doc("pc.prec"). 
#' @param lambda parameter of the type-2 Gumbel, used if hyperB is missing (default 1.428)
#' @param offset A optional offset for the linear predictor 
#' @param offset_par Value of the offset used for predictions.
#' @param xrange A matrix with in rows the min and max of the range of 
#' the B-spline basis for each predictor. Default: \code{\link{extend_range}} giving 10\% extension at both sides
#' @param ngrid Number of grid points  for prediction. Default 100
#' @examples Example program in BayesianPsplines\demo\BayesPspline.r 
#' @export

smooth_inla <- function(data,  Ntrials, family = "gaussian", hyperB, lambda, weights, offset,
     offset_par, xrange, ngrid = 100, nseg = 20 , degree = rep(3,ncol(data)), 
    diff.order = rep(2,ncol(data)), grid_with_x = TRUE, quantiles = c(0.05, 0.5, 0.95), verbose = FALSE ){
# smooth_inla assumes that data consists of the columns (in this order)
#  y and optionally a group factor
#  zero, one or more predictor variables to be used
#   and no other variables!!
# author: Cajo.terBraak@wur.nl 
# copyright: CC-BY (Creative Commons Attribution 3.0 Netherlands)
  basisP = list()
  if (is.null(data$group)) dataX = data[,-1, drop = FALSE] else dataX = data[,-c(1,2), drop = FALSE]
  if (is.null(data$group)) Group = NULL else Group = model.matrix(~-1+group, data)
  Amlist = list(intercept = matrix(1, nrow = nrow(data), ncol = 1))
  for (k in seq_len(ncol(dataX))) {
    if (missing(xrange)) xrange.k = extend_range(dataX[[k]]) else 
      if (is.matrix(xrange)) xrange.k = xrange[k,] else xrange.k = xrange
    basisP[[k]] = prepare_basis_P(dataX[[k]], xrange.k,diff.order[k], ngrid, nseg, degree[k], grid_with_x=grid_with_x)     
    Amlist[[k+1]] = basisP[[k]]$B
  }
  Amlist$Group = Group
  if (missing(hyperB)){
    alpha = 0.01
    if (missing(lambda)){   # lambda = -log(alpha)/U =1.428
      sd_marg_u = 1; U = sd_marg_u /0.31
    } else  U = -log(alpha)/lambda
     hyperB = list(prec = list(prior="pc.prec", param=c(U,0.01)))
  }
  if (length(hyperB)==1) hyperB = list(hyperB,hyperB,hyperB,hyperB,hyperB,hyperB,hyperB,hyperB,hyperB,hyperB)
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
  nB = ncol(A.matrix)
  if (missing(offset)) offset = NULL
  if (missing(Ntrials)) Ntrials = NULL
  datalist$y =data[,1] # response y


B_grid = matrix(1, nrow =1, ncol = 1)
indices = list()
if (grid_with_x) ngridn = ngrid + nrow(dataX) else ngridn = ngrid
for (k in seq_len(ncol(dataX))) {
     nb =  basisP[[k]]$nb          
     B_grid = cBind(B_grid,matrix(0,nrow = nrow(B_grid), ncol = nb ))
     B_grid = rBind(B_grid, cBind(matrix(0,nrow = ngridn, ncol = ncol(B_grid)-nb ), basisP[[k]]$B_grid) )
     indices[[k]] = sapply(basisP[[k]]$indices_B_grid, function(x){x+ngridn*(k-1) + 1 }, simplify = FALSE)
}
B_grid[,1]= 1 # intercept added to all predictions
#
# B_grid is a block diagonal matrix with first column = 1
# not counting the 1.1 entry, the k-th block contains the basis for the k-th predictor and its grid
# the indices are indices
# indices_B_grid[[k]]$xgrid and $x give the rows in B_grid of the grid and the data of the k-th predictor
#
# In the INLA A-matrix approach, the output is c(eta*, eta), so shift B_grid...
B_grid_plus =cBind(Matrix(0,nrow = nrow(B_grid), ncol = nrow(A.matrix)),B_grid)

lc.grid = inla.make.lincombs(Predictor = B_grid_plus)

hyper.group = list(prec = list(prior = "loggamma", param = c(1, 0.01)))
hyper.fixed = list(prec = list(initial = log(0.001), fixed=TRUE))
prior.observation.precision = list(prior = "flat", param = numeric(0))

# make formula
formula.txt = paste(" formula.P = y ~ -1 + f(idx0,  model=\"iid\", hyper = hyper.fixed, constr = FALSE)") # intercept
for (k in seq_len(ncol(dataX))) {
  txt = paste("+ f(idx", k,", model=\"generic\", Cmatrix = basisP[[",k,"]]$P, constr = TRUE,  hyper = hyperB[[",k,"]], diagonal = 1e-4)", sep = "")
  formula.txt= paste(formula.txt, txt,sep = "")
}
if (!is.null(data$group)) {
  txt = paste("+ f(idx", ncol(dataX)+1,", model=\"iid\", constr = TRUE,  hyper = hyper.group)", sep = "")
  formula.txt= paste(formula.txt, txt,sep = "")
}

formula.P = as.formula(formula.txt)
print(formula.P)

contr.comp = inla.set.control.compute.default()
contr.comp$dic = TRUE
mod.P = NULL
start = proc.time()[3]
if (missing(weights)) {
  if (family %in% c("binomial","betabinomial")){
  mod.P = inla(formula.P, family = family, Ntrials = Ntrials, data = datalist, offset = offset,
        control.predictor = list(A = A.matrix, compute = TRUE,link = NULL), quantiles = quantiles, 
                control.compute = contr.comp,
                lincomb = lc.grid, verbose = verbose)
  } else if (family %in% c("gaussian","t","lognormal")){
    mod.P = inla(formula.P, family = family, data = datalist, offset = offset,
                 control.predictor = list(A = A.matrix, compute = TRUE,link = NULL), quantiles = quantiles,
                 control.family = list(hyper = list(prec = prior.observation.precision)),
                 control.compute = contr.comp,
                 lincomb = lc.grid, verbose = verbose)
  } else { 
  mod.P = inla(formula.P, family = family, data = datalist, offset = offset,
               control.predictor = list(A = A.matrix, compute = TRUE,link = NULL), quantiles = quantiles, 
               control.compute = contr.comp,
               lincomb = lc.grid, verbose = verbose) 
  }
} else { 
  inla.setOption(enable.inla.argument.weights=TRUE)
  if (family %in% c("binomial","betabinomial")){
    mod.P = inla(formula.P, family = family, Ntrials = Ntrials, data = datalist, offset = offset,
                 control.predictor = list(A = A.matrix, compute = TRUE,link = NULL), quantiles = quantiles, 
                 control.compute = contr.comp,
                 lincomb = lc.grid, verbose = verbose)
  } else if (family %in% c("gaussian","t","lognormal")){
    mod.P = inla(formula.P, family = family, data = datalist, weights= weights, offset = offset,
                 control.predictor = list(A = A.matrix, compute = TRUE,link = NULL), quantiles = quantiles, 
                 control.family = list(hyper = list(prec = prior.observation.precision)),
                 control.compute = contr.comp,
                 lincomb = lc.grid, verbose = verbose)
  } else { 
    mod.P = inla(formula.P, family = family, data = datalist, weights= weights, offset = offset,
                 control.predictor = list(A = A.matrix, compute = TRUE,link = NULL), quantiles = quantiles, 
                 control.compute = contr.comp,
                 lincomb = lc.grid, verbose = verbose) 
  }       
}
end = proc.time()[3]
time.P = end - start
if (verbose) summary(mod.P)
# get the default quantiles
intercept = mod.P$summary.lincomb.derived[1,]
Pred.rw = mod.P$summary.lincomb.derived

x_grid = NULL
for (k in seq_along(basisP))x_grid = cbind(x_grid, basisP[[k]]$x_grid)
#fitted = mod.P$summary.fitted.values[nB + (seq_len(nrow(dataX))),]
fitted = mod.P$summary.linear.predictor[seq_len(nrow(dataX)),]
rownames(fitted)= rownames(data)
coef = mod.P$summary.linear.predictor[nrow(dataX) + 1:nB,] # inclusive intercept
rownames(coef)= c("(Intercept)", 1:(nB-1))

list(model_inla = mod.P, x = dataX, y=datalist$y, Ntrials=Ntrials, offset = offset, group= data$group, intercept = intercept, coefficients = coef,fitted = fitted, pred = Pred.rw, x_grid = x_grid, B_grid=B_grid, indices_B_grid = indices, time = time.P)
}




