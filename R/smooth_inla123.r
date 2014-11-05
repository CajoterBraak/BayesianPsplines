# tutororial R functions for Bayesian P-splines using INLA
# for smoooth_inla.r in real work
# author: Cajo.terBraak@wur.nl 
# copyright: CC-BY (Creative Commons Attribution 3.0 Netherlands)
smooth_inla0 <- function(x,y, Ntrials, offset, family, hyperB = list(prec = list(prior = "loggamma", param = c(1, 0.01))),    
          xrange, diff.order = 2, ngrid = 100, nseg = 20, degree = 3, q = c(0.05, 0.5, 0.95) ){
  # tutorial code showing how to fit a single spline, not using an intercept
  # Prepare basis and penalty matrix
  if (missing(xrange))xrange = extend_range(x)
  B = bbase(x, xrange[1], xrange[2], nseg, degree)
  nb = ncol(B)
  D = diff(diag(nb), diff = diff.order)
  P = t(D) %*% D 
  x_grid = seq(xrange[1], xrange[2], length = ngrid)
  B_grid = bbase(x_grid, xrange[1], xrange[2], nseg, degree)
  # In the INLA A-matrix approach, the output is c(eta*, eta), so shift B_grid...
  B_grid_plus =cBind(Matrix(0,nrow = nrow(B_grid), ncol = nrow(B)),B_grid)
  
  formula.P = y ~ -1 + f(id.b, model="generic", Cmatrix = P, constr = FALSE, hyper= hyperB)
  datalist = list(y = y, id.b = 1:nb)
  lc.grid = inla.make.lincombs(Predictor = B_grid_plus)
  mod.P = inla(formula.P,  Ntrials = Ntrials, offset = offset,family = family, data = datalist, control.predictor = list(A = B, compute = TRUE), quantiles = q, lincomb = lc.grid)
  Pred = mod.P$summary.lincomb.derived
  list(model_inla = mod.P, pred = Pred, x_grid = x_grid, B_grid=B_grid)
}

smooth_inla1 <- function(x,y, Ntrials, offset, family, hyperB = list(prec = list(prior = "loggamma", param = c(1, 0.01))), 
      xrange, diff.order = 2, ngrid = 100, nseg = 20, degree = 3, q = c(0.05, 0.5, 0.95) ){

  # tutorial code showing how to fit a single spline, using an intercept 
  # Prepare basis and penalty matrix
if (missing(xrange))xrange = extend_range(x)
B = bbase(x, xrange[1], xrange[2], nseg, degree)
nb = ncol(B)
D = diff(diag(nb), diff = diff.order)
P = t(D) %*% D
# add intercept as first column
A = Matrix::cbind2(1,B)
# add ncol(A) at the top of A and to y (and Ntrials)
# used for creating predictions
datalist = list(y= y, idx0 = c(1,rep(NA,nb)), idx1 = c(NA, 1:nb))

# set up an equispaced grid for making prediction 
x_grid = seq(xrange[1], xrange[2], length = ngrid)
B_grid = bbase(x_grid, xrange[1], xrange[2], nseg, degree)
# add intercept as first column
B_grid = Matrix::cbind2(1,B_grid)
# add a row (1,0,0,0,....) for the intercept
B_grid = Matrix(rBind(c(1, rep(0,nb)),B_grid))
# In the INLA A-matrix approach, the output is c(eta*, eta), so shift B_grid...
B_grid_plus =cBind(Matrix(0,nrow = nrow(B_grid), ncol = nrow(A)),B_grid)

# set up INLA call using the A-matrix approach
hyper.fixed = list(prec = list(initial = log(0.001), fixed=TRUE))
formula.P = y ~ -1 +
  f(idx0,  model="iid", hyper = hyper.fixed, constr = FALSE) + # intercept
  f(idx1, model="generic", Cmatrix = P, constr = TRUE,  hyper = hyperB)
lc.grid = inla.make.lincombs(Predictor = B_grid_plus)
mod.P = inla(formula.P,  Ntrials = Ntrials, offset = offset, family = family, data = datalist, control.predictor = list(A = A, compute = TRUE), quantiles = q, lincomb = lc.grid)
Pred = mod.P$summary.lincomb.derived
list(model_inla = mod.P, pred = Pred, x_grid = x_grid, B_grid=B_grid)
}

smooth_inla2 <- function(x,y, Ntrials, offset, family, hyperB = list(prec = list(prior = "loggamma", param = c(1, 0.01))), 
                         xrange, diff.order = 2, ngrid = 100, nseg = 20, degree = 3, q = c(0.05, 0.5, 0.95) ){
  
  # tutorial code showing how to fit a single spline, using an intercept and linear term for x
  # Prepare basis and penalty matrix
  if (missing(xrange)) xrange = extend_range(x)
  B = bbase(x, xrange[1], xrange[2], nseg, degree)
  nb = ncol(B)
  D = diff(diag(nb), diff = diff.order)
  P = t(D) %*% D
  # add intercept and x as first column
  A = Matrix::cbind2(cbind(1,x),B)
  # add ncol(A) at the top of A and to y (and Ntrials)
  # used for creating predictions
  datalist = list(y= y, idx0 = c(1:2,rep(NA,nb)), idx1 = c(NA, NA, 1:nb))
  
  # set up an equispaced grid for making prediction 
  x_grid = seq(xrange[1], xrange[2], length = ngrid)
  B_grid = bbase(x_grid, xrange[1], xrange[2], nseg, degree)
  # add intercept and x_grid as first two columns
  B_grid = Matrix::cbind2(cbind(1,x_grid),B_grid)
  # add a row (1,0,0,0,....) for the intercept
  # add a row (0,1,0,0,....) for linear effect of x
  intx = rbind(c(1, rep(0,nb+1)),c(0,1, rep(0,nb)))
  B_grid = Matrix(rBind(intx,B_grid))
  # In the INLA A-matrix approach, the output is c(eta*, eta), so shift B_grid...
  B_grid_plus =cBind(Matrix(0,nrow = nrow(B_grid), ncol = nrow(A)),B_grid)
  
  # set up INLA call using the A-matrix approach
  hyper.fixed = list(prec = list(initial = log(0.001), fixed=TRUE))
  formula.P = y ~ -1 +
    f(idx0,  model="iid", hyper = hyper.fixed, constr = FALSE) + # intercept
    f(idx1, model="generic", Cmatrix = P, constr = TRUE,  hyper = hyperB)
  lc.grid = inla.make.lincombs(Predictor = B_grid_plus)
  mod.P = inla(formula.P,  Ntrials = Ntrials,offset = offset, family = family, data = datalist, control.predictor = list(A = A, compute = TRUE), quantiles = q, lincomb = lc.grid)
  Pred = mod.P$summary.lincomb.derived
  list(model_inla = mod.P, pred = Pred, x_grid = x_grid, B_grid=B_grid)
}

smooth_inla3 <- function(x,group, y, Ntrials, offset, family, hyperB = list(prec = list(prior = "loggamma", param = c(1, 0.01))), 
                         xrange, diff.order = 2, ngrid = 100, nseg = 20, degree = 3, q = c(0.05, 0.5, 0.95), grid_with_x = TRUE){
  # tutorial code showing how to fit a single spline, using an intercept and 
  # an extra grouping factor
  # Prepare basis and penalty matrix
  if (missing(xrange))xrange = extend_range(x)
  basisP = prepare_basis_P0(x, xrange, ngrid, diff.order, nseg, degree, grid_with_x=grid_with_x)
  nb = basisP$nb
  # add intercept as first column
  A = Matrix::cbind2(1,basisP$B)
  # add group indicator matrix at the end
  ng = nlevels(group)
  A = Matrix::cbind2(A,model.matrix(~-1+group ))
  # add ncol(A) at the top of A and to y (and Ntrials)
  # used for creating predictions
  datalist = list(y= y, idx0 = c(1,rep(NA,nb+ng)), idx1 = c(NA, 1:nb, rep(NA,ng)),
                  idg = c(rep(NA,(nb+1)), rep(1:ng)))
  
  # use the equispaced grid with x for making prediction 
  B_grid =  basisP$B_grid
  # add intercept as first column
  B_grid = Matrix::cbind2(1,B_grid)
  # add a row (1,0,0,0,....) for the intercept
  B_grid = Matrix(rBind(c(1, rep(0,nb)),B_grid))
  B_grid_plus =cBind(Matrix(0,nrow = nrow(B_grid), ncol = nrow(A)),B_grid)
  # set up INLA call using the A-matrix approach
  hyper.fixed = list(prec = list(initial = log(0.001), fixed=TRUE))
  hyper.group = list(prec = list(prior = "loggamma", param = c(1, 0.01)))
  formula.P = y ~ -1 +
    f(idx0,  model="iid", hyper = hyper.fixed, constr = FALSE) + # intercept
    f(idx1, model="generic", Cmatrix = basisP$P, constr = TRUE,  hyper = hyperB) +
    f(idg, model = "iid", constr = TRUE, hyper = hyper.group) 
  lc.grid = inla.make.lincombs(Predictor = B_grid_plus)
  mod.P = inla(formula.P,  Ntrials = Ntrials,offset = offset, family = family, data = datalist, control.predictor = list(A = A, compute = TRUE), quantiles = q, lincomb = lc.grid)
  Pred = mod.P$summary.lincomb.derived
  fitted = mod.P$summary.fitted.values[seq_along(y),]
  list(model_inla = mod.P, fitted = fitted, pred = Pred, x_grid = basisP$x_grid, B_grid=B_grid, indices = basisP$indices)
}

prepare_basis_P0 <- function(x, xrange= c(0,1), ngrid = 100, diff.order = 2, nseg = 20, degree = 3,  eps = 1e-5,  grid_with_x = TRUE){
  # Prepare basis and penalty matrix for data x and a grid
  #  Bgrid_with_x makes a grid of c(x, seq(min,max,length = ngrid)
  B = bbase(x, xrange[1], xrange[2], nseg, degree)
  nb = ncol(B)
  D = diff(diag(nb), diff = diff.order)
  P = t(D) %*% D
  x_grid = seq(xrange[1], xrange[2], length = ngrid)
  B_grid = bbase(x_grid, xrange[1], xrange[2], nseg, degree)
  B_grid[abs(B_grid)<eps] = 0
  B_grid = Matrix(B_grid)
  indices = list(grid = seq_len(ngrid))
  if (grid_with_x) {
    B_grid = rBind(B_grid,B)
    indices$x  = ngrid+seq_along(x)
  }
  list(x = x, B = B, P = P, nb = nb, x_grid = x_grid, B_grid = B_grid, indices = indices)
}
