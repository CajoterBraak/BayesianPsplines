prepare_basis_P <- function(x, xrange = c(min(x),max(x)), diff.order = 2, ngrid = 100, nseg = 20, degree = 3, eps = 1e-5, grid_with_x = TRUE){
  
  # Prepare basis and penalty matrix for data x and a grid
  #  grid_with_x makes a grid of c(x, seq(min,max,length = ngrid)
  ########
  B = bbase(x, xrange[1], xrange[2], nseg, degree)
  B[abs(B)<eps] = 0
  B = Matrix(B)
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
    indices$x = ngrid+seq_along(x)
  }
  list(x = x, B = B, P = P, diff.order = diff.order, nb = nb, x_grid = x_grid, B_grid = B_grid, indices_B_grid = indices, grid_with_x = grid_with_x, xrange =xrange, nseg = nseg, degree= degree)
}