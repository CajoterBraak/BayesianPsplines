# Bayesian P-splines using INLA, 1 dimensional smooth

# Simulate data
n =  50
set.seed(2013)
x = seq(0, 1, length = n)
y0 = sin(14 * x)
sd = 0.1
y = y0 + rnorm(n) * sd
library(INLA)
library(BayesianPsplines)
mydata = data.frame(y = y, x= x)
rs = smooth_inla(mydata)
## the default prior used here is:
# hyperGumbel = list(prec = list(prior="pc.prec", param=c(U=1/0.31,alpha=0.01)))
# rs = smooth_inla(mydata, hyperB=hyperGumbel)
# rs = smooth_inla(mydata, lambda = 1.427603) #
##
summary(rs$model) # rs$model is the inla object
pred = component_plus_residual(rs)
names(pred$pred.grid)# predicted values for gridpoints
names(pred$pred.data) # predicted values for data
# as visible in a fit with few predictor points
rs2 = smooth_inla(mydata, ngrid = 10)
pred2 = component_plus_residual(rs2, mplot = 0)
pred2$pred.grid
