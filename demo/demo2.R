# Bayesian P-splines using INLA, 2 dimensional
# smooth + linear, smooth + smooth

# Simulate data
n =  50 # 
set.seed(2013)
z = seq(0, 1, length = n)
x = runif(n)
x = x - mean(x)
y0 = 5*x + sin(14 * z) 
sd = 0.1
y = y0 + rnorm(n) * sd

library(INLA)
library(BayesianPsplines)

# smooth(z)+smooth(x)
mydata2 = mydata = data.frame(y = y, z= z, x = x)
rs2 = smooth_inla(mydata2)
rs2$model # inla object 
rs2$mod$summary.random$idx0 # intercept
coef(rs2)  # intercept, B-spline coefs.
predz = component_plus_residual(rs2)
predx = component_plus_residual(rs2, k=2)

# smooth(z) + linear(x)
mydata = data.frame(y = y, z= z)
data.linear = data.frame(x = x - mean(x))
# make sure fixed (quantitative) predictors have mean 0 
# for component_plus_residual to have a meaningful response scale
hyperGumbel = list(prec = list(prior="pc.prec", param=c(U=1/0.31,alpha=0.01)))
rs = smooth_inla(mydata, data.linear = data.linear)
mod = rs$model # inla object 
mod$summary.random$idx0
coef(rs)  # intercept, linear coefs, B-spline coefs.
pred = component_plus_residual(rs)
# currently no plots available for linear terms...

# smooth(z) only...
mydataz = data.frame(y = y, z= z)
rs3 = smooth_inla(mydataz)
pred3 = component_plus_residual(rs3)

# smooth(z) only with little smoothing
mydataz = data.frame(y = y, z= z)
rs4 = smooth_inla(mydataz,   lambda = 0.001)
pred4 = component_plus_residual(rs4)








