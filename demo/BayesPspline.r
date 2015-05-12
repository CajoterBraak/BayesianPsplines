# R script for Bayesian P-splines using INLA
# author: Cajo.terBraak@wur.nl 
# copyright: CC-BY (Creative Commons Attribution 3.0 Netherlands
rm(list=ls(all=TRUE))

library(INLA)
library(faraway) # for ilogit 
library(BayesianPsplines)

Data <- read.csv("../data/germination_Ranlin_transect_wl_ntot.csv")[,-1]
summary(Data)
names(Data)

fam = "betabinomial"

sd_marg_u = 1; U = sd_marg_u /0.31; alpha =  0.01
lambda = -log(alpha)/U;
hyperGumbel = list(prec = list(prior="pc.prec", param=c(U,alpha)))

mydata = data.frame(y = Data$y, group = Data$transect, x = Data$wl)
Ntrials = Data$Ntrials
rs = NULL  
rs = smooth_inla(mydata, Ntrials, family = fam, hyperB= hyperGumbel, diff.order = 3)
# intercept
intercept = rs$pred[1,]
dic = rs$model$dic$dic # dic

# transform from proportion to logit scale
transform <- function(x){x = ifelse(abs(x)<0.025,0.025,x);x = ifelse(abs(x-1)<0.025,0.975,x); logit(x)}

# backtransform from logit scale to a percentage scale
backtransform <- function(x){100*ilogit(x)}

species= "Ran lin"; my.ylab = "germinated (%)"; my.xlab = "water level (m)" 
ylim =c(-5,105)
k = 1 # first predictor

pred = component_plus_residual(rs, k, transform, backtransform, ylim= ylim, ylab= my.ylab,xlab = my.xlab , main = species)
## plot data points as percentages by
#pred = component_plus_residual(rs, k, transform, backtransform,  mplot = 1, yscale=100, ylim= ylim, ylab= my.ylab,xlab = my.xlab , main = species)

mydata = data.frame(y = Data$y, group = Data$transect, x1 = Data$wl, x2 = Data$log_ntot)
rs = NULL  
rs = smooth_inla(mydata, Ntrials, family = fam, hyperB= hyperGumbel, diff.order = c(3,3))
# intercept
intercept = rs$pred[1,]
my.xlab = c("water level (m)" ,"ln ntot")
par(mfrow=c(1,2))
k = 1 # first predictor
pred1 = component_plus_residual(rs, k,transform, backtransform, ylim= ylim, ylab= my.ylab,xlab = my.xlab[1], main = species)
k = 2 # second predictor
pred2 = component_plus_residual(rs, k,transform, backtransform, ylim= ylim, ylab= my.ylab,xlab = my.xlab[2], main = species)
rs$model$summary.hyper

#
# smoothing using simple functions that show the ideas of how the wrapper was build.
# See 
# http://onlinelibrary.wiley.com/store/10.1111/1365-2435.12441/asset/supinfo/fec12441-sup-0004-AppendixS3.pdf?v=1&s=ad9af49e363328c26f3353441dc02a63be984ad0
#
x= Data$wl; y = Data$y; Ntrials = Data$Ntrials; group = Data$transect

hyperB = list(prec = list(prior = "loggamma", param = c(1, 0.001)))

rs0 = NULL
#rs0 = smooth_inla0(x,y, Ntrials, family = fam, hyperB= hyperB, xrange= xrange, diff.order = 3)
rs0 = BayesianPsplines:::smooth_inla0(x,y, Ntrials, family = fam, hyperB= hyperB, diff.order = 3)

rs1 = NULL  
rs1 = BayesianPsplines:::smooth_inla1(x,y, Ntrials, family = fam, hyperB= hyperB, diff.order = 3)
rs2 = NULL  
rs2 = BayesianPsplines:::smooth_inla2(x,y, Ntrials, family = fam, hyperB= hyperB, diff.order = 3)


pred0 = rs0$pred
pred1 = rs1$pred[-1,]
pred2 = rs2$pred[-c(1,2),]
# nearly identical fits, as intended; small differences due to???

pmax = max(pred0[, 6])
pmin = min(pred0[, 4])
plim = c(pmin, pmax)
plot(rs0$x_grid, pred0[,5], ty = 'l', ylim = plim)
lines(rs0$x_grid, pred0[,4], lty = 'dashed')
lines(rs0$x_grid, pred0[,6], lty = 'dashed')

col1 = 'red'
lines(rs0$x_grid, pred1[,5], ty = 'l', col= col1)
lines(rs0$x_grid, pred1[,4], lty = 'dashed', col= col1)
lines(rs0$x_grid, pred1[,6], lty = 'dashed', col= col1)

col1 = 'blue'
lines(rs0$x_grid, pred2[,5], ty = 'l', col= col1)
lines(rs0$x_grid, pred2[,4], lty = 'dashed', col= col1)
lines(rs0$x_grid, pred2[,6], lty = 'dashed', col= col1)

intercept = rs2$pred[1,]
rbind(intercept, rs1$pred[1,])
slope = rs2$pred[2,]; slope

rs3 = NULL  
rs3 = BayesianPsplines:::smooth_inla3(x,group, y, Ntrials, family = fam, hyperB= hyperGumbel, diff.order = 3)
# intercept
intercept = rs3$pred[1,]
Pred.rw = rs3$pred[-1,]
x_grid = rs3$x_grid; irange = rs3$indices$grid#

y_pred = Pred.rw[irange,5] # median
y_lo =  Pred.rw[irange,4] # lower quantile
y_hi =  Pred.rw[irange,6] # upper quantile
y_pred_mean = Pred.rw[irange,"mean"] # mean on transformed scale
y_pred_sd = Pred.rw[irange,"sd"] # sd on transformed scale



plot(x_grid, backtransform(y_pred), type = "l", ylim= ylim, ylab= my.ylab,xlab = my.xlab , main = species)
lines(x_grid, backtransform(y_hi), lty="dashed")
lines(x_grid, backtransform(y_lo), lty="dashed")
# points(x, 100*y/Ntrials) # raw data

# from residuals on the transformed scale
yhat = transform(rs3$fitted[,"0.5quant"])
observed = transform(y/Ntrials) # on transformed scale
residual = observed - yhat
irangex = rs3$indices$x#
componentplus = residual + Pred.rw[irangex,"0.5quant"] 
points(x, backtransform(componentplus) ) # component-plus-residual plot

