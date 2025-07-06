rm(list=ls())

source('R/fe_var.R')

dat <- read.csv('VAR3_levels_fixedeffects_IVprctiles.csv', header = F)
dat <- dat[,c(1,2,seq(7,NCOL(dat)))]
fe <- dat[['V1']]
lags <- 5
Y <- as.matrix(dat[,c(3,4,5)])
colnames(Y) <- c('household_debt_to_lagged_gdp', 'nonfinancial_debt_to_lagged_gdp', 'real_log_gdp')
X <- as.matrix(dat[,6:NCOL(dat)])
X <- X[,1:(lags*ncol(Y))]
colnames(X) <- paste0(c('household_debt_to_lagged_gdp', 'nonfinancial_debt_to_lagged_gdp', 'real_log_gdp'),
                      '_lag', rep(1:lags, each = ncol(Y)))
res <- fe_var(Y, X, fe, lags, seed = 1,
              perm = matrix(c(0,0,1,0,1,0,1,0,0), 3,3,byrow=T) # with real log GDP ordered first, followed
              # by nonfinancial debt to lagged GDP, and household debt to lagged GDP
)
res$bet0
res$bet
irf(res$bet0, res$sigma0, perm = matrix(c(0,0,1,0,1,0,1,0,0), 3,3,byrow=T), normalize = T, irf_hor = 11, lags = 5)
res$virf$mean
plot_irf(res$virf, 1)
plot_irf(res$virf, 2)
plot_irf(res$virf, 3)
