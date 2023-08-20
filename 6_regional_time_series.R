# Notes: This script is about analyzing regional time series after averaging CHL.

library(quantreg)
library(raster)
library(oceanmap)

########## Part I. Regional averaging timeseries: North Pacific ##########
### get timeseries
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/np_fit_panel.Rdata")
month <- 1:304
np_ts <- array(data=NA,dim = c(length(month)))
for(i in month){
  print(i)
  montemp <- as.numeric(mean(np_deszn[,,i],na.rm=T))
  np_ts[i] <- montemp
}
np_ts <- as.vector(np_ts)
### quantile levels analysis
# linear function
lm_func <- function(month,y){
  fit <- lm(y ~ month,na.action = na.omit)
  fitted.value <- fit$fitted.values
  return(fitted.value)
}
# quantile function
qr_func <- function(month,y,z){
  fit <- rq(y ~ month,tau=z,na.action = na.omit)
  fitted.value <- fit$fitted.values
  return(fitted.value)
}

ts.lm <- lm_func(month,np_ts)
ts.01 <- qr_func(month,np_ts,0.01)
ts.05 <- qr_func(month,np_ts,0.05)
ts.10 <- qr_func(month,np_ts,0.10)
ts.50 <- qr_func(month,np_ts,0.50)
ts.90 <- qr_func(month,np_ts,0.90)
ts.95 <- qr_func(month,np_ts,0.95)
ts.99 <- qr_func(month,np_ts,0.99)

save(np_ts,ts.lm,ts.01,ts.05,ts.50,ts.95,ts.99,ts.10,ts.90,month,file='np_fit_ts.Rdata')

########## Part II. Regional averaging timeseries: North Atlantic ##########
### get timeseries
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/na_fit_panel.Rdata")
month <- 1:304
na_ts <- array(data=NA,dim = c(length(month)))
for(i in month){
  print(i)
  montemp <- as.numeric(mean(na_deszn[,,i],na.rm=T))
  na_ts[i] <- montemp
}
na_ts <- as.vector(na_ts)
### quantile levels analysis
# linear function
lm_func <- function(month,y){
  fit <- lm(y ~ month,na.action = na.omit)
  fitted.value <- fit$fitted.values
  return(fitted.value)
}
# quantile function
qr_func <- function(month,y,z){
  fit <- rq(y ~ month,tau=z,na.action = na.omit)
  fitted.value <- fit$fitted.values
  return(fitted.value)
}

ts.lm <- lm_func(month,na_ts)
ts.01 <- qr_func(month,na_ts,0.01)
ts.05 <- qr_func(month,na_ts,0.05)
ts.10 <- qr_func(month,na_ts,0.10)
ts.50 <- qr_func(month,na_ts,0.50)
ts.90 <- qr_func(month,na_ts,0.90)
ts.95 <- qr_func(month,na_ts,0.95)
ts.99 <- qr_func(month,na_ts,0.99)

save(na_ts,ts.lm,ts.01,ts.05,ts.50,ts.95,ts.99,ts.10,ts.90,month,file='na_fit_ts.Rdata')

########## Part III. Regional averaging timeseries: Equatorial Pacific ##########
### get timeseries
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/ep_fit_panel.Rdata")
month <- 1:304
ep_ts <- array(data=NA,dim = c(length(month)))
for(i in month){
  print(i)
  montemp <- as.numeric(mean(ep_deszn[,,i],na.rm=T))
  ep_ts[i] <- montemp
}
ep_ts <- as.vector(ep_ts)
### quantile levels analysis
# linear function
lm_func <- function(month,y){
  fit <- lm(y ~ month,na.action = na.omit)
  fitted.value <- fit$fitted.values
  return(fitted.value)
}
# quantile function
qr_func <- function(month,y,z){
  fit <- rq(y ~ month,tau=z,na.action = na.omit)
  fitted.value <- fit$fitted.values
  return(fitted.value)
}

ts.lm <- lm_func(month,ep_ts)
ts.01 <- qr_func(month,ep_ts,0.01)
ts.05 <- qr_func(month,ep_ts,0.05)
ts.10 <- qr_func(month,ep_ts,0.10)
ts.50 <- qr_func(month,ep_ts,0.50)
ts.90 <- qr_func(month,ep_ts,0.90)
ts.95 <- qr_func(month,ep_ts,0.95)
ts.99 <- qr_func(month,ep_ts,0.99)

save(ep_ts,ts.lm,ts.01,ts.05,ts.50,ts.95,ts.99,ts.10,ts.90,month,file='ep_fit_ts.Rdata')

########## Part IV. Regional averaging timeseries: Southern Ocean ##########
### get timeseries
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/so_fit_panel.Rdata")
month <- 1:304
so_ts <- array(data=NA,dim = c(length(month)))
for(i in month){
  print(i)
  montemp <- as.numeric(mean(so_deszn[,,i],na.rm=T))
  so_ts[i] <- montemp
}
so_ts <- as.vector(so_ts)
### quantile levels analysis
# linear function
lm_func <- function(month,y){
  fit <- lm(y ~ month,na.action = na.omit)
  fitted.value <- fit$fitted.values
  return(fitted.value)
}
# quantile function
qr_func <- function(month,y,z){
  fit <- rq(y ~ month,tau=z,na.action = na.omit)
  fitted.value <- fit$fitted.values
  return(fitted.value)
}

ts.lm <- lm_func(month,so_ts)
ts.01 <- qr_func(month,so_ts,0.01)
ts.05 <- qr_func(month,so_ts,0.05)
ts.10 <- qr_func(month,so_ts,0.10)
ts.50 <- qr_func(month,so_ts,0.50)
ts.90 <- qr_func(month,so_ts,0.90)
ts.95 <- qr_func(month,so_ts,0.95)
ts.99 <- qr_func(month,so_ts,0.99)

save(so_ts,ts.lm,ts.01,ts.05,ts.50,ts.95,ts.99,ts.10,ts.90,month,file='so_fit_ts.Rdata')

########## Part V. Regional averaging timeseries: North Pacific Subtropical ##########
### get timeseries
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/npo_fit_panel.Rdata")
month <- 1:304
npo_ts <- array(data=NA,dim = c(length(month)))
for(i in month){
  print(i)
  montemp <- as.numeric(mean(npo_deszn[,,i],na.rm=T))
  npo_ts[i] <- montemp
}
npo_ts <- as.vector(npo_ts)
### quantile levels analysis
# linear function
lm_func <- function(month,y){
  fit <- lm(y ~ month,na.action = na.omit)
  fitted.value <- fit$fitted.values
  return(fitted.value)
}
# quantile function
qr_func <- function(month,y,z){
  fit <- rq(y ~ month,tau=z,na.action = na.omit)
  fitted.value <- fit$fitted.values
  return(fitted.value)
}

ts.lm <- lm_func(month,npo_ts)
ts.01 <- qr_func(month,npo_ts,0.01)
ts.05 <- qr_func(month,npo_ts,0.05)
ts.10 <- qr_func(month,npo_ts,0.10)
ts.50 <- qr_func(month,npo_ts,0.50)
ts.90 <- qr_func(month,npo_ts,0.90)
ts.95 <- qr_func(month,npo_ts,0.95)
ts.99 <- qr_func(month,npo_ts,0.99)

save(npo_ts,ts.lm,ts.01,ts.05,ts.50,ts.95,ts.99,ts.10,ts.90,month,file='npo_fit_ts.Rdata')

########## Part VI. Regional averaging timeseries: North Atlantic Subtropical ##########
### get timeseries
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/nao_fit_panel.Rdata")
month <- 1:304
nao_ts <- array(data=NA,dim = c(length(month)))
for(i in month){
  print(i)
  montemp <- as.numeric(mean(nao_deszn[,,i],na.rm=T))
  nao_ts[i] <- montemp
}
nao_ts <- as.vector(nao_ts)
### quantile levels analysis
# linear function
lm_func <- function(month,y){
  fit <- lm(y ~ month,na.action = na.omit)
  fitted.value <- fit$fitted.values
  return(fitted.value)
}
# quantile function
qr_func <- function(month,y,z){
  fit <- rq(y ~ month,tau=z,na.action = na.omit)
  fitted.value <- fit$fitted.values
  return(fitted.value)
}

ts.lm <- lm_func(month,nao_ts)
ts.01 <- qr_func(month,nao_ts,0.01)
ts.05 <- qr_func(month,nao_ts,0.05)
ts.10 <- qr_func(month,nao_ts,0.10)
ts.50 <- qr_func(month,nao_ts,0.50)
ts.90 <- qr_func(month,nao_ts,0.90)
ts.95 <- qr_func(month,nao_ts,0.95)
ts.99 <- qr_func(month,nao_ts,0.99)

save(nao_ts,ts.lm,ts.01,ts.05,ts.50,ts.95,ts.99,ts.10,ts.90,month,file='nao_fit_ts.Rdata')




