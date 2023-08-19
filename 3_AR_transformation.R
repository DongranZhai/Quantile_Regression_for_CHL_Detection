# Notes: This script is about global projection in quantile regression
#        after taking autocorrelation of residual into account.
#        Using Cochrane-Orcutt and Hildreth-Lu procedures.
###  Part I. Cochrane-Orcutt procedure to transform 
###  Part II. Hildreth Lu estimation (similar with Cochrane) 
### The same script is also used to do the log-transformation of CHL

library(rgdal)
library(quantreg)
library(boot)
library(stats)
library(orcutt)

############### ############### ############### ############### ############### 
###  Part I. Cochrane-Orcutt procedure to transform ####
###### 1% ######
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression/Cochrane_Orcutt")
# setwd("/Volumes/Doris/OCCCI/quantile_regression/log_trans")
load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
# mon_deszn <- as.vector(log(mon_deszn+1)) 

tr <- array(data=NA,dim=c(length(lon),length(lat)))
tr_trans <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
pValue <- array(data=NA,dim=c(length(lon),length(lat)))
se <- array(data=NA,dim=c(length(lon),length(lat)))
dw.cotf <- array(data=NA,dim=c(length(lon),length(lat)))

z=0.01
month <- seq(1,304,1)
n <- length(month)
tick_na <- 0
tick <- 0 

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(monchl_deszn[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit <- rq(monchl_deszn[i,j,] ~ month,tau=z,na.action = na.omit)
      tr[i,j] <- as.vector(fit$coefficients[2])
      # b0 <- fit$coefficients[1]
      # b1 <- fit$coefficients[2]
      fit.cotf <- cochrane.orcutt(fit)
      dw.cotf[i,j] <- fit.cotf$DW[3]
      
      rho <- fit.cotf$rho
      # rhohat = sum(fit$residuals[-1]*fit$residuals[-n])/sum(fit$residuals^2)
      y.trans <- monchl_deszn[i,j,-1]-monchl_deszn[i,j,-n]*rho
      # x.trans <- month[-1]-month[-n]*rho
      fit.trans <- rq(y.trans ~ month[-1],tau=z,na.action = na.omit) # x.trans
      
      tr_trans[i,j] <- as.vector(fit.trans$coefficients[2])
      pValue[i,j] <- summary(fit.trans,se='ker')$coefficients[8]
      se[i,j] <- summary(fit.trans,se='ker')$coefficients[4]
      
      tick <- tick + 1
      
    } else {
      tick_na <- tick_na + 1
    }
  }
}

savename <- paste0("ar_cochrane_01.Rdata")
save(tr,tr_trans,dw.cotf,pValue,se,lon,lat,file = savename)

###### 5% ######
rm(list=ls())
# setwd("/Volumes/Doris/OCCCI/quantile_regression/global/Cochrane_Orcutt")
setwd("/Volumes/Doris/OCCCI/quantile_regression/log_trans")
load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
mon_deszn <- as.vector(log(mon_deszn+1)) 

tr <- array(data=NA,dim=c(length(lon),length(lat)))
tr_trans <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
pValue <- array(data=NA,dim=c(length(lon),length(lat)))
se <- array(data=NA,dim=c(length(lon),length(lat)))
dw.cotf <- array(data=NA,dim=c(length(lon),length(lat)))

z=0.05
month <- seq(1,304,1)
n <- length(month)
tick_na <- 0
tick <- 0 

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(monchl_deszn[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit <- rq(monchl_deszn[i,j,] ~ month,tau=z,na.action = na.omit)
      tr[i,j] <- as.vector(fit$coefficients[2])
      # b0 <- fit$coefficients[1]
      # b1 <- fit$coefficients[2]
      fit.cotf <- cochrane.orcutt(fit)
      dw.cotf[i,j] <- fit.cotf$DW[3]
      
      rho <- fit.cotf$rho
      y.trans <- monchl_deszn[i,j,-1]-monchl_deszn[i,j,-n]*rho
      x.trans <- month[-1]-month[-n]*rho
      fit.trans <- rq(y.trans ~ x.trans,tau=z,na.action = na.omit)
      
      tr_trans[i,j] <- as.vector(fit.trans$coefficients[2])
      pValue[i,j] <- summary(fit.trans,se='ker')$coefficients[8]
      se[i,j] <- summary(fit.trans,se='ker')$coefficients[4]
      
      tick <- tick + 1
      # Box test
      
    } else {
      tick_na <- tick_na + 1
    }
  }
}

savename <- paste0("ar_cochrane_05.Rdata")
save(tr,tr_trans,dw.cotf,pValue,se,lon,lat,file = savename)

###### 10% ######
rm(list=ls())
# setwd("/Volumes/Doris/OCCCI/quantile_regression/global/Cochrane_Orcutt")
setwd("/Volumes/Doris/OCCCI/quantile_regression/log_trans")
load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
mon_deszn <- as.vector(log(mon_deszn+1)) 

tr <- array(data=NA,dim=c(length(lon),length(lat)))
tr_trans <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
pValue <- array(data=NA,dim=c(length(lon),length(lat)))
se <- array(data=NA,dim=c(length(lon),length(lat)))
dw.cotf <- array(data=NA,dim=c(length(lon),length(lat)))

z=0.10
month <- seq(1,304,1)
n <- length(month)
tick_na <- 0
tick <- 0 

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(monchl_deszn[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit <- rq(monchl_deszn[i,j,] ~ month,tau=z,na.action = na.omit)
      tr[i,j] <- as.vector(fit$coefficients[2])
      # b0 <- fit$coefficients[1]
      # b1 <- fit$coefficients[2]
      fit.cotf <- cochrane.orcutt(fit)
      dw.cotf[i,j] <- fit.cotf$DW[3]
      
      rho <- fit.cotf$rho
      y.trans <- monchl_deszn[i,j,-1]-monchl_deszn[i,j,-n]*rho
      x.trans <- month[-1]-month[-n]*rho
      fit.trans <- rq(y.trans ~ x.trans,tau=z,na.action = na.omit)
      
      tr_trans[i,j] <- as.vector(fit.trans$coefficients[2])
      pValue[i,j] <- summary(fit.trans,se='ker')$coefficients[8]
      se[i,j] <- summary(fit.trans,se='ker')$coefficients[4]
      
      tick <- tick + 1
      # Box test
      
    } else {
      tick_na <- tick_na + 1
    }
  }
}

savename <- paste0("ar_cochrane_10.Rdata")
save(tr,tr_trans,dw.cotf,pValue,se,lon,lat,file = savename)

###### 50% ######
rm(list=ls())
# setwd("/Volumes/Doris/OCCCI/quantile_regression/global/Cochrane_Orcutt")
setwd("/Volumes/Doris/OCCCI/quantile_regression/log_trans")
load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
mon_deszn <- as.vector(log(mon_deszn+1)) 

tr <- array(data=NA,dim=c(length(lon),length(lat)))
tr_trans <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
pValue <- array(data=NA,dim=c(length(lon),length(lat)))
se <- array(data=NA,dim=c(length(lon),length(lat)))
dw.cotf <- array(data=NA,dim=c(length(lon),length(lat)))

z=0.50
month <- seq(1,304,1)
n <- length(month)
tick_na <- 0
tick <- 0 

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(monchl_deszn[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit <- rq(monchl_deszn[i,j,] ~ month,tau=z,na.action = na.omit)
      tr[i,j] <- as.vector(fit$coefficients[2])
      # b0 <- fit$coefficients[1]
      # b1 <- fit$coefficients[2]
      fit.cotf <- cochrane.orcutt(fit)
      dw.cotf[i,j] <- fit.cotf$DW[3]
      
      rho <- fit.cotf$rho
      y.trans <- monchl_deszn[i,j,-1]-monchl_deszn[i,j,-n]*rho
      x.trans <- month[-1]-month[-n]*rho
      fit.trans <- rq(y.trans ~ x.trans,tau=z,na.action = na.omit)
      
      tr_trans[i,j] <- as.vector(fit.trans$coefficients[2])
      pValue[i,j] <- summary(fit.trans,se='ker')$coefficients[8]
      se[i,j] <- summary(fit.trans,se='ker')$coefficients[4]
      
      tick <- tick + 1
      # Box test
      
    } else {
      tick_na <- tick_na + 1 
    }
  }
}

savename <- paste0("ar_cochrane_50.Rdata")
save(tr,tr_trans,dw.cotf,pValue,se,lon,lat,file = savename)

###### 90% ######
rm(list=ls())
# setwd("/Volumes/Doris/OCCCI/quantile_regression/global/Cochrane_Orcutt")
setwd("/Volumes/Doris/OCCCI/quantile_regression/log_trans")
load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
mon_deszn <- as.vector(log(mon_deszn+1)) 

tr <- array(data=NA,dim=c(length(lon),length(lat)))
tr_trans <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
pValue <- array(data=NA,dim=c(length(lon),length(lat)))
se <- array(data=NA,dim=c(length(lon),length(lat)))
dw.cotf <- array(data=NA,dim=c(length(lon),length(lat)))

z=0.90
month <- seq(1,304,1)
n <- length(month)
tick_na <- 0
tick <- 0 

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(monchl_deszn[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit <- rq(monchl_deszn[i,j,] ~ month,tau=z,na.action = na.omit)
      tr[i,j] <- as.vector(fit$coefficients[2])
      # b0 <- fit$coefficients[1]
      # b1 <- fit$coefficients[2]
      fit.cotf <- cochrane.orcutt(fit)
      dw.cotf[i,j] <- fit.cotf$DW[3]
      
      rho <- fit.cotf$rho
      y.trans <- monchl_deszn[i,j,-1]-monchl_deszn[i,j,-n]*rho
      x.trans <- month[-1]-month[-n]*rho
      fit.trans <- rq(y.trans ~ x.trans,tau=z,na.action = na.omit)
      
      tr_trans[i,j] <- as.vector(fit.trans$coefficients[2])
      pValue[i,j] <- summary(fit.trans,se='ker')$coefficients[8]
      se[i,j] <- summary(fit.trans,se='ker')$coefficients[4]
      
      tick <- tick + 1
      # Box test
      
    } else {
      tick_na <- tick_na + 1
    }
  }
}

savename <- paste0("ar_cochrane_90.Rdata")
save(tr,tr_trans,dw.cotf,pValue,se,lon,lat,file = savename)

###### 95% ######
rm(list=ls())
# setwd("/Volumes/Doris/OCCCI/quantile_regression/global/Cochrane_Orcutt")
setwd("/Volumes/Doris/OCCCI/quantile_regression/log_trans")
load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
mon_deszn <- as.vector(log(mon_deszn+1)) 

tr <- array(data=NA,dim=c(length(lon),length(lat)))
tr_trans <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
pValue <- array(data=NA,dim=c(length(lon),length(lat)))
se <- array(data=NA,dim=c(length(lon),length(lat)))
dw.cotf <- array(data=NA,dim=c(length(lon),length(lat)))

z=0.95
month <- seq(1,304,1)
n <- length(month)
tick_na <- 0
tick <- 0 

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(monchl_deszn[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit <- rq(monchl_deszn[i,j,] ~ month,tau=z,na.action = na.omit)
      tr[i,j] <- as.vector(fit$coefficients[2])
      # b0 <- fit$coefficients[1]
      # b1 <- fit$coefficients[2]
      fit.cotf <- cochrane.orcutt(fit)
      dw.cotf[i,j] <- fit.cotf$DW[3]
      
      rho <- fit.cotf$rho
      y.trans <- monchl_deszn[i,j,-1]-monchl_deszn[i,j,-n]*rho
      x.trans <- month[-1]-month[-n]*rho
      fit.trans <- rq(y.trans ~ x.trans,tau=z,na.action = na.omit)
      
      tr_trans[i,j] <- as.vector(fit.trans$coefficients[2])
      pValue[i,j] <- summary(fit.trans,se='ker')$coefficients[8]
      se[i,j] <- summary(fit.trans,se='ker')$coefficients[4]
      
      tick <- tick + 1
      # Box test
      
    } else {
      tick_na <- tick_na + 1 
    }
  }
}

savename <- paste0("ar_cochrane_95.Rdata")
save(tr,tr_trans,dw.cotf,pValue,se,lon,lat,file = savename)

###### 99% ######
rm(list=ls())
# setwd("/Volumes/Doris/OCCCI/quantile_regression/global/Cochrane_Orcutt")
setwd("/Volumes/Doris/OCCCI/quantile_regression/log_trans")
load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
mon_deszn <- as.vector(log(mon_deszn+1)) 

tr <- array(data=NA,dim=c(length(lon),length(lat)))
tr_trans <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
pValue <- array(data=NA,dim=c(length(lon),length(lat)))
se <- array(data=NA,dim=c(length(lon),length(lat)))
dw.cotf <- array(data=NA,dim=c(length(lon),length(lat)))

z=0.99
month <- seq(1,304,1)
n <- length(month)
tick_na <- 0
tick <- 0 

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(monchl_deszn[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit <- rq(monchl_deszn[i,j,] ~ month,tau=z,na.action = na.omit)
      tr[i,j] <- as.vector(fit$coefficients[2])
      # b0 <- fit$coefficients[1]
      # b1 <- fit$coefficients[2]
      fit.cotf <- cochrane.orcutt(fit)
      dw.cotf[i,j] <- fit.cotf$DW[3]
      
      rho <- fit.cotf$rho
      y.trans <- monchl_deszn[i,j,-1]-monchl_deszn[i,j,-n]*rho
      x.trans <- month[-1]-month[-n]*rho
      fit.trans <- rq(y.trans ~ x.trans,tau=z,na.action = na.omit)
      
      tr_trans[i,j] <- as.vector(fit.trans$coefficients[2])
      pValue[i,j] <- summary(fit.trans,se='ker')$coefficients[8]
      se[i,j] <- summary(fit.trans,se='ker')$coefficients[4]
      
      tick <- tick + 1
      # Box test
      
    } else {
      tick_na <- tick_na + 1
    }
  }
}

savename <- paste0("ar_cochrane_99.Rdata")
save(tr,tr_trans,dw.cotf,pValue,se,lon,lat,file = savename)

###### simple linear model ######
rm(list=ls())
# setwd("/Volumes/Doris/OCCCI/quantile_regression/global/Cochrane_Orcutt")
setwd("/Volumes/Doris/OCCCI/quantile_regression/log_trans")
load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
mon_deszn <- as.vector(log(mon_deszn+1)) 

month <- 1:304
n <- length(month)
tr <- array(data=NA,dim=c(length(lon),length(lat))) #trend
tr_trans <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
pValue <- array(data=NA,dim=c(length(lon),length(lat)))
se <- array(data=NA,dim=c(length(lon),length(lat)))
dw.cotf <- array(data=NA,dim=c(length(lon),length(lat)))

tick_na <- 0 
tick <- 0 

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(monchl_deszn[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit <- lm(monchl_deszn[i,j,] ~ month,na.action = na.omit)
      tr[i,j] <- as.vector(fit$coefficients[2])
      
      fit.cotf <- cochrane.orcutt(fit)
      dw.cotf[i,j] <- fit.cotf$DW[3]
      
      rho <- fit.cotf$rho
      y.trans <- monchl_deszn[i,j,-1]-monchl_deszn[i,j,-n]*rho
      x.trans <- month[-1]-month[-n]*rho
      fit.trans <- lm(y.trans ~ x.trans,na.action = na.omit)
      
      tr_trans[i,j] <- as.vector(fit.trans$coefficients[2])
      pValue[i,j] <- summary(fit.trans)$coefficients[8]
      se[i,j] <- summary(fit.trans,se='ker')$coefficients[4]
      
      tick <- tick + 1
    } else {
      tick_na <- tick_na + 1
    }
  }
}
savename <- paste0("ar_cochrane_lm.Rdata")
save(tr,tr_trans,dw.cotf,pValue,se,lon,lat,file = savename)

############### ############### ############### ############### ############### 
###  Part II. Hildreth Lu estimation (similar with Cochrane) ####
###### 1% ######
rm(list=ls())
# setwd("/Volumes/Doris/OCCCI/quantile_regression/global/Cochrane_Orcutt")
# load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
setwd("/home/dzhai/quantile_regression")
load("/home/dzhai/quantile_regression/monthly_chl_rm_coastal_occci.Rdata")

tr <- array(data=NA,dim=c(length(lon),length(lat)))
tr_trans <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
pValue <- array(data=NA,dim=c(length(lon),length(lat)))
se <- array(data=NA,dim=c(length(lon),length(lat)))

# Create hildreth function
hildreth.lu.func <- function(r,z,model){
  x <- model.matrix(model)[,-1]
  y <- model.response(model.frame(model))
  n <- length(y)
  t <- 2:n
  y <- y[-1]-r*y[-n]
  x <- x[-1]-r*x[-n]
  new.model <- rq(y~x,tau=z,na.action = na.omit)
  return(new.model)
}

z=0.01
month <- seq(1,304,1)
n <- length(month)
tick_na <- 0
tick <- 0 
r <- seq(0.1,0.9, by= 0.01)

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(monchl_deszn[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit <- rq(monchl_deszn[i,j,] ~ month,tau=z,na.action = na.omit)
      tr[i,j] <- as.vector(fit$coefficients[2])
      
      # Hildreth iteration to find the best rho
      hl.tab <- data.frame("rho" = r, 
                           "SSE" = sapply(r,function(k){sum((fitted(hildreth.lu.func(k,z,fit))-hildreth.lu.func(k,z,fit)$y)^2)}))
      optrho <- which.min(round(hl.tab,4)[,2])
      rho <- round(hl.tab, 4)[optrho,1]
      
      # Hildreth transform original model using that rho
      fit.hildreth <- hildreth.lu.func(rho,z,fit)
      tr_trans[i,j] <- as.vector(fit.hildreth$coefficients[2])
      pValue[i,j] <- summary(fit.hildreth,se='ker')$coefficients[8]
      se[i,j] <- summary(fit.hildreth,se='ker')$coefficients[4]
      
      tick <- tick + 1
      
    } else {
      tick_na <- tick_na + 1 
    }
  }
}

savename <- paste0("ar_hildreth_01.Rdata")
save(tr,tr_trans,pValue,se,lon,lat,file = savename)

###### 5% ######
rm(list=ls())
# setwd("/Volumes/Doris/OCCCI/quantile_regression/global/Cochrane_Orcutt")
# load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
setwd("/home/dzhai/quantile_regression")
load("/home/dzhai/quantile_regression/monthly_chl_rm_coastal_occci.Rdata")

tr <- array(data=NA,dim=c(length(lon),length(lat)))
tr_trans <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
pValue <- array(data=NA,dim=c(length(lon),length(lat)))
se <- array(data=NA,dim=c(length(lon),length(lat)))

# Create hildreth function
hildreth.lu.func <- function(r,z,model){
  x <- model.matrix(model)[,-1]
  y <- model.response(model.frame(model))
  n <- length(y)
  t <- 2:n
  y <- y[-1]-r*y[-n]
  x <- x[-1]-r*x[-n]
  new.model <- rq(y~x,tau=z,na.action = na.omit)
  return(new.model)
}

z=0.05
month <- seq(1,304,1)
n <- length(month)
tick_na <- 0
tick <- 0 
r <- seq(0.1,0.9, by= 0.01)

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(monchl_deszn[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit <- rq(monchl_deszn[i,j,] ~ month,tau=z,na.action = na.omit)
      tr[i,j] <- as.vector(fit$coefficients[2])
      
      # Hildreth iteration to find the best rho
      hl.tab <- data.frame("rho" = r, 
                           "SSE" = sapply(r,function(k){sum((fitted(hildreth.lu.func(k,z,fit))-hildreth.lu.func(k,z,fit)$y)^2)}))
      optrho <- which.min(round(hl.tab,4)[,2])
      rho <- round(hl.tab, 4)[optrho,1]
      
      # Hildreth transform original model using that rho
      fit.hildreth <- hildreth.lu.func(rho,z,fit)
      tr_trans[i,j] <- as.vector(fit.hildreth$coefficients[2])
      pValue[i,j] <- summary(fit.hildreth,se='ker')$coefficients[8]
      se[i,j] <- summary(fit.hildreth,se='ker')$coefficients[4]
      
      tick <- tick + 1
      
    } else {
      tick_na <- tick_na + 1 
    }
  }
}

savename <- paste0("ar_hildreth_05.Rdata")
save(tr,tr_trans,pValue,se,lon,lat,file = savename)

###### 10% ######
rm(list=ls())
# setwd("/Volumes/Doris/OCCCI/quantile_regression/global/Cochrane_Orcutt")
# load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
setwd("/home/dzhai/quantile_regression")
load("/home/dzhai/quantile_regression/monthly_chl_rm_coastal_occci.Rdata")

tr <- array(data=NA,dim=c(length(lon),length(lat)))
tr_trans <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
pValue <- array(data=NA,dim=c(length(lon),length(lat)))
se <- array(data=NA,dim=c(length(lon),length(lat)))

# Create hildreth function
hildreth.lu.func <- function(r,z,model){
  x <- model.matrix(model)[,-1]
  y <- model.response(model.frame(model))
  n <- length(y)
  t <- 2:n
  y <- y[-1]-r*y[-n]
  x <- x[-1]-r*x[-n]
  new.model <- rq(y~x,tau=z,na.action = na.omit)
  return(new.model)
}

z=0.10
month <- seq(1,304,1)
n <- length(month)
tick_na <- 0
tick <- 0 
r <- seq(0.1,0.9, by= 0.01)

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(monchl_deszn[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit <- rq(monchl_deszn[i,j,] ~ month,tau=z,na.action = na.omit)
      tr[i,j] <- as.vector(fit$coefficients[2])
      
      # Hildreth iteration to find the best rho
      hl.tab <- data.frame("rho" = r, 
                           "SSE" = sapply(r,function(k){sum((fitted(hildreth.lu.func(k,z,fit))-hildreth.lu.func(k,z,fit)$y)^2)}))
      optrho <- which.min(round(hl.tab,4)[,2])
      rho <- round(hl.tab, 4)[optrho,1]
      
      # Hildreth transform original model using that rho
      fit.hildreth <- hildreth.lu.func(rho,z,fit)
      tr_trans[i,j] <- as.vector(fit.hildreth$coefficients[2])
      pValue[i,j] <- summary(fit.hildreth,se='ker')$coefficients[8]
      se[i,j] <- summary(fit.hildreth,se='ker')$coefficients[4]
      
      tick <- tick + 1
      
    } else {
      tick_na <- tick_na + 1 
    }
  }
}

savename <- paste0("ar_hildreth_10.Rdata")
save(tr,tr_trans,pValue,se,lon,lat,file = savename)

###### 50% ######
rm(list=ls())
# setwd("/Volumes/Doris/OCCCI/quantile_regression/global/Cochrane_Orcutt")
# load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
setwd("/home/dzhai/quantile_regression")
load("/home/dzhai/quantile_regression/monthly_chl_rm_coastal_occci.Rdata")

tr <- array(data=NA,dim=c(length(lon),length(lat)))
tr_trans <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
pValue <- array(data=NA,dim=c(length(lon),length(lat)))
se <- array(data=NA,dim=c(length(lon),length(lat)))

# Create hildreth function
hildreth.lu.func <- function(r,z,model){
  x <- model.matrix(model)[,-1]
  y <- model.response(model.frame(model))
  n <- length(y)
  t <- 2:n
  y <- y[-1]-r*y[-n]
  x <- x[-1]-r*x[-n]
  new.model <- rq(y~x,tau=z,na.action = na.omit)
  return(new.model)
}

z=0.50
month <- seq(1,304,1)
n <- length(month)
tick_na <- 0
tick <- 0 
r <- seq(0.1,0.9, by= 0.01)

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(monchl_deszn[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit <- rq(monchl_deszn[i,j,] ~ month,tau=z,na.action = na.omit)
      tr[i,j] <- as.vector(fit$coefficients[2])
      
      # Hildreth iteration to find the best rho
      hl.tab <- data.frame("rho" = r, 
                           "SSE" = sapply(r,function(k){sum((fitted(hildreth.lu.func(k,z,fit))-hildreth.lu.func(k,z,fit)$y)^2)}))
      optrho <- which.min(round(hl.tab,4)[,2])
      rho <- round(hl.tab, 4)[optrho,1]
      
      # Hildreth transform original model using that rho
      fit.hildreth <- hildreth.lu.func(rho,z,fit)
      tr_trans[i,j] <- as.vector(fit.hildreth$coefficients[2])
      pValue[i,j] <- summary(fit.hildreth,se='ker')$coefficients[8]
      se[i,j] <- summary(fit.hildreth,se='ker')$coefficients[4]
      
      tick <- tick + 1
      
    } else {
      tick_na <- tick_na + 1 
    }
  }
}

savename <- paste0("ar_hildreth_50.Rdata")
save(tr,tr_trans,pValue,se,lon,lat,file = savename)

###### 90% ######
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression/global/Cochrane_Orcutt")
load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
# setwd("/home/dzhai/quantile_regression")
# load("/home/dzhai/quantile_regression/monthly_chl_rm_coastal_occci.Rdata")

tr <- array(data=NA,dim=c(length(lon),length(lat)))
tr_trans <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
pValue <- array(data=NA,dim=c(length(lon),length(lat)))
se <- array(data=NA,dim=c(length(lon),length(lat)))

# Create hildreth function
hildreth.lu.func <- function(r,z,model){
  x <- model.matrix(model)[,-1]
  y <- model.response(model.frame(model))
  n <- length(y)
  t <- 2:n
  y <- y[-1]-r*y[-n]
  x <- x[-1]-r*x[-n]
  new.model <- rq(y~x,tau=z,na.action = na.omit)
  return(new.model)
}

z=0.90
month <- seq(1,304,1)
n <- length(month)
tick_na <- 0
tick <- 0 
r <- seq(0.1,0.9, by= 0.01)

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(monchl_deszn[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit <- rq(monchl_deszn[i,j,] ~ month,tau=z,na.action = na.omit)
      tr[i,j] <- as.vector(fit$coefficients[2])
      
      # Hildreth iteration to find the best rho
      hl.tab <- data.frame("rho" = r, 
                           "SSE" = sapply(r,function(k){sum((fitted(hildreth.lu.func(k,z,fit))-hildreth.lu.func(k,z,fit)$y)^2)}))
      optrho <- which.min(round(hl.tab,4)[,2])
      rho <- round(hl.tab, 4)[optrho,1]
      
      # Hildreth transform original model using that rho
      fit.hildreth <- hildreth.lu.func(rho,z,fit)
      tr_trans[i,j] <- as.vector(fit.hildreth$coefficients[2])
      pValue[i,j] <- summary(fit.hildreth,se='ker')$coefficients[8]
      se[i,j] <- summary(fit.hildreth,se='ker')$coefficients[4]
      
      tick <- tick + 1
      
    } else {
      tick_na <- tick_na + 1 
    }
  }
}

savename <- paste0("ar_hildreth_90.Rdata")
save(tr,tr_trans,pValue,se,lon,lat,file = savename)

###### 95% ######
rm(list=ls())
# setwd("/Volumes/Doris/OCCCI/quantile_regression/global/Cochrane_Orcutt")
# load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
setwd("/home/dzhai/quantile_regression")
load("/home/dzhai/quantile_regression/monthly_chl_rm_coastal_occci.Rdata")

tr <- array(data=NA,dim=c(length(lon),length(lat)))
tr_trans <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
pValue <- array(data=NA,dim=c(length(lon),length(lat)))
se <- array(data=NA,dim=c(length(lon),length(lat)))

# Create hildreth function
hildreth.lu.func <- function(r,z,model){
  x <- model.matrix(model)[,-1]
  y <- model.response(model.frame(model))
  n <- length(y)
  t <- 2:n
  y <- y[-1]-r*y[-n]
  x <- x[-1]-r*x[-n]
  new.model <- rq(y~x,tau=z,na.action = na.omit)
  return(new.model)
}

z=0.95
month <- seq(1,304,1)
n <- length(month)
tick_na <- 0
tick <- 0 
r <- seq(0.1,0.9, by= 0.01)

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(monchl_deszn[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit <- rq(monchl_deszn[i,j,] ~ month,tau=z,na.action = na.omit)
      tr[i,j] <- as.vector(fit$coefficients[2])
      
      # Hildreth iteration to find the best rho
      hl.tab <- data.frame("rho" = r, 
                           "SSE" = sapply(r,function(k){sum((fitted(hildreth.lu.func(k,z,fit))-hildreth.lu.func(k,z,fit)$y)^2)}))
      optrho <- which.min(round(hl.tab,4)[,2])
      rho <- round(hl.tab, 4)[optrho,1]
      
      # Hildreth transform original model using that rho
      fit.hildreth <- hildreth.lu.func(rho,z,fit)
      tr_trans[i,j] <- as.vector(fit.hildreth$coefficients[2])
      pValue[i,j] <- summary(fit.hildreth,se='ker')$coefficients[8]
      se[i,j] <- summary(fit.hildreth,se='ker')$coefficients[4]
      
      tick <- tick + 1
      
    } else {
      tick_na <- tick_na + 1 
    }
  }
}

savename <- paste0("ar_hildreth_95.Rdata")
save(tr,tr_trans,pValue,se,lon,lat,file = savename)

###### 99% ######
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression/global/Hildreth_lu")
load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
# setwd("/home/dzhai/quantile_regression/")
# load("/home/dzhai/quantile_regression/monthly_chl_rm_coastal_occci.Rdata")

tr <- array(data=NA,dim=c(length(lon),length(lat)))
tr_trans <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
pValue <- array(data=NA,dim=c(length(lon),length(lat)))
se <- array(data=NA,dim=c(length(lon),length(lat)))

# Create hildreth function
hildreth.lu.func <- function(r,z,model){
  x <- model.matrix(model)[,-1]
  y <- model.response(model.frame(model))
  n <- length(y)
  t <- 2:n
  y <- y[-1]-r*y[-n]
  x <- x[-1]-r*x[-n]
  new.model <- rq(y~x,tau=z,na.action = na.omit)
  return(new.model)
}

z=0.99
month <- seq(1,304,1)
n <- length(month)
tick_na <- 0
tick <- 0 
r <- seq(0.1,0.9, by= 0.01)

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(monchl_deszn[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit <- rq(monchl_deszn[i,j,] ~ month,tau=z,na.action = na.omit)
      tr[i,j] <- as.vector(fit$coefficients[2])
      
      # Hildreth iteration to find the best rho
      hl.tab <- data.frame("rho" = r, 
                           "SSE" = sapply(r,function(k){sum((fitted(hildreth.lu.func(k,z,fit))-hildreth.lu.func(k,z,fit)$y)^2)}))
      optrho <- which.min(round(hl.tab,4)[,2])
      rho <- round(hl.tab, 4)[optrho,1]
      
      # Hildreth transform original model using that rho
      fit.hildreth <- hildreth.lu.func(rho,z,fit)
      tr_trans[i,j] <- as.vector(fit.hildreth$coefficients[2])
      pValue[i,j] <- summary(fit.hildreth,se='ker')$coefficients[8]
      se[i,j] <- summary(fit.hildreth,se='ker')$coefficients[4]
      
      tick <- tick + 1
      
    } else {
      tick_na <- tick_na + 1 
    }
  }
}

savename <- paste0("ar_hildreth_99.Rdata")
save(tr,tr_trans,pValue,se,lon,lat,file = savename)

###### simple linear model ######
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression/global/Hildreth_lu")
load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
# setwd("/home/dzhai/quantile_regression")
# load("/home/dzhai/quantile_regression/monthly_chl_rm_coastal_occci.Rdata")

tr <- array(data=NA,dim=c(length(lon),length(lat)))
tr_trans <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
pValue <- array(data=NA,dim=c(length(lon),length(lat)))
se <- array(data=NA,dim=c(length(lon),length(lat)))

# Create hildreth function
hildreth.lu.func <- function(r,model){
  x <- model.matrix(model)[,-1]
  y <- model.response(model.frame(model))
  n <- length(y)
  t <- 2:n
  y <- y[-1]-r*y[-n]
  x <- x[-1]-r*x[-n]
  new.model <- lm(y~x,na.action = na.omit)
  return(new.model)
}

month <- seq(1,304,1)
n <- length(month)
tick_na <- 0
tick <- 0 
r <- seq(0.1,0.9, by= 0.01)

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(monchl_deszn[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit <- lm(monchl_deszn[i,j,] ~ month,na.action = na.omit)
      tr[i,j] <- as.vector(fit$coefficients[2])
      
      # Hildreth iteration to find the best rho
      hl.tab <- data.frame("rho" = r, 
                           "SSE" = sapply(r,function(k){sum((fitted(hildreth.lu.func(k,fit))-hildreth.lu.func(k,fit)$y)^2)}))
      optrho <- which.min(round(hl.tab,4)[,2])
      rho <- round(hl.tab, 4)[optrho,1]
      
      # Hildreth transform original model using that rho
      fit.hildreth <- hildreth.lu.func(rho,fit)
      tr_trans[i,j] <- as.vector(fit.hildreth$coefficients[2])
      pValue[i,j] <- summary(fit.hildreth,se='ker')$coefficients[8]
      se[i,j] <- summary(fit.hildreth,se='ker')$coefficients[4]
      
      tick <- tick + 1
      
    } else {
      tick_na <- tick_na + 1 
    }
  }
}

savename <- paste0("ar_hildreth_lm.Rdata")
save(tr,tr_trans,pValue,se,lon,lat,file = savename)

