# Note: This script is about getting confidence intervals of quantile
#       levels.
# Reference:
#       The meboot R Package: Maximum Entropy Bootstrap for Time Series
#       boot package

library(quantreg) # quantile regression
library(meboot) # bootstrap
library(boot)
library(boot) # bootstrap
library(hdrcde) # highest density regions
library(aod) # wald test

##### Part I: Confidence Interval of quantile levels ---- North Pacific #####
# TWO (GOOD!!!)
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression/")
load("/Volumes/Doris/OCCCI/quantile_regression/np_fit_ts.Rdata")
n <- length(month)

# Functions
# quantile regression coefficient
theta <- function(y, x, tau){
  fit <- rq(y~x,tau=tau,na.action = na.omit)
  thet <- as.numeric(fit$coefficients[2])
  return(thet)
}

# generates a large number of bootstrap parameter estimates
bstar <- function(y, x, tau, theta, level = 0.95, bigJ = 999, seed = 123){
  set.seed(123)
  semy <- meboot(x = y, reps = bigJ)$ensemble
  semx <- meboot(x = x, reps = bigJ)$ensemble
  n <- NROW(y)
  m <- length(theta)
  if(m!=1){
    stop("too many parameters in theta")
  }else{
    bb <- matrix(NA, bigJ) 
    for(j in 1:bigJ){
      yy <- semy[,j]
      xx <- semx[,j] 
      bb[j] <- theta(yy, xx, tau) 
    }
  }
  return(bb)
}

# matrix ts
tau <- c(0.05,seq(0.1, 0.9, by = 0.1),0.95)
ci_df <- matrix(data=NA,nrow=5,ncol=1)
bstar_levels <- matrix(data=NA,nrow=999,ncol=1)

for(i in 1:11){ # 11 quantile levels
  print(i)
  tick <- 2*i
  # trend
  fit <- rq(np_ts~month,tau=tau[i],na.action = na.omit)
  sd <- summary(fit,se='ker')$coefficients[4]
  
  # ci
  tr <- theta(np_ts,month,tau[i])
  coef.star <- bstar(y=np_ts,x=month,tau=tau[i],theta=tr)
  
  ci.m1 <- quantile(coef.star, c(0.025, 0.975), type = 8)
  ci.m2 <- null.ci(coef.star)
  out <- list(t = coef.star, t0 = tr, var.t0 = sd^2, R = 999)
  class(out) <- "boot"
  
  boot.percentile <- boot.ci(out, type = "perc")$percent[4:5] # two sides: (0.025,0.975)
  boot.norm <- boot.ci(out, type = "norm")$normal[2:3] # 0.95
  boot.basic <- boot.ci(out, type = "basic")$basic[4:5] # two sides: (0.025,0.975)
  temp <- rbind(ci.m1, ci.m2, boot.percentile, boot.norm, boot.basic)
  temp_df <- data.frame(tr,sd,temp)*1000
  quant_label <- c(paste0(tau[i],'_25'),paste0(tau[i],'_975'))
  colnames(temp_df) <- c(paste0(tau[i],'_tr'),paste0(tau[i],'_sd'),paste0(tau[i],'_25'),paste0(tau[i],'_975'))
  # dataframe of quantile levels with confidence interval
  ci_df <- cbind(ci_df,temp_df)
  bstar_levels <- cbind(bstar_levels,coef.star)
}
ci_df <- ci_df[,-1]
bstar_levels <- bstar_levels[,-1]

save(ci_df,bstar_levels,month,file='np_quants_ci.Rdata')
# # Find highest density regions with graphics for theta_hat under the null.
# tau <- c(0.01,0.05,seq(0.1, 0.9, by = 0.1),0.95,0.99)
# jpeg(filename = 'HDR_np.png',width=2430,height=2052,res=144)
# par(mfrow=c(4,3))
# for(j in c(1,2,4:13)){
#   print(j)
#   # manually enhance 1000
#   temp <- bstar_levels[,j]*1000
#   hdr.den(temp,main=paste0('HDR of ',tau[j],' quantile'))
# }
# # dev.off()

##### Part II: Confidence Interval of quantile levels ---- North Atlantic #####
# TWO (GOOD!!!)
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/na_fit_ts.Rdata")
n <- length(month)

# Functions
# quantile regression coefficient
theta <- function(y, x, tau){
  fit <- rq(y~x,tau=tau,na.action = na.omit)
  thet <- as.numeric(fit$coefficients[2])
  return(thet)
}

# generates a large number of bootstrap parameter estimates
bstar <- function(y, x, tau, theta, level = 0.95, bigJ = 999, seed = 123){
  set.seed(123)
  semy <- meboot(x = y, reps = bigJ)$ensemble
  semx <- meboot(x = x, reps = bigJ)$ensemble
  n <- NROW(y)
  m <- length(theta)
  if(m!=1){
    stop("too many parameters in theta")
  }else{
    bb <- matrix(NA, bigJ) 
    for(j in 1:bigJ){
      yy <- semy[,j]
      xx <- semx[,j] 
      bb[j] <- theta(yy, xx, tau) 
    }
  }
  return(bb)
}

# matrix ts
tau <- c(0.05,seq(0.1, 0.9, by = 0.1),0.95)
ci_df <- matrix(data=NA,nrow=5,ncol=1)
bstar_levels <- matrix(data=NA,nrow=999,ncol=1)
for(i in 1:11){ # 11 quantile levels
  print(i)
  tick <- 2*i
  # trend
  fit <- rq(na_ts~month,tau=tau[i],na.action = na.omit)
  sd <- summary(fit,se='ker')$coefficients[4]
  
  # ci
  tr <- theta(na_ts,month,tau[i])
  coef.star <- bstar(y=na_ts,x=month,tau=tau[i],theta=tr)
  
  ci.m1 <- quantile(coef.star, c(0.025, 0.975), type = 8)
  ci.m2 <- null.ci(coef.star)
  out <- list(t = coef.star, t0 = tr, var.t0 = sd^2, R = 999)
  class(out) <- "boot"
  
  boot.percentile <- boot.ci(out, type = "perc")$percent[4:5] # two sides: (0.025,0.975)
  boot.norm <- boot.ci(out, type = "norm")$normal[2:3] # 0.95
  boot.basic <- boot.ci(out, type = "basic")$basic[4:5] # two sides: (0.025,0.975)
  temp <- rbind(ci.m1, ci.m2, boot.percentile, boot.norm, boot.basic)
  temp_df <- data.frame(tr,sd,temp)*1000
  quant_label <- c(paste0(tau[i],'_25'),paste0(tau[i],'_975'))
  colnames(temp_df) <- c(paste0(tau[i],'_tr'),paste0(tau[i],'_sd'),paste0(tau[i],'_25'),paste0(tau[i],'_975'))
  # dataframe of quantile levels with confidence interval
  ci_df <- cbind(ci_df,temp_df)
  bstar_levels <- cbind(bstar_levels,coef.star)
}
ci_df <- ci_df[,-1]
bstar_levels <- bstar_levels[,-1]

save(ci_df,bstar_levels,month,file='na_quants_ci.Rdata')
# # Find highest density regions with graphics for theta_hat under the null.
# tau <- c(0.01,0.05,seq(0.1, 0.9, by = 0.1),0.95,0.99)
# jpeg(filename = 'HDR_all_na.png',width=2430,height=2052,res=144)
# par(mfrow=c(4,3))
# for(j in c(1,2,4:13)){
#   print(j)
#   # manually enhance 1000
#   temp <- bstar_levels[,j]*1000
#   hdr.den(temp,main=paste0('HDR of ',tau[j],' quantile'))
# }
# dev.off()

##### Part III: Confidence Interval of quantile levels ---- Equatorial Pacific #####
# TWO (GOOD!!!)
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/ep_fit_ts.Rdata")
n <- length(month)

# Functions
# quantile regression coefficient
theta <- function(y, x, tau){
  fit <- rq(y~x,tau=tau,na.action = na.omit)
  thet <- as.numeric(fit$coefficients[2])
  return(thet)
}

# generates a large number of bootstrap parameter estimates
bstar <- function(y, x, tau, theta, level = 0.95, bigJ = 999, seed = 123){
  set.seed(123)
  semy <- meboot(x = y, reps = bigJ)$ensemble
  semx <- meboot(x = x, reps = bigJ)$ensemble
  n <- NROW(y)
  m <- length(theta)
  if(m!=1){
    stop("too many parameters in theta")
  }else{
    bb <- matrix(NA, bigJ) 
    for(j in 1:bigJ){
      yy <- semy[,j]
      xx <- semx[,j] 
      bb[j] <- theta(yy, xx, tau) 
    }
  }
  return(bb)
}

# matrix ts
tau <- c(0.05,seq(0.1, 0.9, by = 0.1),0.95)
ci_df <- matrix(data=NA,nrow=5,ncol=1)
bstar_levels <- matrix(data=NA,nrow=999,ncol=1)

for(i in 1:11){ # 11 quantile levels
  print(i)
  tick <- 2*i
  # trend
  fit <- rq(ep_ts~month,tau=tau[i],na.action = na.omit)
  sd <- summary(fit,se='ker')$coefficients[4]
  
  # ci
  tr <- theta(ep_ts,month,tau[i])
  coef.star <- bstar(y=ep_ts,x=month,tau=tau[i],theta=tr)
  
  ci.m1 <- quantile(coef.star, c(0.025, 0.975), type = 8)
  ci.m2 <- null.ci(coef.star)
  out <- list(t = coef.star, t0 = tr, var.t0 = sd^2, R = 999)
  class(out) <- "boot"
  
  boot.percentile <- boot.ci(out, type = "perc")$percent[4:5] # two sides: (0.025,0.975)
  boot.norm <- boot.ci(out, type = "norm")$normal[2:3] # 0.95
  boot.basic <- boot.ci(out, type = "basic")$basic[4:5] # two sides: (0.025,0.975)
  temp <- rbind(ci.m1, ci.m2, boot.percentile, boot.norm, boot.basic)
  temp_df <- data.frame(tr,sd,temp)*1000
  quant_label <- c(paste0(tau[i],'_25'),paste0(tau[i],'_975'))
  colnames(temp_df) <- c(paste0(tau[i],'_tr'),paste0(tau[i],'_sd'),paste0(tau[i],'_25'),paste0(tau[i],'_975'))
  # dataframe of quantile levels with confidence interval
  ci_df <- cbind(ci_df,temp_df)
  bstar_levels <- cbind(bstar_levels,coef.star)
}
ci_df <- ci_df[,-1]
bstar_levels <- bstar_levels[,-1]

save(ci_df,bstar_levels,month,file='ep_quants_ci.Rdata')
# Find highest density regions with graphics for theta_hat under the null.
# tau <- c(0.01,0.05,seq(0.1, 0.9, by = 0.1),0.95,0.99)
# jpeg(filename = 'HDR_all_ep.png',width=2430,height=2052,res=144)
# par(mfrow=c(4,3))
# for(j in c(1,2,4:13)){
#   print(j)
#   # manually enhance 1000
#   temp <- bstar_levels[,j]*1000
#   hdr.den(temp,main=paste0('HDR of ',tau[j],' quantile'))
# }
# dev.off()

##### Part IV: Confidence Interval of quantile levels ---- Southern Ocean #####
# TWO (GOOD!!!)
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/so_fit_ts.Rdata")
n <- length(month)

# Functions
# quantile regression coefficient
theta <- function(y, x, tau){
  fit <- rq(y~x,tau=tau,na.action = na.omit)
  thet <- as.numeric(fit$coefficients[2])
  return(thet)
}

# generates a large number of bootstrap parameter estimates
bstar <- function(y, x, tau, theta, level = 0.95, bigJ = 999, seed = 123){
  set.seed(123)
  semy <- meboot(x = y, reps = bigJ)$ensemble
  semx <- meboot(x = x, reps = bigJ)$ensemble
  n <- NROW(y)
  m <- length(theta)
  if(m!=1){
    stop("too many parameters in theta")
  }else{
    bb <- matrix(NA, bigJ) 
    for(j in 1:bigJ){
      yy <- semy[,j]
      xx <- semx[,j] 
      bb[j] <- theta(yy, xx, tau) 
    }
  }
  return(bb)
}

# matrix ts
tau <- c(0.05,seq(0.1, 0.9, by = 0.1),0.95)
ci_df <- matrix(data=NA,nrow=5,ncol=1)
bstar_levels <- matrix(data=NA,nrow=999,ncol=1)

for(i in 1:11){ # 11 quantile levels
  print(i)
  tick <- 2*i
  # trend
  fit <- rq(so_ts~month,tau=tau[i],na.action = na.omit)
  sd <- summary(fit,se='ker')$coefficients[4]
  
  # ci
  tr <- theta(so_ts,month,tau[i])
  coef.star <- bstar(y=so_ts,x=month,tau=tau[i],theta=tr)
  
  ci.m1 <- quantile(coef.star, c(0.025, 0.975), type = 8)
  ci.m2 <- null.ci(coef.star)
  out <- list(t = coef.star, t0 = tr, var.t0 = sd^2, R = 999)
  class(out) <- "boot"
  
  boot.percentile <- boot.ci(out, type = "perc")$percent[4:5] # two sides: (0.025,0.975)
  boot.norm <- boot.ci(out, type = "norm")$normal[2:3] # 0.95
  boot.basic <- boot.ci(out, type = "basic")$basic[4:5] # two sides: (0.025,0.975)
  temp <- rbind(ci.m1, ci.m2, boot.percentile, boot.norm, boot.basic)
  temp_df <- data.frame(tr,sd,temp)*1000
  quant_label <- c(paste0(tau[i],'_25'),paste0(tau[i],'_975'))
  colnames(temp_df) <- c(paste0(tau[i],'_tr'),paste0(tau[i],'_sd'),paste0(tau[i],'_25'),paste0(tau[i],'_975'))
  # dataframe of quantile levels with confidence interval
  ci_df <- cbind(ci_df,temp_df)
  bstar_levels <- cbind(bstar_levels,coef.star)
}
ci_df <- ci_df[,-1]
bstar_levels <- bstar_levels[,-1]

save(ci_df,bstar_levels,month,file='so_quants_ci.Rdata')
# # Find highest density regions with graphics for theta_hat under the null.
# tau <- c(0.01,0.05,seq(0.1, 0.9, by = 0.1),0.95,0.99)
# jpeg(filename = 'HDR_all_so.png',width=2430,height=2052,res=144)
# par(mfrow=c(4,3))
# for(j in c(1,2,4:13)){
#   print(j)
#   # manually enhance 1000
#   temp <- bstar_levels[,j]*1000
#   hdr.den(temp,main=paste0('HDR of ',tau[j],' quantile'))
# }
# dev.off()

##### Part V: Confidence Interval of quantile levels ---- North Pacific Oligotrophic#####
# TWO (GOOD!!!)
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/npo_fit_ts.Rdata")
n <- length(month)

# Functions
# quantile regression coefficient
theta <- function(y, x, tau){
  fit <- rq(y~x,tau=tau,na.action = na.omit)
  thet <- as.numeric(fit$coefficients[2])
  return(thet)
}

# generates a large number of bootstrap parameter estimates
bstar <- function(y, x, tau, theta, level = 0.95, bigJ = 999, seed = 123){
  set.seed(123)
  semy <- meboot(x = y, reps = bigJ)$ensemble
  semx <- meboot(x = x, reps = bigJ)$ensemble
  n <- NROW(y)
  m <- length(theta)
  if(m!=1){
    stop("too many parameters in theta")
  }else{
    bb <- matrix(NA, bigJ) 
    for(j in 1:bigJ){
      yy <- semy[,j]
      xx <- semx[,j] 
      bb[j] <- theta(yy, xx, tau) 
    }
  }
  return(bb)
}

# matrix ts
tau <- c(0.05,seq(0.1, 0.9, by = 0.1),0.95)
ci_df <- matrix(data=NA,nrow=5,ncol=1)
bstar_levels <- matrix(data=NA,nrow=999,ncol=1)

for(i in 1:11){ # 11 quantile levels
  print(i)
  tick <- 2*i
  # trend
  fit <- rq(npo_ts~month,tau=tau[i],na.action = na.omit)
  sd <- summary(fit,se='ker')$coefficients[4]
  
  # ci
  tr <- theta(npo_ts,month,tau[i])
  coef.star <- bstar(y=npo_ts,x=month,tau=tau[i],theta=tr)
  
  ci.m1 <- quantile(coef.star, c(0.025, 0.975), type = 8)
  ci.m2 <- null.ci(coef.star)
  out <- list(t = coef.star, t0 = tr, var.t0 = sd^2, R = 999)
  class(out) <- "boot"
  
  boot.percentile <- boot.ci(out, type = "perc")$percent[4:5] # two sides: (0.025,0.975)
  boot.norm <- boot.ci(out, type = "norm")$normal[2:3] # 0.95
  boot.basic <- boot.ci(out, type = "basic")$basic[4:5] # two sides: (0.025,0.975)
  temp <- rbind(ci.m1, ci.m2, boot.percentile, boot.norm, boot.basic)
  temp_df <- data.frame(tr,sd,temp)*1000
  quant_label <- c(paste0(tau[i],'_25'),paste0(tau[i],'_975'))
  colnames(temp_df) <- c(paste0(tau[i],'_tr'),paste0(tau[i],'_sd'),paste0(tau[i],'_25'),paste0(tau[i],'_975'))
  # dataframe of quantile levels with confidence interval
  ci_df <- cbind(ci_df,temp_df)
  bstar_levels <- cbind(bstar_levels,coef.star)
}
ci_df <- ci_df[,-1]
bstar_levels <- bstar_levels[,-1]

save(ci_df,bstar_levels,month,file='npo_quants_ci.Rdata')

##### Part V: Confidence Interval of quantile levels ---- North Atlantic Oligotrophic#####
# TWO (GOOD!!!)
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/nao_fit_ts.Rdata")
n <- length(month)

# Functions
# quantile regression coefficient
theta <- function(y, x, tau){
  fit <- rq(y~x,tau=tau,na.action = na.omit)
  thet <- as.numeric(fit$coefficients[2])
  return(thet)
}

# generates a large number of bootstrap parameter estimates
bstar <- function(y, x, tau, theta, level = 0.95, bigJ = 999, seed = 123){
  set.seed(123)
  semy <- meboot(x = y, reps = bigJ)$ensemble
  semx <- meboot(x = x, reps = bigJ)$ensemble
  n <- NROW(y)
  m <- length(theta)
  if(m!=1){
    stop("too many parameters in theta")
  }else{
    bb <- matrix(NA, bigJ) 
    for(j in 1:bigJ){
      yy <- semy[,j]
      xx <- semx[,j] 
      bb[j] <- theta(yy, xx, tau) 
    }
  }
  return(bb)
}

# matrix ts
tau <- c(0.05,seq(0.1, 0.9, by = 0.1),0.95)
ci_df <- matrix(data=NA,nrow=5,ncol=1)
bstar_levels <- matrix(data=NA,nrow=999,ncol=1)

for(i in 1:11){ # 11 quantile levels
  print(i)
  tick <- 2*i
  # trend
  fit <- rq(nao_ts~month,tau=tau[i],na.action = na.omit)
  sd <- summary(fit,se='ker')$coefficients[4]
  
  # ci
  tr <- theta(nao_ts,month,tau[i])
  coef.star <- bstar(y=nao_ts,x=month,tau=tau[i],theta=tr)
  
  ci.m1 <- quantile(coef.star, c(0.025, 0.975), type = 8)
  ci.m2 <- null.ci(coef.star)
  out <- list(t = coef.star, t0 = tr, var.t0 = sd^2, R = 999)
  class(out) <- "boot"
  
  boot.percentile <- boot.ci(out, type = "perc")$percent[4:5] # two sides: (0.025,0.975)
  boot.norm <- boot.ci(out, type = "norm")$normal[2:3] # 0.95
  boot.basic <- boot.ci(out, type = "basic")$basic[4:5] # two sides: (0.025,0.975)
  temp <- rbind(ci.m1, ci.m2, boot.percentile, boot.norm, boot.basic)
  temp_df <- data.frame(tr,sd,temp)*1000
  quant_label <- c(paste0(tau[i],'_25'),paste0(tau[i],'_975'))
  colnames(temp_df) <- c(paste0(tau[i],'_tr'),paste0(tau[i],'_sd'),paste0(tau[i],'_25'),paste0(tau[i],'_975'))
  # dataframe of quantile levels with confidence interval
  ci_df <- cbind(ci_df,temp_df)
  bstar_levels <- cbind(bstar_levels,coef.star)
}
ci_df <- ci_df[,-1]
bstar_levels <- bstar_levels[,-1]

save(ci_df,bstar_levels,month,file='nao_quants_ci.Rdata')
