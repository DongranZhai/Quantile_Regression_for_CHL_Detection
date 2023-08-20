# Note: Log transformation of CHL to compare CHL sensitivity about upper quantiles.
#       Panals
library(ggplot2)
library(raster)
library(oceanmap)
library(scales)
library(meboot) # bootstrap
library(boot) # bootstrap
library(zoo)
library(R.matlab)
########## Part I: log transform directly on ts ##########
# Results almost same
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/month_whitening_occci.Rdata")

month <- 1:298
n <- length(month)
mon_deszn <- as.vector(log(mon_deszn_trans+1)) # range: 0.025-0.057 to 0.024-0.056
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
tau <- c(0.01,0.05,seq(0.1, 0.9, by = 0.1),0.95,0.99)
ci_df <- matrix(data=NA,nrow=5,ncol=1)
bstar_levels <- matrix(data=NA,nrow=999,ncol=1)

for(i in 1:13){ # 13 quantile levels
  print(i)
  tick <- 2*i
  # trend
  fit <- rq(mon_deszn~month,tau=tau[i],na.action = na.omit)
  sd <- summary(fit,se='ker')$coefficients[4]
  
  # ci
  tr <- theta(mon_deszn,month,tau[i])
  coef.star <- bstar(y=mon_deszn,x=month,tau=tau[i],theta=tr)
  
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

save(ci_df,bstar_levels,month,file='globe_quants_trans_log_ci.Rdata')

##### Plot CI
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/globe_quants_trans_log_ci.Rdata")
# One of the CI method: create df with percentile CI
ci_df <- as.matrix(ci_df)
multi_rq_df <- matrix(data=NA,nrow=1,ncol=4)
for(i in 1:13){
  print(i)
  tick1 <- i*4-3
  tick2 <- i*4
  boot.percentile <- ci_df[3,c(tick1:tick2)]
  multi_rq_df <- rbind(multi_rq_df,boot.percentile)
  
}
multi_rq_df <- as.data.frame(multi_rq_df)
multi_rq_df <- multi_rq_df[-1,]
colnames(multi_rq_df) <- c('trend','sd','lower_ci','upper_ci')
quantiletemp <- round(c(1,5,seq(10, 90, by = 10),95,99),digits=2)
multi_rq_df$quantile <- quantiletemp

# plot II
plot <- ggplot(data=multi_rq_df,aes(x=quantile,y=trend),size=1) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=1.5, position=position_dodge(0.1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.8) +
  scale_x_continuous(breaks = seq(0,100,10)) + #,expand = c(0, 0)
  scale_y_continuous(breaks = seq(0,0.06,0.01)) +
  xlab("Quantile levels (%)") +
  ylab(expression(paste("Trend (10"^"-3", "·mg·m"^" -3","·yr"^"-1",")"))) +
  theme_bw()+
  theme(plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"),
        legend.key.size = unit(0.7, 'cm'))
ggsave("global_log_ci.png",width=7.5, height=3.44, dpi=300)

########## Part II: log transform from CHL ##########
load("/Volumes/Doris/OCCCI/occci_chl_v6.Rdata")
chl <- log(chl+1)
# keep NA values
temp <- readMat("/Volumes/Doris/Projects/sptmodel/Longhurst_180.mat") #used to identify longhursst regions
Longhurst <- temp$Longhurst
nlon <- length(lon)
nlat <- length(lat)
area <- sort(unique(Longhurst[!is.na(Longhurst)]))#list of longhurst areas to iterate through
area <- area[-c(1,2,4,6,7,8,10,11,12,13,17,18,19,22,25,26,27,29,30,31,34,35,36,40,43,41,42,47,50,53,54)] #remove coastal, polar, GoM, Archipelagic Deep basin, Mediterranean

for(i in 1:nlon){
  for(j in 1:nlat){
    if(Longhurst[i,j] %in% area){
      Longhurst[i,j] <- Longhurst[i,j]
    }else{
      Longhurst[i,j] <- NA
    }
  }
}

ind <- which(is.na(Longhurst))
##### ######### ######### Monthly global projection
monchl <- array(data=NA,dim=c(nlon,nlat,298))
for(i in 1:298){
  tempchl <- chl[,,i]
  tempchl[ind] <- NA
  monchl[,,i] <- tempchl
}

##### ######### time series
dim.chl <- dim(monchl)
time <- c(1:dim.chl[3])
monchl_ts <- array(data=NA,dim = c(length(time)))
for(i in time){
  print(i)
  montemp <- as.numeric(mean(monchl[,,i],na.rm=T))
  monchl_ts[i] <- montemp
}
# quant_ts <- data.frame(month=time,chl=monchl)
save(monchl,monchl_ts,lon,lat,file='monthly_chl_rm_coastal_log.Rdata')

# ts plot
# plot(time,mon_glavrg,type='l',col='orange',lwd=2.5,xaxt="n",
     # xlab='Time',ylab=expression(paste('CHL (mg m'^'-3',')')))
plot(time,monchl_ts,type='l',col='#7f8c8d',lwd=2.5,xaxt="n",
     xlab='Time',ylab=expression(paste('CHL (mg m'^'-3',')')))
axis(1, at = seq(1,298,by =24),labels=seq(1997,2022,by=2),las=1,col='black')

