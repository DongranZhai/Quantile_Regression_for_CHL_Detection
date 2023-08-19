# Notes:
# Part I: pre-whitening monthly data
# Part II: Use pre-whitening data to plot global time series
# Part III: Confidence interval of global trends

library(quantreg) # quantile regression
library(orcutt)
library(ggplot2)
library(raster)
library(oceanmap)
library(scales)
library(meboot) # bootstrap
library(boot) # bootstrap
########## Part I: pre-whitening monthly data ################
rm(list=ls())
# set working directory and import initial dataset
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
nlon <- length(lon)
nlat <- length(lat)
month <- 1:304
n <- length(month)
ntime <- length(month)
# rho
fit <- lm(mon_deszn ~ month,na.action = na.omit)
fit.cotf <- cochrane.orcutt(fit)
rho <- fit.cotf$rho
# rhohat = sum(fit$residuals[-1]*fit$residuals[-n])/sum(fit$residuals^2)
# mon_deszn
mon_deszn_trans <- array(NA,dim=c(304))
temp <- mon_deszn[-1]-mon_deszn[-ntime]*rho
mon_deszn_trans[1] <- temp[1]
mon_deszn_trans[2:ntime] <- temp

plot(month,mon_deszn_trans,type='l')
# monchl_deszn
y.trans <- array(NA,dim=c(nlon,nlat,ntime))
for(i in 1:nlon){
  for(j in 1:nlat){
    naNum <- sum(is.na(monchl_deszn[i,j,]))
    if (naNum <= 300){
      temp.1 <- monchl_deszn[i,j,]
      temp.2 <- temp.1[-1]-temp.1[-ntime]*rho
      y.trans[i,j,1] <- temp.2[1]
      y.trans[i,j,c(2:ntime)] <- temp.2
      
    }
  }
}

monchl_deszn_trans <- y.trans
save(mon_deszn_trans,monchl_deszn_trans,lon,lat,file = 'month_whitening_occci.Rdata')

########## Part II: Use pre-whitening data plot global ts ################
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
load("/Volumes/Doris/OCCCI/quantile_regression/month_whitening_occci.Rdata")
# Setting date format
Date <- seq(from = as.Date("1997/09/01",format="%Y/%m/%d"), 
            to = as.Date("2022/6/01",format="%Y/%m/%d"), 
            by = "1 month")
month <- 1:304
n <- length(month)
quan_level <- c(0.01,0.05,0.1,0.5,0.90,0.95,0.99)
data_ts <- data.frame(month=month,date=Date,fitted=mon_deszn_trans)
quantile_df <- as.data.frame(data_ts[,-c(1:2)])
colnames(quantile_df) <- c('fitted')

# quantile
for(i in c(1:7)){
  ind <- 304*i
  fit.trans <- rq(data_ts$fitted ~ month,tau=quan_level[i],na.action = na.omit)
  # fit.cotf <- cochrane.orcutt(fit)
  # rho <- fit.cotf$rho
  # y.trans <- mon_deszn[-1]-mon_deszn[-n]*rho
  # x.trans <- month[-1]-month[-n]*rho
  # fit.trans <- rq(y.trans ~ x.trans,tau=quan_level[i],na.action = na.omit)
  tr_trans <- as.vector(fit.trans$coefficients[2])
  # fitted.value <- as.data.frame(c(fit.trans$fitted.values[1],fit.trans$fitted.values))
  fitted.value <- as.data.frame(fit.trans$fitted.values)
  colnames(fitted.value) <- c('fitted')
  quantile_df <- rbind(quantile_df,fitted.value)
}
# linear
fit.trans <- lm(mon_deszn ~ month,na.action = na.omit)
# fit.cotf <- cochrane.orcutt(fit)
# rho <- fit.cotf$rho
# y.trans <- mon_deszn[-1]-mon_deszn[-n]*rho
# x.trans <- month[-1]-month[-n]*rho
# fit.trans <- lm(y.trans ~ x.trans,na.action = na.omit)
# fitted.value <- as.data.frame(c(fit.trans$fitted.values[1],fit.trans$fitted.values))
fitted.value <- as.data.frame(fit.trans$fitted.values)
colnames(fitted.value) <- c('fitted')
quantile_df <- rbind(quantile_df,fitted.value)

quantile_df$quantile <- rep(c('obs','01st','05th','10th','50th','90th','95th','99th','mean'),each=304)
quantile_df$month <- rep(month,time=9)
quantile_df$date <- rep(Date,time=9)

### plot
plot <- ggplot() +
  geom_line(data = data_ts,aes(x=date,y=mon_deszn),alpha=0.8,position = "identity",col='#bdc3c7')+
  geom_path(data=quantile_df,aes(x=date,y=fitted,group=quantile,colour=as.factor(quantile),size=quantile))+ 
  scale_colour_manual(values=c('#bdc3c7','#2e86de','#48dbfb','#c7ecee','#2ecc71','#fad390','#e67e22','#e74c3c','#6c5ce7'),
                      limits=c('obs','01st','05th','10th','50th','90th','95th','99th','mean'))+
  scale_size_manual(values = c('obs'=0.6,'01st'=1,'05th'=1,'10th'=1,
                               '50th'=1,'90th'=1,'95th'=1,'99th'=1,'mean'=1)) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year", limit = c(as.Date("1997-09-01"),as.Date("2022-06-01"))) + 
  scale_y_continuous(breaks = seq(0.03,0.06,0.01)) +
  coord_cartesian(xlim = c(as.Date("1997-06-01"),as.Date("2022-06-01")),
                  ylim = c(0.025,0.06),expand=F)+ 
  xlab("Time (yr)") +
  ylab(expression(paste("CHL  ","(mg·m"^"-3",")"))) +
  # ggtitle(paste0("1997-2020 Time Series Monthly Chlorophyll a data (.m)")) +
  guides(color=guide_legend(title = 'quantile levels',title.position = "top",title.hjust = .5, show.limits = T,
                            label.position = 'left', direction = 'verticle',ncol = 1, byrow = F),size=F,
         guide_colorbar(barwidth = unit(2, "cm"),barheight = unit(.5, "cm"))) +
  theme_bw()+
  theme(legend.position = 'right',plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"),
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.7, 'cm'),
        axis.text.x=element_text(angle=60, hjust=1))

ggsave("global_pre_ts.png",width=7.5, height=3.44, dpi=300)


############### Part III: Confidence interval of global trends ###############
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
# load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
load("/Volumes/Doris/OCCCI/quantile_regression/month_whitening_occci.Rdata")
month <- 1:304
n <- length(month)
mon_deszn <- as.vector(mon_deszn_trans) # range: 0.15-0.22 to 0.025-0.06
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

save(ci_df,bstar_levels,month,file='globe_quants_trans_ci.Rdata')

##### Plot CI ###########
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/globe_quants_trans_ci.Rdata")
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

# plot I
png(filename = 'globe_quantile_trans_ci.png',width=1400,height=700)
par(mar=c(4,4,2,2))
plot(multi_rq_df$quantile,multi_rq_df$trend,type='o',lty=1,pch=20,lwd=5,
     ylim=c(0,0.06),xlim=c(0,100),xaxt = "n",yaxt = "n",
     xlab=substitute(paste(bold('Quantile levels (%)'))),ylab=expression(bold(paste("Trend (10"^"-3", "·mg·m"^" -3","·yr"^"-1",")"))),
     # main=substitute(paste(bold('Globe'))),
     mgp = c(2, 0.5, 0),tck = 0.02)
axis(1, at = seq(0,100,by =10),labels = seq(0,100,by =10),tck = 0.02,mgp = c(2, 0.5, 0),font = 2)
axis(2, at = seq(0,0.06,by =0.01),labels = seq(0,0.06,by =0.01),tck = 0.02,mgp = c(2, 0.5, 0),font = 2)
arrows(x0=multi_rq_df$quantile, y0=multi_rq_df$lower_ci, x1=multi_rq_df$quantile, y1=multi_rq_df$upper_ci,
       code=3, angle=90, length=0.07, col="black", lwd=3)
abline(h=0,col='#e74c3c',lty=2,lwd=2.5)
dev.off()

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
ggsave("global_ci.png",width=7.5, height=3.44, dpi=300)

