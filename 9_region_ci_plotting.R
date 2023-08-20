# Note: This script is about plotting quantile levels (trend curves)
#       with confidence interval.
#       Two methods for plotting

library(dplyr)
library(ggplot2)
library(raster)
library(oceanmap)
library(scales)
library(viridis)
library(rgdal)

##### Part I: Confidence Interval of quantile levels #####
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
# par(mfrow=c(2,2))
# North Pacific 
load("/Volumes/Doris/OCCCI/quantile_regression/np_quants_ci.Rdata")
# One of the CI method: create df with percentile CI
ci_df <- as.matrix(ci_df)
multi_rq_df <- matrix(data=NA,nrow=1,ncol=4)
for(i in 1:11){
  print(i)
  tick1 <- i*4-3
  tick2 <- i*4
  boot.percentile <- ci_df[3,c(tick1:tick2)]
  multi_rq_df <- rbind(multi_rq_df,boot.percentile)
  
}
multi_rq_df <- as.data.frame(multi_rq_df)
multi_rq_df <- multi_rq_df[-1,]
colnames(multi_rq_df) <- c('trend','sd','lower_ci','upper_ci')
quantiletemp <- round(c(5,seq(10, 90, by = 10),95),digits=2)
multi_rq_df$quantile <- quantiletemp
range(multi_rq_df$lower_ci)
range(multi_rq_df$upper_ci)

# # plot
# png(filename = 'np_ci.png',width=1400,height=700)
# plot(multi_rq_df$quantile,multi_rq_df$trend,type='o',lty=1,pch=20,lwd=3.5,
#      ylim=c(-0.04,0.8),xlim=c(0,1),xaxt = "n",yaxt = "n",cex.main=2,
#      xlab=substitute(paste(bold('Quantile levels (%)'))),
#      ylab=expression(bold(paste("Trend (10"^"-3", "·mg·m"^" -3","·yr"^"-1",")"))),
#      main=substitute(paste(bold('North Pacific Subarctic Gyres Province'))),
#      mgp = c(2,0.5,0),tck = 0.02) # ylim c(0.16,0.86)
# axis(1, at = seq(0,1,by =0.1),labels = seq(0,1,by =0.1),
#      tck = 0.02,mgp = c(2,0.5,0),font = 2)
# # axis(2, at = seq(0.16,0.86,by =0.05),labels = seq(0.16,0.86,by =0.05),tck = 0.02,mgp = c(2, 0.5, 0))
# axis(2, at = seq(0,0.8,by =0.1),labels = seq(0,0.8,by =0.1),
#      tck = 0.02,mgp = c(2,0.5,0),font = 2)
# arrows(x0=multi_rq_df$quantile, y0=multi_rq_df$lower_ci, x1=multi_rq_df$quantile, y1=multi_rq_df$upper_ci,
#        code=3, angle=90, length=0.07, col="black", lwd=3)
# abline(h=0,col='#e74c3c',lty=2,lwd=2.5)
# dev.off()

# plot II
plot <- ggplot(data=multi_rq_df,aes(x=quantile,y=trend),size=1) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=1.5, position=position_dodge(0.1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.8) +
  scale_x_continuous(breaks = seq(0,100,10)) + #,expand = c(0, 0)
  scale_y_continuous(breaks = seq(-0.1,0.4,by =0.1),limits = c(-0.1,0.4)) +
  xlab("Quantile levels (%)") +
  ylab(expression(paste("Trend (10"^"-3", "·mg·m"^" -3","·yr"^"-1",")"))) +
  ggtitle(paste0("North Pacific Subarctic Gyres Province")) +
  theme_bw()+
  theme(plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"),
        legend.title = element_text(size = 10),plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.7, 'cm'))
ggsave("np_ci_new.png",width=7.5, height=3.44, dpi=300)

# North Atlantic
load("/Volumes/Doris/OCCCI/quantile_regression/na_quants_ci.Rdata")
# One of the CI method: create df with percentile CI
ci_df <- as.matrix(ci_df)
multi_rq_df <- matrix(data=NA,nrow=1,ncol=4)
for(i in 1:11){
  print(i)
  tick1 <- i*4-3
  tick2 <- i*4
  boot.percentile <- ci_df[3,c(tick1:tick2)]
  multi_rq_df <- rbind(multi_rq_df,boot.percentile)
  
}
multi_rq_df <- as.data.frame(multi_rq_df)
multi_rq_df <- multi_rq_df[-1,]
colnames(multi_rq_df) <- c('trend','sd','lower_ci','upper_ci')
quantiletemp <- round(c(5,seq(10, 90, by = 10),95),digits=2)
multi_rq_df$quantile <- quantiletemp
range(multi_rq_df$lower_ci)
range(multi_rq_df$upper_ci)
# # plot
# png(filename = 'na_ci.png',width=1400,height=700)
# plot(multi_rq_df$quantile,multi_rq_df$trend,type='o',lty=1,pch=20,lwd=3,
#      ylim=c(-0.1,0.8),xlim=c(0,1),xaxt = "n",yaxt = "n",cex.main=2,
#      xlab=substitute(paste(bold('Quantile levels (%)'))),
#      ylab=expression(bold(paste("Trend (10"^"-3", "·mg·m"^" -3","·yr"^"-1",")"))),
#      main=substitute(paste(bold('North Atlantic Drift Province'))),
#      mgp = c(2,0.5,0),tck = 0.02) # ylim c(0.16,0.86)
# axis(1, at = seq(0,1,by =0.1),labels = seq(0,1,by =0.1),
#      tck = 0.02,mgp = c(2,0.5,0),font = 2)
# # axis(2, at = seq(0.16,0.86,by =0.05),labels = seq(0.16,0.86,by =0.05),tck = 0.02,mgp = c(2, 0.5, 0))
# axis(2, at = seq(0,0.8,by =0.1),labels = seq(0,0.8,by =0.1),
#      tck = 0.02,mgp = c(2,0.5,0),font = 2)
# arrows(x0=multi_rq_df$quantile, y0=multi_rq_df$lower_ci, x1=multi_rq_df$quantile, y1=multi_rq_df$upper_ci,
#        code=3, angle=90, length=0.07, col="black", lwd=3)
# abline(h=0,col='#e74c3c',lty=2,lwd=2.5)
# dev.off()

# plot II
plot <- ggplot(data=multi_rq_df,aes(x=quantile,y=trend),size=1) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=1.5, position=position_dodge(0.1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.8) +
  scale_x_continuous(breaks = seq(0,100,10)) + #,expand = c(0, 0)
  scale_y_continuous(breaks = seq(0,0.4,by =0.1),limits = c(-0.1,0.4)) +
  xlab("Quantile levels (%)") +
  ylab(expression(paste("Trend (10"^"-3", "·mg·m"^" -3","·yr"^"-1",")"))) +
  ggtitle(paste0("North Atlantic Drift Province")) +
  theme_bw()+
  theme(plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"),
        legend.title = element_text(size = 10),plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.7, 'cm'))
ggsave("na_ci_new.png",width=7.5, height=3.44, dpi=300)

# Equatorial Pacific
load("/Volumes/Doris/OCCCI/quantile_regression/ep_quants_ci.Rdata")
# One of the CI method: create df with percentile CI
ci_df <- as.matrix(ci_df)
multi_rq_df <- matrix(data=NA,nrow=1,ncol=4)
for(i in 1:11){
  print(i)
  tick1 <- i*4-3
  tick2 <- i*4
  boot.percentile <- ci_df[3,c(tick1:tick2)]
  multi_rq_df <- rbind(multi_rq_df,boot.percentile)
  
}
multi_rq_df <- as.data.frame(multi_rq_df)
multi_rq_df <- multi_rq_df[-1,]
colnames(multi_rq_df) <- c('trend','sd','lower_ci','upper_ci')
quantiletemp <- round(c(5,seq(10, 90, by = 10),95),digits=2)
multi_rq_df$quantile <- quantiletemp
range(multi_rq_df$lower_ci)
range(multi_rq_df$upper_ci)
# # plot
# png(filename = 'ep_ci.png',width=1400,height=700)
# plot(multi_rq_df$quantile,multi_rq_df$trend,type='o',lty=1,pch=20,lwd=3,
#      ylim=c(-0.04,0.07),xlim=c(0,1),xaxt = "n",yaxt = "n",cex.main=2,
#      xlab=substitute(paste(bold('Quantile levels (%)'))),
#      ylab=expression(bold(paste("Trend (10"^"-3", "·mg·m"^" -3","·yr"^"-1",")"))),
#      main=substitute(paste(bold('Pacific Equatorial Province'))),
#      mgp = c(2,0.5,0),tck = 0.02) # ylim c(0.16,0.86)
# axis(1, at = seq(0,1,by =0.1),labels = seq(0,1,by =0.1),
#      tck = 0.02,mgp = c(2,0.5,0),font = 2)
# # axis(2, at = seq(0.16,0.86,by =0.05),labels = seq(0.16,0.86,by =0.05),tck = 0.02,mgp = c(2, 0.5, 0))
# axis(2, at = seq(-0.04,0.07,by =0.01),labels = seq(-0.04,0.07,by =0.01),
#      tck = 0.02,mgp = c(2,0.5,0),font = 2)
# arrows(x0=multi_rq_df$quantile, y0=multi_rq_df$lower_ci, x1=multi_rq_df$quantile, y1=multi_rq_df$upper_ci,
#        code=3, angle=90, length=0.07, col="black", lwd=3)
# abline(h=0,col='#e74c3c',lty=2,lwd=2.5)
# dev.off()

# plot II
plot <- ggplot(data=multi_rq_df,aes(x=quantile,y=trend),size=1) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=1.5, position=position_dodge(0.1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.8) +
  scale_x_continuous(breaks = seq(0,100,10)) + #,expand = c(0, 0)
  scale_y_continuous(breaks = seq(-0.04,0.04,by =0.01),limits = c(-0.04,0.04)) +
  xlab("Quantile levels (%)") +
  ylab(expression(paste("Trend (10"^"-3", "·mg·m"^" -3","·yr"^"-1",")"))) +
  ggtitle(paste0("Pacific Equatorial Province")) +
  theme_bw()+
  theme(plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"),
        legend.title = element_text(size = 10),plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.7, 'cm'))
ggsave("ep_ci_new.png",width=7.5, height=3.44, dpi=300)

# Southern Ocean
load("/Volumes/Doris/OCCCI/quantile_regression/so_quants_ci.Rdata")
# One of the CI method: create df with percentile CI
ci_df <- as.matrix(ci_df)
multi_rq_df <- matrix(data=NA,nrow=1,ncol=4)
for(i in 1:11){
  print(i)
  tick1 <- i*4-3
  tick2 <- i*4
  boot.percentile <- ci_df[3,c(tick1:tick2)]
  multi_rq_df <- rbind(multi_rq_df,boot.percentile)
  
}
multi_rq_df <- as.data.frame(multi_rq_df)
multi_rq_df <- multi_rq_df[-1,]
colnames(multi_rq_df) <- c('trend','sd','lower_ci','upper_ci')
quantiletemp <- round(c(5,seq(10, 90, by = 10),95),digits=2)
multi_rq_df$quantile <- quantiletemp
range(multi_rq_df$lower_ci)
range(multi_rq_df$upper_ci)
# # plot
# png(filename = 'so_ci.png',width=1400,height=700)
# plot(multi_rq_df$quantile,multi_rq_df$trend,type='o',lty=1,pch=20,lwd=3,
#      ylim=c(0,0.13),xlim=c(0,1),xaxt = "n",yaxt = "n",cex.main=2,
#      xlab=substitute(paste(bold('Quantile levels (%)'))),
#      ylab=expression(bold(paste("Trend (10"^"-3", "·mg·m"^" -3","·yr"^"-1",")"))),
#      main=substitute(paste(bold('Subantarctic Province'))),
#      mgp = c(2,0.5,0),tck = 0.02) # ylim c(0.16,0.86)
# axis(1, at = seq(0,1,by =0.1),labels = seq(0,1,by =0.1),
#      tck = 0.02,mgp = c(2,0.5,0),font = 2)
# # axis(2, at = seq(0.16,0.86,by =0.05),labels = seq(0.16,0.86,by =0.05),tck = 0.02,mgp = c(2, 0.5, 0))
# axis(2, at = seq(0,0.14,by =0.02),labels = seq(0,0.14,by =0.02),
#      tck = 0.02,mgp = c(2,0.5,0),font = 2)
# arrows(x0=multi_rq_df$quantile, y0=multi_rq_df$lower_ci, x1=multi_rq_df$quantile, y1=multi_rq_df$upper_ci,
#        code=3, angle=90, length=0.07, col="black", lwd=3)
# abline(h=0,col='#e74c3c',lty=2,lwd=2.5)
# dev.off()

# plot II
plot <- ggplot(data=multi_rq_df,aes(x=quantile,y=trend),size=1) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=1.5, position=position_dodge(0.1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.8) +
  scale_x_continuous(breaks = seq(0,100,10)) + #,expand = c(0, 0)
  scale_y_continuous(breaks = seq(0,0.1,by =0.02),limits = c(0,0.1)) +
  xlab("Quantile levels (%)") +
  ylab(expression(paste("Trend (10"^"-3", "·mg·m"^" -3","·yr"^"-1",")"))) +
  ggtitle(paste0("Subantarctic Province")) +
  theme_bw()+
  theme(plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"),
        legend.title = element_text(size = 10),plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.7, 'cm'))
ggsave("so_ci_new.png",width=7.5, height=3.44, dpi=300)

# North Pacific Oligotrophic
rm(list=ls())
load("/Volumes/Doris/OCCCI/quantile_regression/npo_quants_ci.Rdata")
# One of the CI method: create df with percentile CI
ci_df <- as.matrix(ci_df)
multi_rq_df <- matrix(data=NA,nrow=1,ncol=4)
for(i in 1:11){
  print(i)
  tick1 <- i*4-3
  tick2 <- i*4
  boot.percentile <- ci_df[3,c(tick1:tick2)]
  multi_rq_df <- rbind(multi_rq_df,boot.percentile)
  
}
multi_rq_df <- as.data.frame(multi_rq_df)
multi_rq_df <- multi_rq_df[-1,]
colnames(multi_rq_df) <- c('trend','sd','lower_ci','upper_ci')
quantiletemp <- round(c(5,seq(10, 90, by = 10),95),digits=2)
multi_rq_df$quantile <- quantiletemp
range(multi_rq_df$lower_ci)
range(multi_rq_df$upper_ci)
# # plot
# png(filename = 'npo_ci.png',width=1400,height=700)
# plot(multi_rq_df$quantile,multi_rq_df$trend,type='o',lty=1,pch=20,lwd=3.5,
#      ylim=c(-0.05,0.01),xlim=c(0,1),xaxt = "n",yaxt = "n",cex.main=2,
#      xlab=substitute(paste(bold('Quantile levels (%)'))),
#      ylab=expression(bold(paste("Trend (10"^"-3", "·mg·m"^" -3","·yr"^"-1",")"))),
#      main=substitute(paste(bold('North Pacific Subtropical Province'))),
#      mgp = c(2,0.5,0),tck = 0.02) # ylim c(0.16,0.86)
# axis(1, at = seq(0,1,by =0.1),labels = seq(0,1,by =0.1),
#      tck = 0.02,mgp = c(2,0.5,0),font = 2)
# # axis(2, at = seq(0.16,0.86,by =0.05),labels = seq(0.16,0.86,by =0.05),tck = 0.02,mgp = c(2, 0.5, 0))
# axis(2, at = seq(-0.05,0.01,by =0.02),labels = seq(-0.05,0.01,by =0.02),
#      tck = 0.02,mgp = c(2,0.5,0),font = 2)
# arrows(x0=multi_rq_df$quantile, y0=multi_rq_df$lower_ci, x1=multi_rq_df$quantile, y1=multi_rq_df$upper_ci,
#        code=3, angle=90, length=0.07, col="black", lwd=3)
# abline(h=0,col='#e74c3c',lty=2,lwd=2.5)
# dev.off()

# plot II
plot <- ggplot(data=multi_rq_df,aes(x=quantile,y=trend),size=1) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=1.5, position=position_dodge(0.1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.8) +
  scale_x_continuous(breaks = seq(0,100,10)) + #,expand = c(0, 0)
  scale_y_continuous(breaks = seq(-0.015,0.015,by =0.01),limits = c(-0.015,0.015)) +
  xlab("Quantile levels (%)") +
  ylab(expression(paste("Trend (10"^"-3", "·mg·m"^" -3","·yr"^"-1",")"))) +
  ggtitle(paste0("North Pacific Subtropical Gyre Province")) +
  theme_bw()+
  theme(plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"),
        legend.title = element_text(size = 10),plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.7, 'cm'))
ggsave("npo_ci_new.png",width=7.5, height=3.44, dpi=300)

# North Atlantic Oligotrophic
rm(list=ls())
load("/Volumes/Doris/OCCCI/quantile_regression/nao_quants_ci.Rdata")
# One of the CI method: create df with percentile CI
ci_df <- as.matrix(ci_df)
multi_rq_df <- matrix(data=NA,nrow=1,ncol=4)
for(i in 1:11){
  print(i)
  tick1 <- i*4-3
  tick2 <- i*4
  boot.percentile <- ci_df[3,c(tick1:tick2)]
  multi_rq_df <- rbind(multi_rq_df,boot.percentile)
  
}
multi_rq_df <- as.data.frame(multi_rq_df)
multi_rq_df <- multi_rq_df[-1,]
colnames(multi_rq_df) <- c('trend','sd','lower_ci','upper_ci')
quantiletemp <- round(c(5,seq(10, 90, by = 10),95),digits=2)
multi_rq_df$quantile <- quantiletemp
range(multi_rq_df$lower_ci)
range(multi_rq_df$upper_ci)
# # plot
# png(filename = 'nao_ci.png',width=1400,height=700)
# plot(multi_rq_df$quantile,multi_rq_df$trend,type='o',lty=1,pch=20,lwd=3.5,
#      ylim=c(-0.14,0.2),xlim=c(0,1),xaxt = "n",yaxt = "n",cex.main=2,
#      xlab=substitute(paste(bold('Quantile levels (%)'))),
#      ylab=expression(bold(paste("Trend (10"^"-3", "·mg·m"^" -3","·yr"^"-1",")"))),
#      main=substitute(paste(bold('North Atlantic Subtropical Province'))),
#      mgp = c(2,0.5,0),tck = 0.02) # ylim c(0.16,0.86)
# axis(1, at = seq(0,1,by =0.1),labels = seq(0,1,by =0.1),
#      tck = 0.02,mgp = c(2,0.5,0),font = 2)
# # axis(2, at = seq(0.16,0.86,by =0.05),labels = seq(0.16,0.86,by =0.05),tck = 0.02,mgp = c(2, 0.5, 0))
# axis(2, at = seq(-0.14,0.2,by =0.05),labels = seq(-0.14,0.2,by =0.05),
#      tck = 0.02,mgp = c(2,0.5,0),font = 2)
# arrows(x0=multi_rq_df$quantile, y0=multi_rq_df$lower_ci, x1=multi_rq_df$quantile, y1=multi_rq_df$upper_ci,
#        code=3, angle=90, length=0.07, col="black", lwd=3)
# abline(h=0,col='#e74c3c',lty=2,lwd=2.5)
# dev.off()

# plot II
plot <- ggplot(data=multi_rq_df,aes(x=quantile,y=trend),size=1) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=1.5, position=position_dodge(0.1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.8) +
  scale_x_continuous(breaks = seq(0,100,10)) + #,expand = c(0, 0)
  scale_y_continuous(breaks = seq(-0.02,0.03,by =0.01),limits = c(-0.02,0.035)) +
  xlab("Quantile levels (%)") +
  ylab(expression(paste("Trend (10"^"-3", "·mg·m"^" -3","·yr"^"-1",")"))) +
  ggtitle(paste0("North Atlantic Subtropical Gyre Province")) +
  theme_bw()+
  theme(plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"),
        legend.title = element_text(size = 10),plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.7, 'cm'))
ggsave("nao_ci_new.png",width=7.5, height=3.44, dpi=300)

