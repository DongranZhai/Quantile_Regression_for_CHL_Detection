# Note: regional time-series
library(dplyr)
library(ggplot2)
library(raster)
library(oceanmap)
library(scales)
library(viridis)
library(rgdal)
library(mapproj)
library(munsell)
library(RColorBrewer)
library(quantreg)
library(meboot)
library(abind)

########## Part I. Timeseries with multiple quantile levels ##########
# NP
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/np_fit_ts.Rdata")
quantile_df <- as.data.frame(abind(np_ts,ts.05,ts.10,ts.50,
                                   ts.90,ts.95,ts.lm))
colnames(quantile_df) <- c('fitted')
Date <- seq(from = as.Date("1997/09/01",format="%Y/%m/%d"), 
            to = as.Date("2022/12/01",format="%Y/%m/%d"), 
            by = "1 month")
quantile_df$date <- rep(Date,time=7)
quantile_df$quantile <- rep(c('obs','05th','10th','50th','90th','95th','mean'),each=304)
range(quantile_df$fitted)
### plot
plot <- ggplot() +
  # geom_line(data = data_ts,aes(x=date,y=mon_deszn),alpha=0.8,position = "identity",col='#bdc3c7')+
  geom_line(data=quantile_df,aes(x=date,y=fitted,group=quantile,colour=as.factor(quantile),size=quantile))+ 
  scale_colour_manual(values=c('#bdc3c7','#2e86de','#c7ecee','#2ecc71','#fad390','#e74c3c','#6c5ce7'), # '#48dbfb',,'#e67e22'
                      limits=c('obs','05th','10th','50th','90th','95th','mean'))+
  scale_size_manual(values = c('obs'=0.6,'05th'=1,'10th'=1,
                               '50th'=1,'90th'=1,'95th'=1,'mean'=1)) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year", limit = c(as.Date("1997-09-01"),as.Date("2022-12-01"))) + 
  scale_y_continuous(breaks = seq(0,0.5,0.1)) +
  coord_cartesian(xlim = c(as.Date("1997-06-01"),as.Date("2022-12-01")),
                  ylim = c(-0.06,0.5),expand=F)+ 
  xlab("Time (yr)") +
  ylab(expression(paste("CHL  ","(mg·m"^"-3",")"))) +
  ggtitle(paste0("North Pacific Subarctic Gyres Province")) +
  guides(color=guide_legend(title = 'quantile levels',title.position = "top",title.hjust = .5, show.limits = T,
                            label.position = 'left', direction = 'verticle',ncol = 1, byrow = F),size=F,
         guide_colorbar(barwidth = unit(2, "cm"),barheight = unit(.5, "cm"))) +
  theme_bw()+
  theme(legend.position = 'right',plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"),
        legend.title = element_text(size = 10),plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.7, 'cm'),
        axis.text.x=element_text(angle=60, hjust=1))

ggsave("np_quants_new.png",width=7.5, height=3.44, dpi=300)


# NA
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/na_fit_ts.Rdata")
quantile_df <- as.data.frame(abind(na_ts,ts.05,ts.10,ts.50,
                                   ts.90,ts.95,ts.lm))
colnames(quantile_df) <- c('fitted')
Date <- seq(from = as.Date("1997/09/01",format="%Y/%m/%d"), 
            to = as.Date("2022/12/01",format="%Y/%m/%d"), 
            by = "1 month")
quantile_df$date <- rep(Date,time=7)
quantile_df$quantile <- rep(c('obs','05th','10th','50th','90th','95th','mean'),each=304)
range(quantile_df$fitted)
### plot
plot <- ggplot() +
  # geom_line(data = data_ts,aes(x=date,y=mon_deszn),alpha=0.8,position = "identity",col='#bdc3c7')+
  geom_line(data=quantile_df,aes(x=date,y=fitted,group=quantile,colour=as.factor(quantile),size=quantile))+ 
  scale_colour_manual(values=c('#bdc3c7','#2e86de','#c7ecee','#2ecc71','#fad390','#e74c3c','#6c5ce7'), # '#48dbfb',,'#e67e22'
                      limits=c('obs','05th','10th','50th','90th','95th','mean'))+
  scale_size_manual(values = c('obs'=0.6,'05th'=1,'10th'=1,
                               '50th'=1,'90th'=1,'95th'=1,'mean'=1)) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year", limit = c(as.Date("1997-09-01"),as.Date("2022-12-01"))) + 
  scale_y_continuous(breaks = seq(0,0.45,0.1)) +
  coord_cartesian(xlim = c(as.Date("1997-06-01"),as.Date("2022-12-01")),
                  ylim = c(-0.1,0.45),expand=F)+ 
  xlab("Time (yr)") +
  ylab(expression(paste("CHL  ","(mg·m"^"-3",")"))) +
  ggtitle(paste0("North Atlantic Drift Province")) +
  guides(color=guide_legend(title = 'quantile levels',title.position = "top",title.hjust = .5, show.limits = T,
                            label.position = 'left', direction = 'verticle',ncol = 1, byrow = F),size=F,
         guide_colorbar(barwidth = unit(2, "cm"),barheight = unit(.5, "cm"))) +
  theme_bw()+
  theme(legend.position = 'right',plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"),
        legend.title = element_text(size = 10),plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.7, 'cm'),
        axis.text.x=element_text(angle=60, hjust=1))

ggsave("na_quants_new.png",width=7.5, height=3.44, dpi=300)


# EP
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/ep_fit_ts.Rdata")
quantile_df <- as.data.frame(abind(ep_ts,ts.05,ts.10,ts.50,
                                   ts.90,ts.95,ts.lm))
colnames(quantile_df) <- c('fitted')
Date <- seq(from = as.Date("1997/09/01",format="%Y/%m/%d"), 
            to = as.Date("2022/12/01",format="%Y/%m/%d"), 
            by = "1 month")
quantile_df$date <- rep(Date,time=7)
quantile_df$quantile <- rep(c('obs','05th','10th','50th','90th','95th','mean'),each=304)
range(quantile_df$fitted)
### plot
plot <- ggplot() +
  # geom_line(data = data_ts,aes(x=date,y=mon_deszn),alpha=0.8,position = "identity",col='#bdc3c7')+
  geom_line(data=quantile_df,aes(x=date,y=fitted,group=quantile,colour=as.factor(quantile),size=quantile))+ 
  scale_colour_manual(values=c('#bdc3c7','#2e86de','#c7ecee','#2ecc71','#fad390','#e74c3c','#6c5ce7'), # '#48dbfb',,'#e67e22'
                      limits=c('obs','05th','10th','50th','90th','95th','mean'))+
  scale_size_manual(values = c('obs'=0.6,'05th'=1,'10th'=1,
                               '50th'=1,'90th'=1,'95th'=1,'mean'=1)) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year", limit = c(as.Date("1997-09-01"),as.Date("2022-12-01"))) + 
  scale_y_continuous(breaks = seq(0,0.1,0.02)) +
  coord_cartesian(xlim = c(as.Date("1997-06-01"),as.Date("2022-12-01")),
                  ylim = c(0.01,0.1),expand=F)+ 
  xlab("Time (yr)") +
  ylab(expression(paste("CHL  ","(mg·m"^"-3",")"))) +
  ggtitle(paste0("Pacific Equatorial Province")) +
  guides(color=guide_legend(title = 'quantile levels',title.position = "top",title.hjust = .5, show.limits = T,
                            label.position = 'left', direction = 'verticle',ncol = 1, byrow = F),size=F,
         guide_colorbar(barwidth = unit(2, "cm"),barheight = unit(.5, "cm"))) +
  theme_bw()+
  theme(legend.position = 'right',plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"),
        legend.title = element_text(size = 10),plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.7, 'cm'),
        axis.text.x=element_text(angle=60, hjust=1))

ggsave("ep_quants_new.png",width=7.5, height=3.44, dpi=300)

# SO
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/so_fit_ts.Rdata")
quantile_df <- as.data.frame(abind(so_ts,ts.05,ts.10,ts.50,
                                   ts.90,ts.95,ts.lm))
colnames(quantile_df) <- c('fitted')
Date <- seq(from = as.Date("1997/09/01",format="%Y/%m/%d"), 
            to = as.Date("2022/12/01",format="%Y/%m/%d"), 
            by = "1 month")
quantile_df$date <- rep(Date,time=7)
quantile_df$quantile <- rep(c('obs','05th','10th','50th','90th','95th','mean'),each=304)
range(quantile_df$fitted)
### plot
plot <- ggplot() +
  # geom_line(data = data_ts,aes(x=date,y=mon_deszn),alpha=0.8,position = "identity",col='#bdc3c7')+
  geom_line(data=quantile_df,aes(x=date,y=fitted,group=quantile,colour=as.factor(quantile),size=quantile))+ 
  scale_colour_manual(values=c('#bdc3c7','#2e86de','#c7ecee','#2ecc71','#fad390','#e74c3c','#6c5ce7'), # '#48dbfb',,'#e67e22'
                      limits=c('obs','05th','10th','50th','90th','95th','mean'))+
  scale_size_manual(values = c('obs'=0.6,'05th'=1,'10th'=1,
                               '50th'=1,'90th'=1,'95th'=1,'mean'=1)) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year", limit = c(as.Date("1997-09-01"),as.Date("2022-12-01"))) + 
  scale_y_continuous(breaks = seq(0,0.12,0.02)) +
  coord_cartesian(xlim = c(as.Date("1997-06-01"),as.Date("2022-12-01")),
                  ylim = c(0,0.12),expand=F)+ 
  xlab("Time (yr)") +
  ylab(expression(paste("CHL  ","(mg·m"^"-3",")"))) +
  ggtitle(paste0("Subantarctic Province")) +
  guides(color=guide_legend(title = 'quantile levels',title.position = "top",title.hjust = .5, show.limits = T,
                            label.position = 'left', direction = 'verticle',ncol = 1, byrow = F),size=F,
         guide_colorbar(barwidth = unit(2, "cm"),barheight = unit(.5, "cm"))) +
  theme_bw()+
  theme(legend.position = 'right',plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"),
        legend.title = element_text(size = 10),plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.7, 'cm'),
        axis.text.x=element_text(angle=60, hjust=1))

ggsave("so_quants_new.png",width=7.5, height=3.44, dpi=300)

# NPO
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/npo_fit_ts.Rdata")
quantile_df <- as.data.frame(abind(npo_ts,ts.05,ts.10,ts.50,
                                   ts.90,ts.95,ts.lm))
colnames(quantile_df) <- c('fitted')
Date <- seq(from = as.Date("1997/09/01",format="%Y/%m/%d"), 
            to = as.Date("2022/12/01",format="%Y/%m/%d"), 
            by = "1 month")
quantile_df$date <- rep(Date,time=7)
quantile_df$quantile <- rep(c('obs','05th','10th','50th','90th','95th','mean'),each=304)
range(quantile_df$fitted)
### plot
plot <- ggplot() +
  # geom_line(data = data_ts,aes(x=date,y=mon_deszn),alpha=0.8,position = "identity",col='#bdc3c7')+
  geom_line(data=quantile_df,aes(x=date,y=fitted,group=quantile,colour=as.factor(quantile),size=quantile))+ 
  scale_colour_manual(values=c('#bdc3c7','#2e86de','#c7ecee','#2ecc71','#fad390','#e74c3c','#6c5ce7'), # '#48dbfb',,'#e67e22'
                      limits=c('obs','05th','10th','50th','90th','95th','mean'))+
  scale_size_manual(values = c('obs'=0.6,'05th'=1,'10th'=1,
                               '50th'=1,'90th'=1,'95th'=1,'mean'=1)) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year", limit = c(as.Date("1997-09-01"),as.Date("2022-12-01"))) + 
  scale_y_continuous(breaks = seq(0,0.03,0.01)) +
  coord_cartesian(xlim = c(as.Date("1997-06-01"),as.Date("2022-12-01")),
                  ylim = c(-0.000,0.032),expand=F)+ 
  xlab("Time (yr)") +
  ylab(expression(paste("CHL  ","(mg·m"^"-3",")"))) +
  ggtitle(paste0("North Pacific Subtropical Gyre Province")) +
  guides(color=guide_legend(title = 'quantile levels',title.position = "top",title.hjust = .5, show.limits = T,
                            label.position = 'left', direction = 'verticle',ncol = 1, byrow = F),size=F,
         guide_colorbar(barwidth = unit(2, "cm"),barheight = unit(.5, "cm"))) +
  theme_bw()+
  theme(legend.position = 'right',plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"),
        legend.title = element_text(size = 10),plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.7, 'cm'),
        axis.text.x=element_text(angle=60, hjust=1))

ggsave("npo_quants_new.png",width=7.5, height=3.44, dpi=300)

# NAO
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/nao_fit_ts.Rdata")
quantile_df <- as.data.frame(abind(nao_ts,ts.05,ts.10,ts.50,
                                   ts.90,ts.95,ts.lm))
colnames(quantile_df) <- c('fitted')
Date <- seq(from = as.Date("1997/09/01",format="%Y/%m/%d"), 
            to = as.Date("2022/12/01",format="%Y/%m/%d"), 
            by = "1 month")
quantile_df$date <- rep(Date,time=7)
quantile_df$quantile <- rep(c('obs','05th','10th','50th','90th','95th','mean'),each=304)
range(quantile_df$fitted)
### plot
plot <- ggplot() +
  # geom_line(data = data_ts,aes(x=date,y=mon_deszn),alpha=0.8,position = "identity",col='#bdc3c7')+
  geom_line(data=quantile_df,aes(x=date,y=fitted,group=quantile,colour=as.factor(quantile),size=quantile))+ 
  scale_colour_manual(values=c('#bdc3c7','#2e86de','#c7ecee','#2ecc71','#fad390','#e74c3c','#6c5ce7'), # '#48dbfb',,'#e67e22'
                      limits=c('obs','05th','10th','50th','90th','95th','mean'))+
  scale_size_manual(values = c('obs'=0.6,'05th'=1,'10th'=1,
                               '50th'=1,'90th'=1,'95th'=1,'mean'=1)) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year", limit = c(as.Date("1997-09-01"),as.Date("2022-12-01"))) +  
  scale_y_continuous(breaks = seq(-0.01,0.06,0.01)) +
  coord_cartesian(xlim = c(as.Date("1997-06-01"),as.Date("2022-12-01")),
                  ylim = c(-0.01,0.07),expand=F)+ 
  xlab("Time (yr)") +
  ylab(expression(paste("CHL  ","(mg·m"^"-3",")"))) +
  ggtitle(paste0("North Atlantic Subtropical Gyre Province")) +
  guides(color=guide_legend(title = 'quantile levels',title.position = "top",title.hjust = .5, show.limits = T,
                            label.position = 'left', direction = 'verticle',ncol = 1, byrow = F),size=F,
         guide_colorbar(barwidth = unit(2, "cm"),barheight = unit(.5, "cm"))) +
  theme_bw()+
  theme(legend.position = 'right',plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"),
        legend.title = element_text(size = 10),plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.7, 'cm'),
        axis.text.x=element_text(angle=60, hjust=1))

ggsave("nao_quants_new.png",width=7.5, height=3.44, dpi=300)


########## Part III. Timeseries with multiple quantile levels facet ##########


