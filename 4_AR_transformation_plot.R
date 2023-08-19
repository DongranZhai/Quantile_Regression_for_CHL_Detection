# Note: Plotting
#  Part I. Test statistic index projection
#          Cochrane_Orcutt
#.         Hildreth-Lu

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

############### ############### ############### ############### ############### 
###  Part I. Test statistic index projection  ####
setwd("/Volumes/Doris/OCCCI/quantile_regression/global/Cochrane_Orcutt")
setwd("/Volumes/Doris/OCCCI/quantile_regression/global/Hildreth_Lu")
setwd("/Volumes/Doris/OCCCI/quantile_regression/log_trans")
###### 1% ######
# significant trend II (grey points shadowed) <GOOD!!!!>
rm(list=ls())
# load("/Volumes/Doris/OCCCI/quantile_regression/global/Cochrane_Orcutt/ar_cochrane_01.Rdata")
# load('/Volumes/Doris/OCCCI/quantile_regression/global/Hildreth_Lu/ar_hildreth_01.Rdata')
load('/Volumes/Doris/OCCCI/quantile_regression/log_trans/ar_cochrane_01.Rdata')
wmap<-readOGR(dsn="/Volumes/Doris/Projects/spt_wc/wmap", layer="ne_110m_land")
wmap_df <- fortify(wmap)
wmap_df <-wmap_df[-(which(wmap_df$group==112.2)),]#remove caspian sea
com <- rev(c("firebrick4","firebrick2","white","dodgerblue2","dodgerblue4"))#color palette

trendvalue <- tr_trans*1000 # magnitude (100 or 1000) depend on actual values but need to be consistent
trend_raster <- matrix2raster(trendvalue, x = lon, y = lat, layer = 2)
trend_df <- as.data.frame(trend_raster,xy = T)
colnames(trend_df) <-c('lon','lat','trend')
rm(list=c('trend_raster'))

ptemp <- pValue[,ncol(pValue):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$pvalue <- pValuetemp
trend_df$sigtrend <- NA

ind <- which(trend_df$pvalue <= 0.05)
trend_df$sigtrend[ind] <- "Yes"
ind <- which(is.na(trend_df$pvalue))
trend_df$sigtrend[ind] <- "No" # NA

trend_df$sigtrend <- as.factor(trend_df$sigtrend)

plot <- ggplot() +
  geom_raster(data = trend_df, aes(x = lon,y = lat,fill = trend)) +
  scale_fill_gradientn(colors = com,limits =c(-2,2),na.value = 'black',oob=squish) +
  geom_point(data =trend_df, aes(x = lon,y = lat,col=sigtrend),size=0.00000001,show.legend = F) +
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) + #
  geom_polygon(data=wmap_df, aes(x=long, y=lat, group = group),colour="black", fill="gainsboro",size=0.25 )+
  coord_equal() + 
  coord_cartesian(ylim = c(-75, 75),xlim = c(-180,180),expand=F)+
  scale_x_continuous(breaks = seq(-180,180,60)) +
  scale_y_continuous(breaks = seq(-90,90,30)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-3"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))

# ggsave(filename = paste0("01_cochrane.png"), width=8.27, height=3.44, dpi=300)
# ggsave(filename = paste0("01_hildreth.png"), width=8.27, height=3.44, dpi=300)
# ggsave(filename = paste0("01_log_cochrane.png"), width=8.27, height=3.44, dpi=300)
ggsave(filename = paste0("01_cochrane.png"), width=8.27, height=3.44, dpi=300)
###### 5% ######
# significant trend II (grey points shadowed) <GOOD!!!!>
rm(list=ls())
# load("/Volumes/Doris/OCCCI/quantile_regression/global/Cochrane_Orcutt/ar_cochrane_05.Rdata")
load('/Volumes/Doris/OCCCI/quantile_regression/log_trans/ar_cochrane_05.Rdata')
wmap<-readOGR(dsn="/Volumes/Doris/Projects/spt_wc/wmap", layer="ne_110m_land")
wmap_df <- fortify(wmap)
wmap_df <-wmap_df[-(which(wmap_df$group==112.2)),]#remove caspian sea
com <- rev(c("firebrick4","firebrick2","white","dodgerblue2","dodgerblue4"))#color palette

trendvalue <- tr_trans*1000 # magnitude (100 or 1000) depend on actual values but need to be consistent
trend_raster <- matrix2raster(trendvalue, x = lon, y = lat, layer = 2)
trend_df <- as.data.frame(trend_raster,xy = T)
colnames(trend_df) <-c('lon','lat','trend')
rm(list=c('trend_raster'))

ptemp <- pValue[,ncol(pValue):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$pvalue <- pValuetemp
trend_df$sigtrend <- NA

ind <- which(trend_df$pvalue <= 0.05)
trend_df$sigtrend[ind] <- "Yes"
ind <- which(is.na(trend_df$pvalue))
trend_df$sigtrend[ind] <- "No" # NA

trend_df$sigtrend <- as.factor(trend_df$sigtrend)

plot <- ggplot() +
  geom_raster(data = trend_df, aes(x = lon,y = lat,fill = trend)) +
  scale_fill_gradientn(colors = com,limits =c(-2,2),na.value = 'black',oob=squish) +
  geom_point(data =trend_df, aes(x = lon,y = lat,col=sigtrend),size=0.00000001,show.legend = F) +
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) + #
  geom_polygon(data=wmap_df, aes(x=long, y=lat, group = group),colour="black", fill="gainsboro",size=0.25 )+
  coord_equal() + 
  coord_cartesian(ylim = c(-75, 75),xlim = c(-180,180),expand=F)+
  scale_x_continuous(breaks = seq(-180,180,60)) +
  scale_y_continuous(breaks = seq(-90,90,30)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-3"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))

# ggsave(filename = paste0("05_cochrane.png"), width=8.27, height=3.44, dpi=300)
ggsave(filename = paste0("05_log_cochrane.png"), width=8.27, height=3.44, dpi=300)

###### 10% ######
# significant trend II (grey points shadowed) <GOOD!!!!>
rm(list=ls())
# load("/Volumes/Doris/OCCCI/quantile_regression/global/Cochrane_Orcutt/ar_cochrane_10.Rdata")
load('/Volumes/Doris/OCCCI/quantile_regression/log_trans/ar_cochrane_10.Rdata')
wmap<-readOGR(dsn="/Volumes/Doris/Projects/spt_wc/wmap", layer="ne_110m_land")
wmap_df <- fortify(wmap)
wmap_df <-wmap_df[-(which(wmap_df$group==112.2)),]#remove caspian sea
com <- rev(c("firebrick4","firebrick2","white","dodgerblue2","dodgerblue4"))#color palette

trendvalue <- tr_trans*1000 # magnitude (100 or 1000) depend on actual values but need to be consistent
trend_raster <- matrix2raster(trendvalue, x = lon, y = lat, layer = 2)
trend_df <- as.data.frame(trend_raster,xy = T)
colnames(trend_df) <-c('lon','lat','trend')
rm(list=c('trend_raster'))

ptemp <- pValue[,ncol(pValue):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$pvalue <- pValuetemp
trend_df$sigtrend <- NA

ind <- which(trend_df$pvalue <= 0.05)
trend_df$sigtrend[ind] <- "Yes"
ind <- which(is.na(trend_df$pvalue))
trend_df$sigtrend[ind] <- "No" # NA

trend_df$sigtrend <- as.factor(trend_df$sigtrend)

plot <- ggplot() +
  geom_raster(data = trend_df, aes(x = lon,y = lat,fill = trend)) +
  scale_fill_gradientn(colors = com,limits =c(-2,2),na.value = 'black',oob=squish) +
  geom_point(data =trend_df, aes(x = lon,y = lat,col=sigtrend),size=0.00000001,show.legend = F) +
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) + #
  geom_polygon(data=wmap_df, aes(x=long, y=lat, group = group),colour="black", fill="gainsboro",size=0.25 )+
  coord_equal() + 
  coord_cartesian(ylim = c(-75, 75),xlim = c(-180,180),expand=F)+
  scale_x_continuous(breaks = seq(-180,180,60)) +
  scale_y_continuous(breaks = seq(-90,90,30)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-3"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))

# ggsave(filename = paste0("10_cochrane.png"), width=8.27, height=3.44, dpi=300)
ggsave(filename = paste0("10_log_cochrane.png"), width=8.27, height=3.44, dpi=300)

###### 50% ######
# significant trend II (grey points shadowed) <GOOD!!!!>
rm(list=ls())
# load("/Volumes/Doris/OCCCI/quantile_regression/global/Cochrane_Orcutt/ar_cochrane_50.Rdata")
load('/Volumes/Doris/OCCCI/quantile_regression/log_trans/ar_cochrane_50.Rdata')
wmap<-readOGR(dsn="/Volumes/Doris/Projects/spt_wc/wmap", layer="ne_110m_land")
wmap_df <- fortify(wmap)
wmap_df <-wmap_df[-(which(wmap_df$group==112.2)),]#remove caspian sea
com <- rev(c("firebrick4","firebrick2","white","dodgerblue2","dodgerblue4"))#color palette

trendvalue <- tr_trans*1000 # magnitude (100 or 1000) depend on actual values but need to be consistent
trend_raster <- matrix2raster(trendvalue, x = lon, y = lat, layer = 2)
trend_df <- as.data.frame(trend_raster,xy = T)
colnames(trend_df) <-c('lon','lat','trend')
rm(list=c('trend_raster'))

ptemp <- pValue[,ncol(pValue):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$pvalue <- pValuetemp
trend_df$sigtrend <- NA

ind <- which(trend_df$pvalue <= 0.05)
trend_df$sigtrend[ind] <- "Yes"
ind <- which(is.na(trend_df$pvalue))
trend_df$sigtrend[ind] <- "No" # NA

trend_df$sigtrend <- as.factor(trend_df$sigtrend)

plot <- ggplot() +
  geom_raster(data = trend_df, aes(x = lon,y = lat,fill = trend)) +
  scale_fill_gradientn(colors = com,limits =c(-2,2),na.value = 'black',oob=squish) +
  geom_point(data =trend_df, aes(x = lon,y = lat,col=sigtrend),size=0.00000001,show.legend = F) +
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) + #
  geom_polygon(data=wmap_df, aes(x=long, y=lat, group = group),colour="black", fill="gainsboro",size=0.25 )+
  coord_equal() + 
  coord_cartesian(ylim = c(-75, 75),xlim = c(-180,180),expand=F)+
  scale_x_continuous(breaks = seq(-180,180,60)) +
  scale_y_continuous(breaks = seq(-90,90,30)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-3"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))

# ggsave(filename = paste0("50_cochrane.png"), width=8.27, height=3.44, dpi=300)
ggsave(filename = paste0("50_log_cochrane.png"), width=8.27, height=3.44, dpi=300)

###### 90% ######
# significant trend II (grey points shadowed) <GOOD!!!!>
rm(list=ls())
# load("/Volumes/Doris/OCCCI/quantile_regression/global/Cochrane_Orcutt/ar_cochrane_90.Rdata")
load('/Volumes/Doris/OCCCI/quantile_regression/log_trans/ar_cochrane_90.Rdata')
wmap<-readOGR(dsn="/Volumes/Doris/Projects/spt_wc/wmap", layer="ne_110m_land")
wmap_df <- fortify(wmap)
wmap_df <-wmap_df[-(which(wmap_df$group==112.2)),]#remove caspian sea
com <- rev(c("firebrick4","firebrick2","white","dodgerblue2","dodgerblue4"))#color palette

trendvalue <- tr_trans*1000 # magnitude (100 or 1000) depend on actual values but need to be consistent
trend_raster <- matrix2raster(trendvalue, x = lon, y = lat, layer = 2)
trend_df <- as.data.frame(trend_raster,xy = T)
colnames(trend_df) <-c('lon','lat','trend')
rm(list=c('trend_raster'))

ptemp <- pValue[,ncol(pValue):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$pvalue <- pValuetemp
trend_df$sigtrend <- NA

ind <- which(trend_df$pvalue <= 0.05)
trend_df$sigtrend[ind] <- "Yes"
ind <- which(is.na(trend_df$pvalue))
trend_df$sigtrend[ind] <- "No" # NA

trend_df$sigtrend <- as.factor(trend_df$sigtrend)

plot <- ggplot() +
  geom_raster(data = trend_df, aes(x = lon,y = lat,fill = trend)) +
  scale_fill_gradientn(colors = com,limits =c(-2,2),na.value = 'black',oob=squish) +
  geom_point(data =trend_df, aes(x = lon,y = lat,col=sigtrend),size=0.00000001,show.legend = F) +
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) + #
  geom_polygon(data=wmap_df, aes(x=long, y=lat, group = group),colour="black", fill="gainsboro",size=0.25 )+
  coord_equal() + 
  coord_cartesian(ylim = c(-75, 75),xlim = c(-180,180),expand=F)+
  scale_x_continuous(breaks = seq(-180,180,60)) +
  scale_y_continuous(breaks = seq(-90,90,30)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-3"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))

# ggsave(filename = paste0("90_cochrane.png"), width=8.27, height=3.44, dpi=300)
ggsave(filename = paste0("90_log_cochrane.png"), width=8.27, height=3.44, dpi=300)

###### 95% ######
# significant trend II (grey points shadowed) <GOOD!!!!>
rm(list=ls())
# load("/Volumes/Doris/OCCCI/quantile_regression/global/Cochrane_Orcutt/ar_cochrane_95.Rdata")
load('/Volumes/Doris/OCCCI/quantile_regression/log_trans/ar_cochrane_95.Rdata')
wmap<-readOGR(dsn="/Volumes/Doris/Projects/spt_wc/wmap", layer="ne_110m_land")
wmap_df <- fortify(wmap)
wmap_df <-wmap_df[-(which(wmap_df$group==112.2)),]#remove caspian sea
com <- rev(c("firebrick4","firebrick2","white","dodgerblue2","dodgerblue4"))#color palette

trendvalue <- tr_trans*1000 # magnitude (100 or 1000) depend on actual values but need to be consistent
trend_raster <- matrix2raster(trendvalue, x = lon, y = lat, layer = 2)
trend_df <- as.data.frame(trend_raster,xy = T)
colnames(trend_df) <-c('lon','lat','trend')
rm(list=c('trend_raster'))

ptemp <- pValue[,ncol(pValue):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$pvalue <- pValuetemp
trend_df$sigtrend <- NA

ind <- which(trend_df$pvalue <= 0.05)
trend_df$sigtrend[ind] <- "Yes"
ind <- which(is.na(trend_df$pvalue))
trend_df$sigtrend[ind] <- "No" # NA

trend_df$sigtrend <- as.factor(trend_df$sigtrend)

plot <- ggplot() +
  geom_raster(data = trend_df, aes(x = lon,y = lat,fill = trend)) +
  scale_fill_gradientn(colors = com,limits =c(-2,2),na.value = 'black',oob=squish) +
  geom_point(data =trend_df, aes(x = lon,y = lat,col=sigtrend),size=0.00000001,show.legend = F) +
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) + #
  geom_polygon(data=wmap_df, aes(x=long, y=lat, group = group),colour="black", fill="gainsboro",size=0.25 )+
  coord_equal() + 
  coord_cartesian(ylim = c(-75, 75),xlim = c(-180,180),expand=F)+
  scale_x_continuous(breaks = seq(-180,180,60)) +
  scale_y_continuous(breaks = seq(-90,90,30)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-3"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))

# ggsave(filename = paste0("95_cochrane.png"), width=8.27, height=3.44, dpi=300)
ggsave(filename = paste0("95_log_cochrane.png"), width=8.27, height=3.44, dpi=300)

###### 99% ######
# significant trend II (grey points shadowed) <GOOD!!!!>
rm(list=ls())
# load("/Volumes/Doris/OCCCI/quantile_regression/global/Cochrane_Orcutt/ar_cochrane_99.Rdata")
load('/Volumes/Doris/OCCCI/quantile_regression/log_trans/ar_cochrane_99.Rdata')
wmap<-readOGR(dsn="/Volumes/Doris/Projects/spt_wc/wmap", layer="ne_110m_land")
wmap_df <- fortify(wmap)
wmap_df <-wmap_df[-(which(wmap_df$group==112.2)),]#remove caspian sea
com <- rev(c("firebrick4","firebrick2","white","dodgerblue2","dodgerblue4"))#color palette

trendvalue <- tr_trans*1000 # magnitude (100 or 1000) depend on actual values but need to be consistent
trend_raster <- matrix2raster(trendvalue, x = lon, y = lat, layer = 2)
trend_df <- as.data.frame(trend_raster,xy = T)
colnames(trend_df) <-c('lon','lat','trend')
rm(list=c('trend_raster'))

ptemp <- pValue[,ncol(pValue):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$pvalue <- pValuetemp
trend_df$sigtrend <- NA

ind <- which(trend_df$pvalue <= 0.05)
trend_df$sigtrend[ind] <- "Yes"
ind <- which(is.na(trend_df$pvalue))
trend_df$sigtrend[ind] <- "No" # NA

trend_df$sigtrend <- as.factor(trend_df$sigtrend)

plot <- ggplot() +
  geom_raster(data = trend_df, aes(x = lon,y = lat,fill = trend)) +
  scale_fill_gradientn(colors = com,limits =c(-2,2),na.value = 'black',oob=squish) +
  geom_point(data =trend_df, aes(x = lon,y = lat,col=sigtrend),size=0.00000001,show.legend = F) +
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) + #
  geom_polygon(data=wmap_df, aes(x=long, y=lat, group = group),colour="black", fill="gainsboro",size=0.25 )+
  coord_equal() + 
  coord_cartesian(ylim = c(-75, 75),xlim = c(-180,180),expand=F)+
  scale_x_continuous(breaks = seq(-180,180,60)) +
  scale_y_continuous(breaks = seq(-90,90,30)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-3"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))

# ggsave(filename = paste0("99_cochrane.png"), width=8.27, height=3.44, dpi=300)
ggsave(filename = paste0("99_log_cochrane.png"), width=8.27, height=3.44, dpi=300)

###### simple linear model ######
# significant trend II (grey points shadowed) <GOOD!!!!>
rm(list=ls())
load("/Volumes/Doris/OCCCI/quantile_regression/log_trans/ar_cochrane_lm.Rdata")
# load("/Volumes/Doris/OCCCI/quantile_regression/global/Hildreth_Lu/ar_hildreth_lm.Rdata")
wmap<-readOGR(dsn="/Volumes/Doris/Projects/spt_wc/wmap", layer="ne_110m_land")
wmap_df <- fortify(wmap)
wmap_df <-wmap_df[-(which(wmap_df$group==112.2)),]#remove caspian sea
com <- rev(c("firebrick4","firebrick2","white","dodgerblue2","dodgerblue4"))#color palette

trendvalue <- tr_trans*1000 # magnitude (100 or 1000) depend on actual values but need to be consistent
trend_raster <- matrix2raster(trendvalue, x = lon, y = lat, layer = 2)
trend_df <- as.data.frame(trend_raster,xy = T)
colnames(trend_df) <-c('lon','lat','trend')
rm(list=c('trend_raster'))

ptemp <- pValue[,ncol(pValue):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$pvalue <- pValuetemp
trend_df$sigtrend <- NA

ind <- which(trend_df$pvalue <= 0.05)
trend_df$sigtrend[ind] <- "Yes"
ind <- which(is.na(trend_df$pvalue))
trend_df$sigtrend[ind] <- "No" # NA

trend_df$sigtrend <- as.factor(trend_df$sigtrend)

plot <- ggplot() +
  geom_raster(data = trend_df, aes(x = lon,y = lat,fill = trend)) +
  scale_fill_gradientn(colors = com,limits =c(-2,2),na.value = 'black',oob=squish) +
  geom_point(data =trend_df, aes(x = lon,y = lat,col=sigtrend),size=0.00000001,show.legend = F) +
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) + #
  geom_polygon(data=wmap_df, aes(x=long, y=lat, group = group),colour="black", fill="gainsboro",size=0.25 )+
  coord_equal() + 
  coord_cartesian(ylim = c(-75, 75),xlim = c(-180,180),expand=F)+
  scale_x_continuous(breaks = seq(-180,180,60)) +
  scale_y_continuous(breaks = seq(-90,90,30)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-3"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))

ggsave(filename = paste0("lm_log_cochrane.png"), width=8.27, height=3.44, dpi=300)
# ggsave(filename = paste0("lm_hildreth.png"), width=8.27, height=3.44, dpi=300)

