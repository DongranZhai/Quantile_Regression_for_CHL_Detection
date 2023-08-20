# Note: seasonality plot
library(quantreg)
library(ggplot2)
library(raster)
library(oceanmap)
library(scales)
library(rgdal)
library(mapproj)
library(munsell)

rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression/seasonality")
########## Part I: Winter NP ##########
load("/Volumes/Doris/OCCCI/quantile_regression/winter_regions_quant.Rdata")
mapp_df <- fortify(map("world2",plot=FALSE,fill=TRUE))
com <- rev(c("firebrick4","firebrick2","white","dodgerblue2","dodgerblue4"))#color palette

##### 50
trend_df <- np_centre_df_50 # df change
trend_df$trend <- trend_df$trend*100
trend_df$sigtrend <- NA
ind <- which(trend_df$pvalue <= 0.05)
trend_df$sigtrend[ind] <- "Yes"
ind <- which(is.na(trend_df$pvalue))
trend_df$sigtrend[ind] <- "No" # NA
trend_df$sigtrend <- as.factor(trend_df$sigtrend)

plot <- ggplot() +
  geom_raster(data = trend_df, aes(x = lon,y = lat,fill = trend)) +
  scale_fill_gradientn(colors = com,limits =c(-2,2),na.value = 'black',oob=squish) +
  geom_point(data =trend_df, aes(x = lon,y = lat,col=sigtrend),size=0.8,show.legend = F) +
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) +
  geom_polygon( data=mapp_df, aes(x=long, y=lat, group = group),colour="gainsboro", fill="gainsboro",size=0.25 ) +
  coord_cartesian(ylim = c(30,60),xlim=c(120,240),expand=F) +
  # scale_x_continuous(breaks = seq(-180,180,60)) + # npsg
  scale_x_continuous(breaks = seq(120,240,30)) + # centre
  scale_y_continuous(breaks = seq(30,60,10)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-2"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))
savename <- paste0("np_winter_50.png")
ggsave(savename, width=8.27, height=3.44, dpi=300)


##### 95th
trend_df <- np_centre_df_95
trend_df$trend <- trend_df$trend*100
trend_df$sigtrend <- NA
ind <- which(trend_df$pvalue <= 0.05)
trend_df$sigtrend[ind] <- "Yes"
ind <- which(is.na(trend_df$pvalue))
trend_df$sigtrend[ind] <- "No" # NA
trend_df$sigtrend <- as.factor(trend_df$sigtrend)

plot <- ggplot() +
  geom_raster(data = trend_df, aes(x = lon,y = lat,fill = trend)) +
  scale_fill_gradientn(colors = com,limits =c(-2,2),na.value = 'black',oob=squish) +
  geom_point(data =trend_df, aes(x = lon,y = lat,col=sigtrend),size=0.8,show.legend = F) +
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) +
  geom_polygon( data=mapp_df, aes(x=long, y=lat, group = group),colour="gainsboro", fill="gainsboro",size=0.25 ) +
  coord_cartesian(ylim = c(30,60),xlim=c(120,240),expand=F) +
  # scale_x_continuous(breaks = seq(-180,180,60)) + # npsg
  scale_x_continuous(breaks = seq(120,240,30)) + # centre
  scale_y_continuous(breaks = seq(30,60,10)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-2"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))

savename <- paste0("np_winter_95.png")
ggsave(savename, width=8.27, height=3.44, dpi=300)

########## Part II: Winter NA ##########
rm(list=ls())
load("/Volumes/Doris/OCCCI/quantile_regression/winter_regions_quant.Rdata")
# check plot
mapp_df <- fortify(map("world2",plot=FALSE,fill=TRUE))
com <- rev(c("firebrick4","firebrick2","white","dodgerblue2","dodgerblue4"))#color palette

##### 50
trend_df <- na_centre_df_50 # df change
trend_df$trend <- trend_df$trend*100
trend_df$sigtrend <- NA
ind <- which(trend_df$pvalue <= 0.05)
trend_df$sigtrend[ind] <- "Yes"
ind <- which(is.na(trend_df$pvalue))
trend_df$sigtrend[ind] <- "No" # NA
trend_df$sigtrend <- as.factor(trend_df$sigtrend)

plot <- ggplot() +
  geom_raster(data = trend_df, aes(x = lon,y = lat,fill = trend)) +
  scale_fill_gradientn(colors = com,limits =c(-2,2),na.value = 'black',oob=squish) +
  geom_point(data =trend_df, aes(x = lon,y = lat,col=sigtrend),size=0.8,show.legend = F) +
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) +
  geom_polygon( data=mapp_df, aes(x=long, y=lat, group = group),colour="gainsboro", fill="gainsboro",size=0.25 ) +
  coord_cartesian(ylim = c(30,60),xlim=c(300,360),expand=F) +
  # scale_x_continuous(breaks = seq(-180,180,60)) + # npsg
  scale_x_continuous(breaks = seq(300,360,30)) + # centre
  scale_y_continuous(breaks = seq(30,60,10)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-3"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))
savename <- paste0("na_winter_50.png")
ggsave(savename, width=8.27, height=3.44, dpi=300)


##### 95th
trend_df <- na_centre_df_95
trend_df$trend <- trend_df$trend*100
trend_df$sigtrend <- NA
ind <- which(trend_df$pvalue <= 0.05)
trend_df$sigtrend[ind] <- "Yes"
ind <- which(is.na(trend_df$pvalue))
trend_df$sigtrend[ind] <- "No" # NA
trend_df$sigtrend <- as.factor(trend_df$sigtrend)

plot <- ggplot() +
  geom_raster(data = trend_df, aes(x = lon,y = lat,fill = trend)) +
  scale_fill_gradientn(colors = com,limits =c(-2,2),na.value = 'black',oob=squish) +
  geom_point(data =trend_df, aes(x = lon,y = lat,col=sigtrend),size=0.8,show.legend = F) +
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) +
  geom_polygon( data=mapp_df, aes(x=long, y=lat, group = group),colour="gainsboro", fill="gainsboro",size=0.25 ) +
  coord_cartesian(ylim = c(30,60),xlim=c(300,360),expand=F) +
  # scale_x_continuous(breaks = seq(-180,180,60)) + # npsg
  scale_x_continuous(breaks = seq(300,360,30)) + # centre
  scale_y_continuous(breaks = seq(30,60,10)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-3"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))

savename <- paste0("na_winter_95.png")
ggsave(savename, width=8.27, height=3.44, dpi=300)

########## Part III: Winter EP ##########
rm(list=ls())
load("/Volumes/Doris/OCCCI/quantile_regression/winter_regions_quant.Rdata")
# check plot
mapp_df <- fortify(map("world2",plot=FALSE,fill=TRUE))
com <- rev(c("firebrick4","firebrick2","white","dodgerblue2","dodgerblue4"))#color palette

##### 50
trend_df <- ep_centre_df_50 # df change
trend_df$trend <- trend_df$trend*100
trend_df$sigtrend <- NA
ind <- which(trend_df$pvalue <= 0.05)
trend_df$sigtrend[ind] <- "Yes"
ind <- which(is.na(trend_df$pvalue))
trend_df$sigtrend[ind] <- "No" # NA
trend_df$sigtrend <- as.factor(trend_df$sigtrend)

plot <- ggplot() +
  geom_raster(data = trend_df, aes(x = lon,y = lat,fill = trend)) +
  scale_fill_gradientn(colors = com,limits =c(-2,2),na.value = 'black',oob=squish) +
  geom_point(data =trend_df, aes(x = lon,y = lat,col=sigtrend),size=0.8,show.legend = F) +
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) +
  geom_polygon( data=mapp_df, aes(x=long, y=lat, group = group),colour="gainsboro", fill="gainsboro",size=0.25 ) +
  coord_cartesian(ylim = c(-30,30),xlim=c(180,300),expand=F) +
  # scale_x_continuous(breaks = seq(-180,180,60)) + # npsg
  scale_x_continuous(breaks = seq(180,300,30)) + # centre
  scale_y_continuous(breaks = seq(-30,30,10)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-3"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))
savename <- paste0("ep_winter_50.png")
ggsave(savename, width=8.27, height=3.44, dpi=300)


##### 95th
trend_df <- ep_centre_df_95
trend_df$trend <- trend_df$trend*100
trend_df$sigtrend <- NA
ind <- which(trend_df$pvalue <= 0.05)
trend_df$sigtrend[ind] <- "Yes"
ind <- which(is.na(trend_df$pvalue))
trend_df$sigtrend[ind] <- "No" # NA
trend_df$sigtrend <- as.factor(trend_df$sigtrend)

plot <- ggplot() +
  geom_raster(data = trend_df, aes(x = lon,y = lat,fill = trend)) +
  scale_fill_gradientn(colors = com,limits =c(-2,2),na.value = 'black',oob=squish) +
  geom_point(data =trend_df, aes(x = lon,y = lat,col=sigtrend),size=0.8,show.legend = F) +
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) +
  geom_polygon( data=mapp_df, aes(x=long, y=lat, group = group),colour="gainsboro", fill="gainsboro",size=0.25 ) +
  coord_cartesian(ylim = c(-30,30),xlim=c(180,300),expand=F) +
  # scale_x_continuous(breaks = seq(-180,180,60)) + # npsg
  scale_x_continuous(breaks = seq(180,300,30)) + # centre
  scale_y_continuous(breaks = seq(-30,30,10)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-3"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))

savename <- paste0("ep_winter_95.png")
ggsave(savename, width=8.27, height=3.44, dpi=300)

########## Part I: Summer NP ##########
load("/Volumes/Doris/OCCCI/quantile_regression/summer_regions_quant.Rdata")
mapp_df <- fortify(map("world2",plot=FALSE,fill=TRUE))
com <- rev(c("firebrick4","firebrick2","white","dodgerblue2","dodgerblue4"))#color palette

##### 50
trend_df <- np_centre_df_50 # df change
trend_df$trend <- trend_df$trend*100
trend_df$sigtrend <- NA
ind <- which(trend_df$pvalue <= 0.05)
trend_df$sigtrend[ind] <- "Yes"
ind <- which(is.na(trend_df$pvalue))
trend_df$sigtrend[ind] <- "No" # NA
trend_df$sigtrend <- as.factor(trend_df$sigtrend)

plot <- ggplot() +
  geom_raster(data = trend_df, aes(x = lon,y = lat,fill = trend)) +
  scale_fill_gradientn(colors = com,limits =c(-2,2),na.value = 'black',oob=squish) +
  geom_point(data =trend_df, aes(x = lon,y = lat,col=sigtrend),size=0.8,show.legend = F) +
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) +
  geom_polygon( data=mapp_df, aes(x=long, y=lat, group = group),colour="gainsboro", fill="gainsboro",size=0.25 ) +
  coord_cartesian(ylim = c(30,60),xlim=c(120,240),expand=F) +
  # scale_x_continuous(breaks = seq(-180,180,60)) + # npsg
  scale_x_continuous(breaks = seq(120,240,30)) + # centre
  scale_y_continuous(breaks = seq(30,60,10)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-2"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))
savename <- paste0("np_summer_50.png")
ggsave(savename, width=8.27, height=3.44, dpi=300)

##### 95th
trend_df <- np_centre_df_95
trend_df$trend <- trend_df$trend*100
trend_df$sigtrend <- NA
ind <- which(trend_df$pvalue <= 0.05)
trend_df$sigtrend[ind] <- "Yes"
ind <- which(is.na(trend_df$pvalue))
trend_df$sigtrend[ind] <- "No" # NA
trend_df$sigtrend <- as.factor(trend_df$sigtrend)

plot <- ggplot() +
  geom_raster(data = trend_df, aes(x = lon,y = lat,fill = trend)) +
  scale_fill_gradientn(colors = com,limits =c(-2,2),na.value = 'black',oob=squish) +
  geom_point(data =trend_df, aes(x = lon,y = lat,col=sigtrend),size=0.8,show.legend = F) +
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) +
  geom_polygon( data=mapp_df, aes(x=long, y=lat, group = group),colour="gainsboro", fill="gainsboro",size=0.25 ) +
  coord_cartesian(ylim = c(30,60),xlim=c(120,240),expand=F) +
  # scale_x_continuous(breaks = seq(-180,180,60)) + # npsg
  scale_x_continuous(breaks = seq(120,240,30)) + # centre
  scale_y_continuous(breaks = seq(30,60,10)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-2"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))

savename <- paste0("np_summer_95.png")
ggsave(savename, width=8.27, height=3.44, dpi=300)

########## Part II: Summer NA ##########
rm(list=ls())
load("/Volumes/Doris/OCCCI/quantile_regression/summer_regions_quant.Rdata")
# check plot
mapp_df <- fortify(map("world2",plot=FALSE,fill=TRUE))
com <- rev(c("firebrick4","firebrick2","white","dodgerblue2","dodgerblue4"))#color palette

##### 50
trend_df <- na_centre_df_50 # df change
trend_df$trend <- trend_df$trend*100
trend_df$sigtrend <- NA
ind <- which(trend_df$pvalue <= 0.05)
trend_df$sigtrend[ind] <- "Yes"
ind <- which(is.na(trend_df$pvalue))
trend_df$sigtrend[ind] <- "No" # NA
trend_df$sigtrend <- as.factor(trend_df$sigtrend)

plot <- ggplot() +
  geom_raster(data = trend_df, aes(x = lon,y = lat,fill = trend)) +
  scale_fill_gradientn(colors = com,limits =c(-2,2),na.value = 'black',oob=squish) +
  geom_point(data =trend_df, aes(x = lon,y = lat,col=sigtrend),size=0.8,show.legend = F) +
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) +
  geom_polygon( data=mapp_df, aes(x=long, y=lat, group = group),colour="gainsboro", fill="gainsboro",size=0.25 ) +
  coord_cartesian(ylim = c(30,60),xlim=c(300,360),expand=F) +
  # scale_x_continuous(breaks = seq(-180,180,60)) + # npsg
  scale_x_continuous(breaks = seq(300,360,30)) + # centre
  scale_y_continuous(breaks = seq(30,60,10)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-3"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))
savename <- paste0("na_summer_50.png")
ggsave(savename, width=8.27, height=3.44, dpi=300)


##### 95th
trend_df <- na_centre_df_95
trend_df$trend <- trend_df$trend*100
trend_df$sigtrend <- NA
ind <- which(trend_df$pvalue <= 0.05)
trend_df$sigtrend[ind] <- "Yes"
ind <- which(is.na(trend_df$pvalue))
trend_df$sigtrend[ind] <- "No" # NA
trend_df$sigtrend <- as.factor(trend_df$sigtrend)

plot <- ggplot() +
  geom_raster(data = trend_df, aes(x = lon,y = lat,fill = trend)) +
  scale_fill_gradientn(colors = com,limits =c(-2,2),na.value = 'black',oob=squish) +
  geom_point(data =trend_df, aes(x = lon,y = lat,col=sigtrend),size=0.8,show.legend = F) +
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) +
  geom_polygon( data=mapp_df, aes(x=long, y=lat, group = group),colour="gainsboro", fill="gainsboro",size=0.25 ) +
  coord_cartesian(ylim = c(30,60),xlim=c(300,360),expand=F) +
  # scale_x_continuous(breaks = seq(-180,180,60)) + # npsg
  scale_x_continuous(breaks = seq(300,360,30)) + # centre
  scale_y_continuous(breaks = seq(30,60,10)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-3"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))

savename <- paste0("na_summer_95.png")
ggsave(savename, width=8.27, height=3.44, dpi=300)

########## Part III: Summer EP ##########
rm(list=ls())
load("/Volumes/Doris/OCCCI/quantile_regression/summer_regions_quant.Rdata")
# check plot
mapp_df <- fortify(map("world2",plot=FALSE,fill=TRUE))
com <- rev(c("firebrick4","firebrick2","white","dodgerblue2","dodgerblue4"))#color palette

##### 50
trend_df <- ep_centre_df_50 # df change
trend_df$trend <- trend_df$trend*100
trend_df$sigtrend <- NA
ind <- which(trend_df$pvalue <= 0.05)
trend_df$sigtrend[ind] <- "Yes"
ind <- which(is.na(trend_df$pvalue))
trend_df$sigtrend[ind] <- "No" # NA
trend_df$sigtrend <- as.factor(trend_df$sigtrend)

plot <- ggplot() +
  geom_raster(data = trend_df, aes(x = lon,y = lat,fill = trend)) +
  scale_fill_gradientn(colors = com,limits =c(-2,2),na.value = 'black',oob=squish) +
  geom_point(data =trend_df, aes(x = lon,y = lat,col=sigtrend),size=0.8,show.legend = F) +
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) +
  geom_polygon( data=mapp_df, aes(x=long, y=lat, group = group),colour="gainsboro", fill="gainsboro",size=0.25 ) +
  coord_cartesian(ylim = c(-30,30),xlim=c(180,300),expand=F) +
  # scale_x_continuous(breaks = seq(-180,180,60)) + # npsg
  scale_x_continuous(breaks = seq(180,300,30)) + # centre
  scale_y_continuous(breaks = seq(-30,30,10)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-3"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))
savename <- paste0("ep_summer_50.png")
ggsave(savename, width=8.27, height=3.44, dpi=300)


##### 95th
trend_df <- ep_centre_df_95
trend_df$trend <- trend_df$trend*100
trend_df$sigtrend <- NA
ind <- which(trend_df$pvalue <= 0.05)
trend_df$sigtrend[ind] <- "Yes"
ind <- which(is.na(trend_df$pvalue))
trend_df$sigtrend[ind] <- "No" # NA
trend_df$sigtrend <- as.factor(trend_df$sigtrend)

plot <- ggplot() +
  geom_raster(data = trend_df, aes(x = lon,y = lat,fill = trend)) +
  scale_fill_gradientn(colors = com,limits =c(-2,2),na.value = 'black',oob=squish) +
  geom_point(data =trend_df, aes(x = lon,y = lat,col=sigtrend),size=0.8,show.legend = F) +
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) +
  geom_polygon( data=mapp_df, aes(x=long, y=lat, group = group),colour="gainsboro", fill="gainsboro",size=0.25 ) +
  coord_cartesian(ylim = c(-30,30),xlim=c(180,300),expand=F) +
  # scale_x_continuous(breaks = seq(-180,180,60)) + # npsg
  scale_x_continuous(breaks = seq(180,300,30)) + # centre
  scale_y_continuous(breaks = seq(-30,30,10)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-3"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))

savename <- paste0("ep_summer_95.png")
ggsave(savename, width=8.27, height=3.44, dpi=300)


########## Part I: Winter Global ##########
load("/Volumes/Doris/OCCCI/quantile_regression/winter_ts_df.Rdata")
wmap<-readOGR(dsn="/Volumes/Doris/Projects/spt_wc/wmap", layer="ne_110m_land")
wmap_df <- fortify(wmap)
wmap_df <-wmap_df[-(which(wmap_df$group==112.2)),]#remove caspian sea
com <- rev(c("firebrick4","firebrick2","white","dodgerblue2","dodgerblue4"))#color palette

##### 50
trend_df <- centre_df_50 # df change
trend_df$trend <- trend_df$trend*1000
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
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) +
  geom_polygon(data=wmap_df, aes(x=long, y=lat, group = group),colour="black", fill="gainsboro",size=0.25)+
  coord_equal() + 
  coord_cartesian(ylim = c(-75, 75),xlim = c(-180,180),expand=F)+
  scale_x_continuous(breaks = seq(-180,180,60)) +
  scale_y_continuous(breaks = seq(-90,90,30)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-3"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))
savename <- paste0("winter_50.png")
ggsave(savename, width=8.27, height=3.44, dpi=300)


##### 95th
trend_df <- centre_df_95
trend_df$trend <- trend_df$trend*1000
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
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) +
  geom_polygon(data=wmap_df, aes(x=long, y=lat, group = group),colour="black", fill="gainsboro",size=0.25)+
  coord_equal() + 
  coord_cartesian(ylim = c(-75, 75),xlim = c(-180,180),expand=F)+
  scale_x_continuous(breaks = seq(-180,180,60)) +
  scale_y_continuous(breaks = seq(-90,90,30)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-3"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))
savename <- paste0("winter_95.png")
ggsave(savename, width=8.27, height=3.44, dpi=300)

########## Part II: Summer Global ##########
load("/Volumes/Doris/OCCCI/quantile_regression/summer_ts_df.Rdata")
wmap<-readOGR(dsn="/Volumes/Doris/Projects/spt_wc/wmap", layer="ne_110m_land")
wmap_df <- fortify(wmap)
wmap_df <-wmap_df[-(which(wmap_df$group==112.2)),]#remove caspian sea
com <- rev(c("firebrick4","firebrick2","white","dodgerblue2","dodgerblue4"))#color palette

##### 50
trend_df <- centre_df_50 # df change
trend_df$trend <- trend_df$trend*1000
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
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) +
  geom_polygon(data=wmap_df, aes(x=long, y=lat, group = group),colour="black", fill="gainsboro",size=0.25)+
  coord_equal() + 
  coord_cartesian(ylim = c(-75, 75),xlim = c(-180,180),expand=F)+
  scale_x_continuous(breaks = seq(-180,180,60)) +
  scale_y_continuous(breaks = seq(-90,90,30)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-3"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))
savename <- paste0("summer_50.png")
ggsave(savename, width=8.27, height=3.44, dpi=300)


##### 95th
trend_df <- centre_df_95
trend_df$trend <- trend_df$trend*1000
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
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) +
  geom_polygon(data=wmap_df, aes(x=long, y=lat, group = group),colour="black", fill="gainsboro",size=0.25)+
  coord_equal() + 
  coord_cartesian(ylim = c(-75, 75),xlim = c(-180,180),expand=F)+
  scale_x_continuous(breaks = seq(-180,180,60)) +
  scale_y_continuous(breaks = seq(-90,90,30)) +
  guides(fill = guide_colorbar(title=expression(paste("Trend  ", "(10"^"-2"," %yr"^" -1",")")),
                               title.position = "top",
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))
savename <- paste0("summer_95.png")
ggsave(savename, width=8.27, height=3.44, dpi=300)
