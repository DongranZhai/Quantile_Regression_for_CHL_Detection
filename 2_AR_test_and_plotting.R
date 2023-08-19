# Notes:
# Part I: Using QF test for significance test.
# Part II: Autocorrelation plot

library(dplyr)
library(quantreg)
library(stats)

###  Part I. QF test: quantile accounting autocorrelation in residual ####
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression/global")
load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")

tr <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
QF <- array(data=NA,dim=c(length(lon),length(lat)))
pValue <- array(data=NA,dim=c(length(lon),length(lat)))
ar <- array(data=NA,dim=c(length(lon),length(lat)))

month <- seq(1,304,1)
p=1 # also represent df
t=length(month)
k=1 # of explanatory variable
tick_na <- 0
tick <- 0 

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(monchl_deszn[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit <- lm(monchl_deszn[i,j,] ~ month,na.action = na.omit)
      b0 <- fit$coefficients[1]
      b1 <- fit$coefficients[2]
      tr[i,j] <- as.vector(fit$coefficients[2])
      
      # residuals in unrestricted auxiliary regression
      res <- fit$residuals
      res.lag1 <- c(0,fit$residuals[-t])
      aux.fit <- lm(res~month+res.lag1)
      v_hat <- aux.fit$residuals
      
      # residuals in restricted auxiliary regression imposing H_0
      aux.h0.fit <- lm(res~month)
      v_tilde <- aux.h0.fit$residuals
      
      # QF test
      QF[i,j] <- (sum(v_tilde^2)-sum(v_hat^2))/(sum(v_hat^2)/(t-p-k))
      pValue[i,j] <- pchisq(QF[i,j],df=1,lower.tail=F)
      
      # AR value
      ar.model <- arima(res,order=c(1,0,0))
      ar[i,j] <- as.vector(ar.model$coef[1])
      
      tick <- tick + 1
      
    } else {
      tick_na <- tick_na + 1 
    }
  }
}

savename <- paste0("qf_test.Rdata")
save(tr,QF,ar,pValue,lon,lat,file = savename)

########### Part II: Autocorrelation plot ###########
# linear model
rm(list=ls())
load("/Volumes/Doris/OCCCI/quantile_regression/qf_test.Rdata")
# land map
wmap<-readOGR(dsn="/Volumes/Doris/Projects/spt_wc/wmap", layer="ne_110m_land")
wmap_df <- fortify(wmap)
wmap_df <-wmap_df[-(which(wmap_df$group==112.2)),]#remove caspian sea
# pallette
com_pallette <- viridis_pal()(100)

ar_raster <- matrix2raster(ar, x = lon, y = lat, layer = 2)
ar_df <- as.data.frame(ar_raster,xy = T)
colnames(ar_df) <-c('lon','lat','ar')
rm('ar_raster')

ptemp <- pValue[,ncol(pValue):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
ar_df$p <- pValuetemp
ar_df$sig<- NA

ind <- which(ar_df$p <= 0.05)
ar_df$sig[ind] <- "Yes"
ind <- which(is.na(ar_df$p))
ar_df$sig[ind] <- "No" # NA
ar_df$sig <- as.factor(ar_df$sig)

plot <- ggplot() +
  geom_raster(data = ar_df, aes(x = lon,y = lat,fill = ar)) +
  scale_fill_gradientn(colors = com_pallette,limits =c(-1,1),na.value = 'black',oob=squish) +
  geom_point(data =ar_df, aes(x = lon,y = lat,col=sig),size=0.00000001,show.legend = F) +
  scale_color_manual(values = c(Yes = alpha("white", 0),No = alpha("#353b48", 0.5))) +
  geom_polygon(data=wmap_df, aes(x=long, y=lat, group = group),colour="black", fill="gainsboro",size=0.25)+
  coord_equal() + 
  coord_cartesian(ylim = c(-75, 75),xlim = c(-180,180),expand=F)+
  scale_x_continuous(breaks = seq(-180,180,60)) +
  scale_y_continuous(breaks = seq(-90,90,30)) +
  guides(fill = guide_colorbar(title=paste('AR coefficient'),title.position = "top",
                               title.theme = element_text(angle = 90),
                               barwidth = unit(.4, "cm"),barheight = unit(3.5, "cm"),title.vjust =1)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = 'right',legend.direction = "vertical",
        legend.title = element_text(size = 10,angle = 90),
        plot.margin = margin(t = .3, r = .3, b = .3, l = .3, unit = "cm"))

ggsave(filename = paste0("qr_ar.png"), width=8.27, height=3.44, dpi=300)

