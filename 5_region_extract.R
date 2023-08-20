# Notes: This script is about extracting regions with high CHL signals,
#        and analyzing regional time series after averaging CHL.
#        I. Extract region based on Longhurst
#        II. Regional simple linear trends: North Pacific
#        III. Regional simple linear trends: North Atlantic
#        IV. Regional simple linear trends: Equatorial Pacific
#        V. Regional simple linear trends: Southern Ocean
#        VI. Regional simple linear trends: North Pacific Subtropical
#        VII. Regional simple linear trends: North Atlantic Subtropical

library(R.matlab)
library(quantreg)
library(raster)
library(oceanmap)

########## Part I. Extract region based on Longhurst ##########
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/month_whitening_occci.Rdata")
monchl_deszn <- monchl_deszn_trans

# Longhurst regions index
temp <- readMat("/Volumes/Doris/Projects/sptmodel/Longhurst_180.mat")
Longhurst <- temp$Longhurst
nlon <- length(lon)
nlat <- length(lat)
area <- sort(unique(Longhurst[!is.na(Longhurst)]))#list of longhurst areas to iterate through
area <- area[-c(1,2,4,6,7,8,10,11,12,13,17,18,19,22,25,26,27,29,30,31,34,35,36,40,43,41,42,47,50,53,54)] #remove coastal, polar, GoM, Archipelagic Deep basin, Mediterranean

##### North Pacific (Longhurst is 18 & 19, which corresponding to 24 & 46 in area)
long_np <- array(data=NA,dim=c(nlon,nlat))
np_deszn <- array(data=NA,dim=c(nlon,nlat,304))

for(i in 1:nlon){
  for(j in 1:nlat){
    if(Longhurst[i,j] %in% c(24,46)){
      long_np[i,j] <- Longhurst[i,j]
    }else{
      long_np[i,j] <- NA
    }
  }
}
ind <- which(is.na(long_np))

for(i in c(1:304)){
  tempchl <- monchl_deszn[,,i]
  tempchl[ind] <- NA
  np_deszn[,,i] <- tempchl
}

##### Equatorial Pacific (Longhurst is 5 & 7, which corresponding to 28 & 14 in area)
long_ep <- array(data=NA,dim=c(nlon,nlat))
ep_deszn <- array(data=NA,dim=c(nlon,nlat,304))

for(i in 1:nlon){
  for(j in 1:nlat){
    if(Longhurst[i,j] %in% c(28,14)){
      long_ep[i,j] <- Longhurst[i,j]
    }else{
      long_ep[i,j] <- NA
    }
  }
}
ind <- which(is.na(long_ep))

for(i in c(1:304)){
  tempchl <- monchl_deszn[,,i]
  tempchl[ind] <- NA
  ep_deszn[,,i] <- tempchl
}

##### North Atlantic (Longhurst is 13, which corresponding to 39 in area)
long_na <- array(data=NA,dim=c(nlon,nlat))
na_deszn <- array(data=NA,dim=c(nlon,nlat,304))

for(i in 1:nlon){
  for(j in 1:nlat){
    if(Longhurst[i,j] %in% c(39)){
      long_na[i,j] <- Longhurst[i,j]
    }else{
      long_na[i,j] <- NA
    }
  }
}
ind <- which(is.na(long_na))

for(i in c(1:304)){
  tempchl <- monchl_deszn[,,i]
  tempchl[ind] <- NA
  na_deszn[,,i] <- tempchl
}

##### Southern Ocean (Longhurst is 22, which corresponding to 51 in area)
long_so <- array(data=NA,dim=c(nlon,nlat))
so_deszn <- array(data=NA,dim=c(nlon,nlat,304))

for(i in 1:nlon){
  for(j in 1:nlat){
    if(Longhurst[i,j] %in% c(51)){
      long_so[i,j] <- Longhurst[i,j]
    }else{
      long_so[i,j] <- NA
    }
  }
}
ind <- which(is.na(long_so))

for(i in c(1:304)){
  tempchl <- monchl_deszn[,,i]
  tempchl[ind] <- NA
  so_deszn[,,i] <- tempchl
}

##### North Pacific Subtropical (Longhurst is 17, which corresponding to 49 in area)
long_npo <- array(data=NA,dim=c(nlon,nlat))
npo_deszn <- array(data=NA,dim=c(nlon,nlat,304))

for(i in 1:nlon){
  for(j in 1:nlat){
    if(Longhurst[i,j] %in% c(49)){
      long_npo[i,j] <- Longhurst[i,j]
    }else{
      long_npo[i,j] <- NA
    }
  }
}
ind <- which(is.na(long_npo))

for(i in c(1:304)){
  tempchl <- monchl_deszn[,,i]
  tempchl[ind] <- NA
  npo_deszn[,,i] <- tempchl
}

##### North Atlantic Subtropical (Longhurst is 15, which corresponding to 37 in area)
long_nao <- array(data=NA,dim=c(nlon,nlat))
nao_deszn <- array(data=NA,dim=c(nlon,nlat,304))

for(i in 1:nlon){
  for(j in 1:nlat){
    if(Longhurst[i,j] %in% c(37)){
      long_nao[i,j] <- Longhurst[i,j]
    }else{
      long_nao[i,j] <- NA
    }
  }
}
ind <- which(is.na(long_nao))

for(i in c(1:304)){
  tempchl <- monchl_deszn[,,i]
  tempchl[ind] <- NA
  nao_deszn[,,i] <- tempchl
}
save(long_np,np_deszn,long_ep,ep_deszn,long_na,na_deszn,long_so,so_deszn,
     long_npo,npo_deszn,long_nao,nao_deszn,
     lon,lat,file='regions_chl_deszn_trans.Rdata')

########## Part II. Regional simple linear trends: North Pacific ##########
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/regions_chl_deszn_trans.Rdata")
month <- c(1:304)
# center Pacific
np_centre_df<- data.frame(matrix(ncol = 4, nrow = 1)) # array(matrix) to raster
colnames(np_centre_df) <- c('x','y','chl_deszn','layer')
for(i in month){
  print(i)
  temp_r <- matrix2raster(np_deszn[,,i], x = lon, y = lat, layer = i)
  
  temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
  temp2 <- crop(temp_r, extent(0, 180, -90, 90))
  extent(temp1) <- c(180, 360, -90, 90)
  
  temp <- merge(temp1, temp2)
  temp <- as.data.frame(temp,xy=T)
  temp$chl_deszn <- i
  np_centre_df <- rbind(np_centre_df,temp) # np_centre is a df with Pacific as centre
}
colnames(np_centre_df) <- c('x','y','layer','chl_deszn')
np_centre_df <- np_centre_df[-1,]
np_centre <- array(data=NA,dim=c(360,180,304))

for(i in month){
  print(i)
  temp <- np_centre_df[which(np_centre_df$layer== i),c(1:2,4)]
  coordinates(temp) <- ~ x + y
  gridded(temp) <- T
  temp_r <- raster(temp)
  np_centre[,,i] <- raster2matrix(temp_r)
}
# save(np_centre,np_centre_df,file='np_chl_centre.Rdata')

# trends fit
tr.lm <- array(data=NA,dim=c(length(lon),length(lat)))
tr.01 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.05 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.10 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.50 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.90 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.95 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.99 <- array(data=NA,dim=c(length(lon),length(lat)))
p.lm <- array(data=NA,dim=c(length(lon),length(lat)))
p.01 <- array(data=NA,dim=c(length(lon),length(lat)))
p.05 <- array(data=NA,dim=c(length(lon),length(lat)))
p.10 <- array(data=NA,dim=c(length(lon),length(lat)))
p.50 <- array(data=NA,dim=c(length(lon),length(lat)))
p.90 <- array(data=NA,dim=c(length(lon),length(lat)))
p.95 <- array(data=NA,dim=c(length(lon),length(lat)))
p.99 <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
tick_na <- 0 
tick <- 0 

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(np_centre[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit.lm <- lm(np_centre[i,j,] ~ month,na.action = na.omit)
      fit.01 <- rq(np_centre[i,j,] ~ month,tau=.01,na.action = na.omit)
      fit.05 <- rq(np_centre[i,j,] ~ month,tau=.05,na.action = na.omit)
      fit.10 <- rq(np_centre[i,j,] ~ month,tau=.10,na.action = na.omit)
      fit.50 <- rq(np_centre[i,j,] ~ month,tau=.50,na.action = na.omit)
      fit.90 <- rq(np_centre[i,j,] ~ month,tau=.90,na.action = na.omit)
      fit.95 <- rq(np_centre[i,j,] ~ month,tau=.95,na.action = na.omit)
      fit.99 <- rq(np_centre[i,j,] ~ month,tau=.99,na.action = na.omit)
      tr.lm[i,j] <- fit.lm$coefficients[2]
      tr.01[i,j] <- fit.01$coefficients[2]
      tr.05[i,j] <- fit.05$coefficients[2]
      tr.10[i,j] <- fit.10$coefficients[2]
      tr.50[i,j] <- fit.50$coefficients[2]
      tr.90[i,j] <- fit.90$coefficients[2]
      tr.95[i,j] <- fit.95$coefficients[2]
      tr.99[i,j] <- fit.99$coefficients[2]
      p.lm[i,j] <- summary(fit.lm)$coefficients[8]
      p.01[i,j] <- summary(fit.01,se='ker')$coefficients[8]
      p.05[i,j] <- summary(fit.05,se='ker')$coefficients[8]
      p.10[i,j] <- summary(fit.10,se='ker')$coefficients[8]
      p.50[i,j] <- summary(fit.50,se='ker')$coefficients[8]
      p.90[i,j] <- summary(fit.90,se='ker')$coefficients[8]
      p.95[i,j] <- summary(fit.95,se='ker')$coefficients[8]
      p.99[i,j] <- summary(fit.99,se='ker')$coefficients[8]
      tick <- tick + 1 
      
    } else {
      tick_na <- tick_na + 1
    }
  }
}

# trend value
trendvalue <- tr.lm*1000
trend_raster <- matrix2raster(trendvalue, x = seq(1,360,1), y = lat, layer = 2)
trend_df <- as.data.frame(trend_raster,xy = T)
colnames(trend_df) <-c('lon','lat','tr.lm')
rm(list=c('trend_raster'))

trendvalue <- tr.01[,ncol(tr.01):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.01 <- trendtemp

trendvalue <- tr.05[,ncol(tr.05):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.05 <- trendtemp

trendvalue <- tr.10[,ncol(tr.10):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.10 <- trendtemp

trendvalue <- tr.50[,ncol(tr.50):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.50 <- trendtemp

trendvalue <- tr.90[,ncol(tr.90):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.90 <- trendtemp

trendvalue <- tr.95[,ncol(tr.95):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.95 <- trendtemp

trendvalue <- tr.99[,ncol(tr.99):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.99 <- trendtemp

# p value
ptemp <- p.lm[,ncol(p.lm):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.lm <- NA
ind <- which(pValuetemp <= 0.05)
trend_df$sig.lm[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.lm[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.lm[ind] <- "No" # NA
trend_df$sig.lm <- as.factor(trend_df$sig.lm)

ptemp <- p.01[,ncol(p.01):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.01 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.01[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.01[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.01[ind] <- "No" # NA
trend_df$sig.01 <- as.factor(trend_df$sig.01)

ptemp <- p.05[,ncol(p.05):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.05 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.05[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.05[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.05[ind] <- "No" # NA
trend_df$sig.05 <- as.factor(trend_df$sig.05)

ptemp <- p.10[,ncol(p.10):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.10 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.10[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.10[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.10[ind] <- "No" # NA
trend_df$sig.10 <- as.factor(trend_df$sig.10)

ptemp <- p.50[,ncol(p.50):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.50 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.50[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.50[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.50[ind] <- "No" # NA
trend_df$sig.50 <- as.factor(trend_df$sig.50)

ptemp <- p.90[,ncol(p.90):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.90 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.90[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.90[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.90[ind] <- "No" # NA
trend_df$sig.90 <- as.factor(trend_df$sig.90)

ptemp <- p.95[,ncol(p.95):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.95 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.95[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.95[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.95[ind] <- "No" # NA
trend_df$sig.95 <- as.factor(trend_df$sig.95)

ptemp <- p.99[,ncol(p.99):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.99 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.99[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.99[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.99[ind] <- "No" # NA
trend_df$sig.99 <- as.factor(trend_df$sig.99)

save(tr.lm,tr.01,tr.05,tr.50,tr.95,tr.99,tr.10,tr.90,
     p.lm,p.01,p.05,p.50,p.95,p.99,p.10,p.90,
     np_centre_df,np_centre,
     trend_df,np_deszn,lon,lat,file='np_fit_panel.Rdata')

########## Part III. Regional simple linear trends: North Atlantic ##########
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/regions_chl_deszn_trans.Rdata")
# setwd('/home/dzhai/quantile_regression')
# load('/home/dzhai/quantile_regression/regions_chl_deszn.Rdata')
month <- c(1:304)

# trends fit (no need centre)
tr.lm <- array(data=NA,dim=c(length(lon),length(lat)))
tr.01 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.05 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.10 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.50 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.90 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.95 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.99 <- array(data=NA,dim=c(length(lon),length(lat)))
p.lm <- array(data=NA,dim=c(length(lon),length(lat)))
p.01 <- array(data=NA,dim=c(length(lon),length(lat)))
p.05 <- array(data=NA,dim=c(length(lon),length(lat)))
p.10 <- array(data=NA,dim=c(length(lon),length(lat)))
p.50 <- array(data=NA,dim=c(length(lon),length(lat)))
p.90 <- array(data=NA,dim=c(length(lon),length(lat)))
p.95 <- array(data=NA,dim=c(length(lon),length(lat)))
p.99 <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
tick_na <- 0 
tick <- 0 

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(na_deszn[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit.lm <- lm(na_deszn[i,j,] ~ month,na.action = na.omit)
      fit.01 <- rq(na_deszn[i,j,] ~ month,tau=.01,na.action = na.omit)
      fit.05 <- rq(na_deszn[i,j,] ~ month,tau=.05,na.action = na.omit)
      fit.10 <- rq(na_deszn[i,j,] ~ month,tau=.10,na.action = na.omit)
      fit.50 <- rq(na_deszn[i,j,] ~ month,tau=.50,na.action = na.omit)
      fit.90 <- rq(na_deszn[i,j,] ~ month,tau=.90,na.action = na.omit)
      fit.95 <- rq(na_deszn[i,j,] ~ month,tau=.95,na.action = na.omit)
      fit.99 <- rq(na_deszn[i,j,] ~ month,tau=.99,na.action = na.omit)
      tr.lm[i,j] <- fit.lm$coefficients[2]
      tr.01[i,j] <- fit.01$coefficients[2]
      tr.05[i,j] <- fit.05$coefficients[2]
      tr.10[i,j] <- fit.10$coefficients[2]
      tr.50[i,j] <- fit.50$coefficients[2]
      tr.90[i,j] <- fit.90$coefficients[2]
      tr.95[i,j] <- fit.95$coefficients[2]
      tr.99[i,j] <- fit.99$coefficients[2]
      p.lm[i,j] <- summary(fit.lm)$coefficients[8]
      p.01[i,j] <- summary(fit.01,se='ker')$coefficients[8]
      p.05[i,j] <- summary(fit.05,se='ker')$coefficients[8]
      p.10[i,j] <- summary(fit.10,se='ker')$coefficients[8]
      p.50[i,j] <- summary(fit.50,se='ker')$coefficients[8]
      p.90[i,j] <- summary(fit.90,se='ker')$coefficients[8]
      p.95[i,j] <- summary(fit.95,se='ker')$coefficients[8]
      p.99[i,j] <- summary(fit.99,se='ker')$coefficients[8]
      tick <- tick + 1 
      
    } else {
      tick_na <- tick_na + 1
    }
  }
}

# trend value
trendvalue <- tr.lm*1000
trend_raster <- matrix2raster(trendvalue, x = lon, y = lat, layer = 2) # no need centre
trend_df <- as.data.frame(trend_raster,xy = T)
colnames(trend_df) <-c('lon','lat','tr.lm')
rm(list=c('trend_raster'))

trendvalue <- tr.01[,ncol(tr.01):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.01 <- trendtemp

trendvalue <- tr.05[,ncol(tr.05):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.05 <- trendtemp

trendvalue <- tr.10[,ncol(tr.10):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.10 <- trendtemp

trendvalue <- tr.50[,ncol(tr.50):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.50 <- trendtemp

trendvalue <- tr.90[,ncol(tr.90):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.90 <- trendtemp

trendvalue <- tr.95[,ncol(tr.95):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.95 <- trendtemp

trendvalue <- tr.99[,ncol(tr.99):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.99 <- trendtemp

# p value
ptemp <- p.lm[,ncol(p.lm):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.lm <- NA
ind <- which(pValuetemp <= 0.05)
trend_df$sig.lm[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.lm[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.lm[ind] <- "No" # NA
trend_df$sig.lm <- as.factor(trend_df$sig.lm)

ptemp <- p.01[,ncol(p.01):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.01 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.01[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.01[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.01[ind] <- "No" # NA
trend_df$sig.01 <- as.factor(trend_df$sig.01)

ptemp <- p.05[,ncol(p.05):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.05 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.05[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.05[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.05[ind] <- "No" # NA
trend_df$sig.05 <- as.factor(trend_df$sig.05)

ptemp <- p.10[,ncol(p.10):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.10 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.10[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.10[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.10[ind] <- "No" # NA
trend_df$sig.10 <- as.factor(trend_df$sig.10)

ptemp <- p.50[,ncol(p.50):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.50 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.50[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.50[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.50[ind] <- "No" # NA
trend_df$sig.50 <- as.factor(trend_df$sig.50)

ptemp <- p.90[,ncol(p.90):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.90 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.90[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.90[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.90[ind] <- "No" # NA
trend_df$sig.90 <- as.factor(trend_df$sig.90)

ptemp <- p.95[,ncol(p.95):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.95 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.95[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.95[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.95[ind] <- "No" # NA
trend_df$sig.95 <- as.factor(trend_df$sig.95)

ptemp <- p.99[,ncol(p.99):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.99 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.99[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.99[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.99[ind] <- "No" # NA
trend_df$sig.99 <- as.factor(trend_df$sig.99)

save(tr.lm,tr.01,tr.05,tr.50,tr.95,tr.99,tr.10,tr.90,
     p.lm,p.01,p.05,p.50,p.95,p.99,p.10,p.90,
     trend_df,na_deszn,lon,lat,file='na_fit_panel.Rdata')

########## Part IV. Regional simple linear trends: Equatorial Pacific ##########
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/regions_chl_deszn_trans.Rdata")
# setwd('/home/dzhai/quantile_regression')
# load('/home/dzhai/quantile_regression/regions_chl_deszn.Rdata')
month <- c(1:304)
# center Pacific
ep_centre_df<- data.frame(matrix(ncol = 4, nrow = 1)) # array(matrix) to raster
colnames(ep_centre_df) <- c('x','y','chl_deszn','layer')

for(i in month){
  print(i)
  temp_r <- matrix2raster(ep_deszn[,,i], x = lon, y = lat, layer = i)
  
  temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
  temp2 <- crop(temp_r, extent(0, 180, -90, 90))
  extent(temp1) <- c(180, 360, -90, 90)
  
  temp <- merge(temp1, temp2)
  temp <- as.data.frame(temp,xy=T)
  temp$chl_deszn <- i
  ep_centre_df <- rbind(ep_centre_df,temp) # ep_centre is a df with Pacific as centre
}
colnames(ep_centre_df) <- c('x','y','layer','chl_deszn')
ep_centre_df <- ep_centre_df[-1,]
ep_centre <- array(data=NA,dim=c(360,180,304))

for(i in month){
  print(i)
  temp <- ep_centre_df[which(ep_centre_df$layer== i),c(1:2,4)]
  coordinates(temp) <- ~ x + y
  gridded(temp) <- T
  temp_r <- raster(temp)
  ep_centre[,,i] <- raster2matrix(temp_r)
}

# trends fit
tr.lm <- array(data=NA,dim=c(length(lon),length(lat)))
tr.01 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.05 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.10 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.50 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.90 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.95 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.99 <- array(data=NA,dim=c(length(lon),length(lat)))
p.lm <- array(data=NA,dim=c(length(lon),length(lat)))
p.01 <- array(data=NA,dim=c(length(lon),length(lat)))
p.05 <- array(data=NA,dim=c(length(lon),length(lat)))
p.10 <- array(data=NA,dim=c(length(lon),length(lat)))
p.50 <- array(data=NA,dim=c(length(lon),length(lat)))
p.90 <- array(data=NA,dim=c(length(lon),length(lat)))
p.95 <- array(data=NA,dim=c(length(lon),length(lat)))
p.99 <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
tick_na <- 0 
tick <- 0 

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(ep_centre[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit.lm <- lm(ep_centre[i,j,] ~ month,na.action = na.omit)
      fit.01 <- rq(ep_centre[i,j,] ~ month,tau=.01,na.action = na.omit)
      fit.05 <- rq(ep_centre[i,j,] ~ month,tau=.05,na.action = na.omit)
      fit.10 <- rq(ep_centre[i,j,] ~ month,tau=.10,na.action = na.omit)
      fit.50 <- rq(ep_centre[i,j,] ~ month,tau=.50,na.action = na.omit)
      fit.90 <- rq(ep_centre[i,j,] ~ month,tau=.90,na.action = na.omit)
      fit.95 <- rq(ep_centre[i,j,] ~ month,tau=.95,na.action = na.omit)
      fit.99 <- rq(ep_centre[i,j,] ~ month,tau=.99,na.action = na.omit)
      tr.lm[i,j] <- fit.lm$coefficients[2]
      tr.01[i,j] <- fit.01$coefficients[2]
      tr.05[i,j] <- fit.05$coefficients[2]
      tr.10[i,j] <- fit.10$coefficients[2]
      tr.50[i,j] <- fit.50$coefficients[2]
      tr.90[i,j] <- fit.90$coefficients[2]
      tr.95[i,j] <- fit.95$coefficients[2]
      tr.99[i,j] <- fit.99$coefficients[2]
      p.lm[i,j] <- summary(fit.lm)$coefficients[8]
      p.01[i,j] <- summary(fit.01,se='ker')$coefficients[8]
      p.05[i,j] <- summary(fit.05,se='ker')$coefficients[8]
      p.10[i,j] <- summary(fit.10,se='ker')$coefficients[8]
      p.50[i,j] <- summary(fit.50,se='ker')$coefficients[8]
      p.90[i,j] <- summary(fit.90,se='ker')$coefficients[8]
      p.95[i,j] <- summary(fit.95,se='ker')$coefficients[8]
      p.99[i,j] <- summary(fit.99,se='ker')$coefficients[8]
      tick <- tick + 1 
      
    } else {
      tick_na <- tick_na + 1
    }
  }
}

# trend value
trendvalue <- tr.lm*1000
trend_raster <- matrix2raster(trendvalue, x = seq(1,360,1), y = lat, layer = 2)
trend_df <- as.data.frame(trend_raster,xy = T)
colnames(trend_df) <-c('lon','lat','tr.lm')
rm(list=c('trend_raster'))

trendvalue <- tr.01[,ncol(tr.01):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.01 <- trendtemp

trendvalue <- tr.05[,ncol(tr.05):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.05 <- trendtemp

trendvalue <- tr.10[,ncol(tr.10):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.10 <- trendtemp

trendvalue <- tr.50[,ncol(tr.50):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.50 <- trendtemp

trendvalue <- tr.90[,ncol(tr.90):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.90 <- trendtemp

trendvalue <- tr.95[,ncol(tr.95):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.95 <- trendtemp

trendvalue <- tr.99[,ncol(tr.99):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.99 <- trendtemp

# p value
ptemp <- p.lm[,ncol(p.lm):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.lm <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.lm[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.lm[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.lm[ind] <- "No" # NA
trend_df$sig.lm <- as.factor(trend_df$sig.lm)

ptemp <- p.01[,ncol(p.01):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.01 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.01[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.01[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.01[ind] <- "No" # NA
trend_df$sig.01 <- as.factor(trend_df$sig.01)

ptemp <- p.05[,ncol(p.05):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.05 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.05[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.05[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.05[ind] <- "No" # NA
trend_df$sig.05 <- as.factor(trend_df$sig.05)

ptemp <- p.10[,ncol(p.10):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.10 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.10[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.10[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.10[ind] <- "No" # NA
trend_df$sig.10 <- as.factor(trend_df$sig.10)

ptemp <- p.50[,ncol(p.50):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.50 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.50[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.50[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.50[ind] <- "No" # NA
trend_df$sig.50 <- as.factor(trend_df$sig.50)

ptemp <- p.90[,ncol(p.90):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.90 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.90[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.90[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.90[ind] <- "No" # NA
trend_df$sig.90 <- as.factor(trend_df$sig.90)

ptemp <- p.95[,ncol(p.95):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.95 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.95[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.95[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.95[ind] <- "No" # NA
trend_df$sig.95 <- as.factor(trend_df$sig.95)

ptemp <- p.99[,ncol(p.99):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.99 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.99[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.99[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.99[ind] <- "No" # NA
trend_df$sig.99 <- as.factor(trend_df$sig.99)

save(tr.lm,tr.01,tr.05,tr.50,tr.95,tr.99,tr.10,tr.90,
     p.lm,p.01,p.05,p.50,p.95,p.99,p.10,p.90,
     ep_centre_df,ep_centre,
     trend_df,ep_deszn,file='ep_fit_panel.Rdata')

########## Part V. Regional simple linear trends: Southern Ocean ########## 
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/regions_chl_deszn_trans.Rdata")
# setwd('/home/dzhai/quantile_regression')
# load('/home/dzhai/quantile_regression/regions_chl_deszn.Rdata')
month <- c(1:304)

# trends fit (no need centre)
tr.lm <- array(data=NA,dim=c(length(lon),length(lat)))
tr.01 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.05 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.10 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.50 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.90 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.95 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.99 <- array(data=NA,dim=c(length(lon),length(lat)))
p.lm <- array(data=NA,dim=c(length(lon),length(lat)))
p.01 <- array(data=NA,dim=c(length(lon),length(lat)))
p.05 <- array(data=NA,dim=c(length(lon),length(lat)))
p.10 <- array(data=NA,dim=c(length(lon),length(lat)))
p.50 <- array(data=NA,dim=c(length(lon),length(lat)))
p.90 <- array(data=NA,dim=c(length(lon),length(lat)))
p.95 <- array(data=NA,dim=c(length(lon),length(lat)))
p.99 <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
tick_na <- 0 
tick <- 0 

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(so_deszn[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit.lm <- lm(so_deszn[i,j,] ~ month,na.action = na.omit)
      fit.01 <- rq(so_deszn[i,j,] ~ month,tau=.01,na.action = na.omit)
      fit.05 <- rq(so_deszn[i,j,] ~ month,tau=.05,na.action = na.omit)
      fit.10 <- rq(so_deszn[i,j,] ~ month,tau=.10,na.action = na.omit)
      fit.50 <- rq(so_deszn[i,j,] ~ month,tau=.50,na.action = na.omit)
      fit.90 <- rq(so_deszn[i,j,] ~ month,tau=.90,na.action = na.omit)
      fit.95 <- rq(so_deszn[i,j,] ~ month,tau=.95,na.action = na.omit)
      fit.99 <- rq(so_deszn[i,j,] ~ month,tau=.99,na.action = na.omit)
      tr.lm[i,j] <- fit.lm$coefficients[2]
      tr.01[i,j] <- fit.01$coefficients[2]
      tr.05[i,j] <- fit.05$coefficients[2]
      tr.10[i,j] <- fit.10$coefficients[2]
      tr.50[i,j] <- fit.50$coefficients[2]
      tr.90[i,j] <- fit.90$coefficients[2]
      tr.95[i,j] <- fit.95$coefficients[2]
      tr.99[i,j] <- fit.99$coefficients[2]
      p.lm[i,j] <- summary(fit.lm)$coefficients[8]
      p.01[i,j] <- summary(fit.01,se='ker')$coefficients[8]
      p.05[i,j] <- summary(fit.05,se='ker')$coefficients[8]
      p.10[i,j] <- summary(fit.10,se='ker')$coefficients[8]
      p.50[i,j] <- summary(fit.50,se='ker')$coefficients[8]
      p.90[i,j] <- summary(fit.90,se='ker')$coefficients[8]
      p.95[i,j] <- summary(fit.95,se='ker')$coefficients[8]
      p.99[i,j] <- summary(fit.99,se='ker')$coefficients[8]
      tick <- tick + 1
      
    } else {
      tick_na <- tick_na + 1
    }
  }
}

# trend value
trendvalue <- tr.lm*1000
trend_raster <- matrix2raster(trendvalue, x = lon, y = lat, layer = 2) # no need centre
trend_df <- as.data.frame(trend_raster,xy = T)
colnames(trend_df) <-c('lon','lat','tr.lm')
rm(list=c('trend_raster'))

trendvalue <- tr.01[,ncol(tr.01):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.01 <- trendtemp

trendvalue <- tr.05[,ncol(tr.05):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.05 <- trendtemp

trendvalue <- tr.10[,ncol(tr.10):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.10 <- trendtemp

trendvalue <- tr.50[,ncol(tr.50):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.50 <- trendtemp

trendvalue <- tr.90[,ncol(tr.90):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.90 <- trendtemp

trendvalue <- tr.95[,ncol(tr.95):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.95 <- trendtemp

trendvalue <- tr.99[,ncol(tr.99):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.99 <- trendtemp

# p value
ptemp <- p.lm[,ncol(p.lm):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.lm <- NA
ind <- which(pValuetemp <= 0.05)
trend_df$sig.lm[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.lm[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.lm[ind] <- "No" # NA
trend_df$sig.lm <- as.factor(trend_df$sig.lm)

ptemp <- p.01[,ncol(p.01):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.01 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.01[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.01[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.01[ind] <- "No" # NA
trend_df$sig.01 <- as.factor(trend_df$sig.01)

ptemp <- p.05[,ncol(p.05):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.05 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.05[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.05[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.05[ind] <- "No" # NA
trend_df$sig.05 <- as.factor(trend_df$sig.05)

ptemp <- p.10[,ncol(p.10):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.10 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.10[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.10[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.10[ind] <- "No" # NA
trend_df$sig.10 <- as.factor(trend_df$sig.10)

ptemp <- p.50[,ncol(p.50):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.50 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.50[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.50[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.50[ind] <- "No" # NA
trend_df$sig.50 <- as.factor(trend_df$sig.50)

ptemp <- p.90[,ncol(p.90):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.90 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.90[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.90[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.90[ind] <- "No" # NA
trend_df$sig.90 <- as.factor(trend_df$sig.90)

ptemp <- p.95[,ncol(p.95):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.95 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.95[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.95[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.95[ind] <- "No" # NA
trend_df$sig.95 <- as.factor(trend_df$sig.95)

ptemp <- p.99[,ncol(p.99):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.99 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.99[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.99[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.99[ind] <- "No" # NA
trend_df$sig.99 <- as.factor(trend_df$sig.99)

save(tr.lm,tr.01,tr.05,tr.50,tr.95,tr.99,tr.10,tr.90,
     p.lm,p.01,p.05,p.50,p.95,p.99,p.10,p.90,
     trend_df,so_deszn,lon,lat,file='so_fit_panel.Rdata')


########## Part VI. Regional simple linear trends: North Pacific Subtropical ########## 
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/regions_chl_deszn_trans.Rdata")
# setwd('/home/dzhai/quantile_regression')
# load('/home/dzhai/quantile_regression/regions_chl_deszn.Rdata')
month <- c(1:304)

# trends fit (no need centre)
tr.lm <- array(data=NA,dim=c(length(lon),length(lat)))
tr.01 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.05 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.10 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.50 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.90 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.95 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.99 <- array(data=NA,dim=c(length(lon),length(lat)))
p.lm <- array(data=NA,dim=c(length(lon),length(lat)))
p.01 <- array(data=NA,dim=c(length(lon),length(lat)))
p.05 <- array(data=NA,dim=c(length(lon),length(lat)))
p.10 <- array(data=NA,dim=c(length(lon),length(lat)))
p.50 <- array(data=NA,dim=c(length(lon),length(lat)))
p.90 <- array(data=NA,dim=c(length(lon),length(lat)))
p.95 <- array(data=NA,dim=c(length(lon),length(lat)))
p.99 <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
tick_na <- 0 
tick <- 0 

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(npo_deszn[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit.lm <- lm(npo_deszn[i,j,] ~ month,na.action = na.omit)
      fit.01 <- rq(npo_deszn[i,j,] ~ month,tau=.01,na.action = na.omit)
      fit.05 <- rq(npo_deszn[i,j,] ~ month,tau=.05,na.action = na.omit)
      fit.10 <- rq(npo_deszn[i,j,] ~ month,tau=.10,na.action = na.omit)
      fit.50 <- rq(npo_deszn[i,j,] ~ month,tau=.50,na.action = na.omit)
      fit.90 <- rq(npo_deszn[i,j,] ~ month,tau=.90,na.action = na.omit)
      fit.95 <- rq(npo_deszn[i,j,] ~ month,tau=.95,na.action = na.omit)
      fit.99 <- rq(npo_deszn[i,j,] ~ month,tau=.99,na.action = na.omit)
      tr.lm[i,j] <- fit.lm$coefficients[2]
      tr.01[i,j] <- fit.01$coefficients[2]
      tr.05[i,j] <- fit.05$coefficients[2]
      tr.10[i,j] <- fit.10$coefficients[2]
      tr.50[i,j] <- fit.50$coefficients[2]
      tr.90[i,j] <- fit.90$coefficients[2]
      tr.95[i,j] <- fit.95$coefficients[2]
      tr.99[i,j] <- fit.99$coefficients[2]
      p.lm[i,j] <- summary(fit.lm)$coefficients[8]
      p.01[i,j] <- summary(fit.01,se='ker')$coefficients[8]
      p.05[i,j] <- summary(fit.05,se='ker')$coefficients[8]
      p.10[i,j] <- summary(fit.10,se='ker')$coefficients[8]
      p.50[i,j] <- summary(fit.50,se='ker')$coefficients[8]
      p.90[i,j] <- summary(fit.90,se='ker')$coefficients[8]
      p.95[i,j] <- summary(fit.95,se='ker')$coefficients[8]
      p.99[i,j] <- summary(fit.99,se='ker')$coefficients[8]
      tick <- tick + 1
      
    } else {
      tick_na <- tick_na + 1
    }
  }
}

# trend value
trendvalue <- tr.lm*1000
trend_raster <- matrix2raster(trendvalue, x = lon, y = lat, layer = 2) # no need centre
trend_df <- as.data.frame(trend_raster,xy = T)
colnames(trend_df) <-c('lon','lat','tr.lm')
rm(list=c('trend_raster'))

trendvalue <- tr.01[,ncol(tr.01):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.01 <- trendtemp

trendvalue <- tr.05[,ncol(tr.05):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.05 <- trendtemp

trendvalue <- tr.10[,ncol(tr.10):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.10 <- trendtemp

trendvalue <- tr.50[,ncol(tr.50):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.50 <- trendtemp

trendvalue <- tr.90[,ncol(tr.90):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.90 <- trendtemp

trendvalue <- tr.95[,ncol(tr.95):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.95 <- trendtemp

trendvalue <- tr.99[,ncol(tr.99):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.99 <- trendtemp

# p value
ptemp <- p.lm[,ncol(p.lm):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.lm <- NA
ind <- which(pValuetemp <= 0.05)
trend_df$sig.lm[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.lm[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.lm[ind] <- "No" # NA
trend_df$sig.lm <- as.factor(trend_df$sig.lm)

ptemp <- p.01[,ncol(p.01):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.01 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.01[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.01[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.01[ind] <- "No" # NA
trend_df$sig.01 <- as.factor(trend_df$sig.01)

ptemp <- p.05[,ncol(p.05):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.05 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.05[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.05[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.05[ind] <- "No" # NA
trend_df$sig.05 <- as.factor(trend_df$sig.05)

ptemp <- p.10[,ncol(p.10):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.10 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.10[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.10[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.10[ind] <- "No" # NA
trend_df$sig.10 <- as.factor(trend_df$sig.10)

ptemp <- p.50[,ncol(p.50):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.50 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.50[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.50[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.50[ind] <- "No" # NA
trend_df$sig.50 <- as.factor(trend_df$sig.50)

ptemp <- p.90[,ncol(p.90):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.90 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.90[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.90[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.90[ind] <- "No" # NA
trend_df$sig.90 <- as.factor(trend_df$sig.90)

ptemp <- p.95[,ncol(p.95):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.95 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.95[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.95[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.95[ind] <- "No" # NA
trend_df$sig.95 <- as.factor(trend_df$sig.95)

ptemp <- p.99[,ncol(p.99):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.99 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.99[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.99[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.99[ind] <- "No" # NA
trend_df$sig.99 <- as.factor(trend_df$sig.99)

save(tr.lm,tr.01,tr.05,tr.50,tr.95,tr.99,tr.10,tr.90,
     p.lm,p.01,p.05,p.50,p.95,p.99,p.10,p.90,
     trend_df,npo_deszn,lon,lat,file='npo_fit_panel.Rdata')

########## Part VI. Regional simple linear trends: North Atlantic Subtropical ########## 
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/regions_chl_deszn_trans.Rdata")
# setwd('/home/dzhai/quantile_regression')
# load('/home/dzhai/quantile_regression/regions_chl_deszn.Rdata')
month <- c(1:304)

# trends fit (no need centre)
tr.lm <- array(data=NA,dim=c(length(lon),length(lat)))
tr.01 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.05 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.10 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.50 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.90 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.95 <- array(data=NA,dim=c(length(lon),length(lat)))
tr.99 <- array(data=NA,dim=c(length(lon),length(lat)))
p.lm <- array(data=NA,dim=c(length(lon),length(lat)))
p.01 <- array(data=NA,dim=c(length(lon),length(lat)))
p.05 <- array(data=NA,dim=c(length(lon),length(lat)))
p.10 <- array(data=NA,dim=c(length(lon),length(lat)))
p.50 <- array(data=NA,dim=c(length(lon),length(lat)))
p.90 <- array(data=NA,dim=c(length(lon),length(lat)))
p.95 <- array(data=NA,dim=c(length(lon),length(lat)))
p.99 <- array(data=NA,dim=c(length(lon),length(lat)))
naNum <- array(data=0,dim=c(length(lon),length(lat)))
tick_na <- 0 
tick <- 0 

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(nao_deszn[i,j,]))
    
    if (naNum[i,j] <= 300){
      fit.lm <- lm(nao_deszn[i,j,] ~ month,na.action = na.omit)
      fit.01 <- rq(nao_deszn[i,j,] ~ month,tau=.01,na.action = na.omit)
      fit.05 <- rq(nao_deszn[i,j,] ~ month,tau=.05,na.action = na.omit)
      fit.10 <- rq(nao_deszn[i,j,] ~ month,tau=.10,na.action = na.omit)
      fit.50 <- rq(nao_deszn[i,j,] ~ month,tau=.50,na.action = na.omit)
      fit.90 <- rq(nao_deszn[i,j,] ~ month,tau=.90,na.action = na.omit)
      fit.95 <- rq(nao_deszn[i,j,] ~ month,tau=.95,na.action = na.omit)
      fit.99 <- rq(nao_deszn[i,j,] ~ month,tau=.99,na.action = na.omit)
      tr.lm[i,j] <- fit.lm$coefficients[2]
      tr.01[i,j] <- fit.01$coefficients[2]
      tr.05[i,j] <- fit.05$coefficients[2]
      tr.10[i,j] <- fit.10$coefficients[2]
      tr.50[i,j] <- fit.50$coefficients[2]
      tr.90[i,j] <- fit.90$coefficients[2]
      tr.95[i,j] <- fit.95$coefficients[2]
      tr.99[i,j] <- fit.99$coefficients[2]
      p.lm[i,j] <- summary(fit.lm)$coefficients[8]
      p.01[i,j] <- summary(fit.01,se='ker')$coefficients[8]
      p.05[i,j] <- summary(fit.05,se='ker')$coefficients[8]
      p.10[i,j] <- summary(fit.10,se='ker')$coefficients[8]
      p.50[i,j] <- summary(fit.50,se='ker')$coefficients[8]
      p.90[i,j] <- summary(fit.90,se='ker')$coefficients[8]
      p.95[i,j] <- summary(fit.95,se='ker')$coefficients[8]
      p.99[i,j] <- summary(fit.99,se='ker')$coefficients[8]
      tick <- tick + 1
      
    } else {
      tick_na <- tick_na + 1
    }
  }
}

# trend value
trendvalue <- tr.lm*1000
trend_raster <- matrix2raster(trendvalue, x = lon, y = lat, layer = 2) # no need centre
trend_df <- as.data.frame(trend_raster,xy = T)
colnames(trend_df) <-c('lon','lat','tr.lm')
rm(list=c('trend_raster'))

trendvalue <- tr.01[,ncol(tr.01):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.01 <- trendtemp

trendvalue <- tr.05[,ncol(tr.05):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.05 <- trendtemp

trendvalue <- tr.10[,ncol(tr.10):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.10 <- trendtemp

trendvalue <- tr.50[,ncol(tr.50):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.50 <- trendtemp

trendvalue <- tr.90[,ncol(tr.90):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.90 <- trendtemp

trendvalue <- tr.95[,ncol(tr.95):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.95 <- trendtemp

trendvalue <- tr.99[,ncol(tr.99):1]
trendvalue <- trendvalue*1000
trendtemp <- as.vector(aperm(trendvalue,c(1,2)))
trend_df$tr.99 <- trendtemp

# p value
ptemp <- p.lm[,ncol(p.lm):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.lm <- NA
ind <- which(pValuetemp <= 0.05)
trend_df$sig.lm[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.lm[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.lm[ind] <- "No" # NA
trend_df$sig.lm <- as.factor(trend_df$sig.lm)

ptemp <- p.01[,ncol(p.01):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.01 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.01[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.01[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.01[ind] <- "No" # NA
trend_df$sig.01 <- as.factor(trend_df$sig.01)

ptemp <- p.05[,ncol(p.05):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.05 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.05[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.05[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.05[ind] <- "No" # NA
trend_df$sig.05 <- as.factor(trend_df$sig.05)

ptemp <- p.10[,ncol(p.10):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.10 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.10[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.10[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.10[ind] <- "No" # NA
trend_df$sig.10 <- as.factor(trend_df$sig.10)

ptemp <- p.50[,ncol(p.50):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.50 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.50[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.50[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.50[ind] <- "No" # NA
trend_df$sig.50 <- as.factor(trend_df$sig.50)

ptemp <- p.90[,ncol(p.90):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.90 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.90[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.90[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.90[ind] <- "No" # NA
trend_df$sig.90 <- as.factor(trend_df$sig.90)

ptemp <- p.95[,ncol(p.95):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.95 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.95[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.95[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.95[ind] <- "No" # NA
trend_df$sig.95 <- as.factor(trend_df$sig.95)

ptemp <- p.99[,ncol(p.99):1]
pValuetemp <- as.vector(aperm(ptemp,c(1,2)))
trend_df$sig.99 <- NA
ind <- which(pValuetemp < 0.05)
trend_df$sig.99[ind] <- "Yes"
ind <- which(pValuetemp >= 0.05)
trend_df$sig.99[ind] <- "No"
ind <- which(is.na(pValuetemp))
trend_df$sig.99[ind] <- "No" # NA
trend_df$sig.99 <- as.factor(trend_df$sig.99)

save(tr.lm,tr.01,tr.05,tr.50,tr.95,tr.99,tr.10,tr.90,
     p.lm,p.01,p.05,p.50,p.95,p.99,p.10,p.90,
     trend_df,nao_deszn,lon,lat,file='nao_fit_panel.Rdata')

######################################################################





