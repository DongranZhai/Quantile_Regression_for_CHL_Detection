# Notes: Seasonality of CHL
#        50th and 95th
#        Winter (Dec, Jan, Feb) and Summer (Jun, Jul, Aug), with averaged 3-month
#        raw data without deseasonality
#        on the global scale
library(dplyr)
library(R.matlab)
library(raster)
library(oceanmap)
library(orcutt)
library(quantreg)
########## Part I: Extract regions: high latitude in NP EP, and NA ##########
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
# monchl

temp <- readMat("/Volumes/Doris/Projects/sptmodel/Longhurst_180.mat")
Longhurst <- temp$Longhurst
nlon <- length(lon)
nlat <- length(lat)
area <- sort(unique(Longhurst[!is.na(Longhurst)]))#list of longhurst areas to iterate through
area <- area[-c(1,2,4,6,7,8,10,11,12,13,17,18,19,22,25,26,27,29,30,31,34,35,36,40,43,41,42,47,50,53,54)] #remove coastal, polar, GoM, Archipelagic Deep basin, Mediterranean

##### North Pacific (Longhurst is 18 & 19, which corresponding to 24 & 46 in area)
long_np <- array(data=NA,dim=c(nlon,nlat))
npchl <- array(data=NA,dim=c(nlon,nlat,298))

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

for(i in c(1:298)){
  tempchl <- monchl[,,i]
  tempchl[ind] <- NA
  npchl[,,i] <- tempchl
}
##### North Atlantic (Longhurst is 13, which corresponding to 39 in area)
long_na <- array(data=NA,dim=c(nlon,nlat))
nachl <- array(data=NA,dim=c(nlon,nlat,298))

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

for(i in c(1:298)){
  tempchl <- monchl[,,i]
  tempchl[ind] <- NA
  nachl[,,i] <- tempchl
}

##### Equatorial Pacific (Longhurst is 5 & 7, which corresponding to 28 & 14 in area)
long_ep <- array(data=NA,dim=c(nlon,nlat))
epchl <- array(data=NA,dim=c(nlon,nlat,298))

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

for(i in c(1:298)){
  tempchl <- monchl[,,i]
  tempchl[ind] <- NA
  epchl[,,i] <- tempchl
}

save(npchl,epchl,nachl,lon,lat,file='seasonality_chl_regions.Rdata')

########## Part II: dataframe generalization ##########
# np
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/seasonality_chl_regions.Rdata")
month <- c(1:298)
# center Pacific
np_centre_df<- data.frame(matrix(ncol = 4, nrow = 1)) # array(matrix) to raster
colnames(np_centre_df) <- c('x','y','npchl','layer')
for(i in month){
  print(i)
  temp_r <- matrix2raster(npchl[,,i], x = lon, y = lat, layer = i)
  
  temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
  temp2 <- crop(temp_r, extent(0, 180, -90, 90))
  extent(temp1) <- c(180, 360, -90, 90)
  
  temp <- merge(temp1, temp2)
  temp <- as.data.frame(temp,xy=T)
  temp$npchl <- i
  np_centre_df <- rbind(np_centre_df,temp) # np_centre is a df with Pacific as centre
}
colnames(np_centre_df) <- c('x','y','layer','npchl')
np_centre_df <- np_centre_df[-1,]
np_centre <- array(data=NA,dim=c(360,180,298))

for(i in month){
  print(i)
  temp <- np_centre_df[which(np_centre_df$layer== i),c(1:2,4)]
  coordinates(temp) <- ~ x + y
  gridded(temp) <- T
  temp_r <- raster(temp)
  np_centre[,,i] <- raster2matrix(temp_r)
}
# np_centre_df and np_centre

# ep
# center Pacific
ep_centre_df<- data.frame(matrix(ncol = 4, nrow = 1)) # array(matrix) to raster
colnames(ep_centre_df) <- c('x','y','epchl','layer')

for(i in month){
  print(i)
  temp_r <- matrix2raster(epchl[,,i], x = lon, y = lat, layer = i)
  
  temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
  temp2 <- crop(temp_r, extent(0, 180, -90, 90))
  extent(temp1) <- c(180, 360, -90, 90)
  
  temp <- merge(temp1, temp2)
  temp <- as.data.frame(temp,xy=T)
  temp$epchl <- i
  ep_centre_df <- rbind(ep_centre_df,temp) # ep_centre is a df with Pacific as centre
}
colnames(ep_centre_df) <- c('x','y','layer','epchl')
ep_centre_df <- ep_centre_df[-1,]
ep_centre <- array(data=NA,dim=c(360,180,298))

for(i in month){
  print(i)
  temp <- ep_centre_df[which(ep_centre_df$layer== i),c(1:2,4)]
  coordinates(temp) <- ~ x + y
  gridded(temp) <- T
  temp_r <- raster(temp)
  ep_centre[,,i] <- raster2matrix(temp_r)
}

# NA
na_df<- data.frame(matrix(ncol = 4, nrow = 1)) # array(matrix) to raster
colnames(na_df) <- c('x','y','nachl','layer')

for(i in month){
  print(i)
  temp_r <- matrix2raster(nachl[,,i], x = lon, y = lat, layer = i)

  temp_r <- as.data.frame(temp_r,xy=T)
  temp_r$nachl <- i
  na_df <- rbind(na_df,temp_r)
}
colnames(na_df) <- c('x','y','layer','nachl')
na_df <- na_df[-1,]

save(np_centre,np_centre_df,ep_centre,ep_centre_df,na_df,
     npchl,epchl,nachl,lon,lat,file='seasonality_chl_regions.Rdata')

########## Part III: Winter ts ##########
# NP Winter
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/seasonality_chl_regions.Rdata")
year <- 1998:2022
nyear <- length(year)
np_winter_ts <- array(data=NA,dim=c(length(lon),length(lat),nyear))

indl <- 4
for(i in 1:nyear){
  print(i)
  indr <- indl+2
  wintertemp <- npchl[,,c(indl:indr)]
  print(paste0('chl.dim: ',indl,'--',indr))
  wintertemp <- apply(wintertemp,c(1,2),mean)
  np_winter_ts[,,i] <- wintertemp
  indl <- indr+10
}

# NA Winter
na_winter_ts <- array(data=NA,dim=c(length(lon),length(lat),nyear))

indl <- 4
for(i in 1:nyear){
  print(i)
  indr <- indl+2
  wintertemp <- nachl[,,c(indl:indr)]
  print(paste0('chl.dim: ',indl,'--',indr))
  wintertemp <- apply(wintertemp,c(1,2),mean)
  na_winter_ts[,,i] <- wintertemp
  indl <- indr+10
}

# EP Winter
ep_winter_ts <- array(data=NA,dim=c(length(lon),length(lat),nyear))

indl <- 4
for(i in 1:nyear){
  print(i)
  indr <- indl+2
  wintertemp <- epchl[,,c(indl:indr)]
  print(paste0('chl.dim: ',indl,'--',indr))
  wintertemp <- apply(wintertemp,c(1,2),mean)
  ep_winter_ts[,,i] <- wintertemp
  indl <- indr+10
}

save(np_winter_ts,na_winter_ts,ep_winter_ts,lon,lat,file='winter_ts.Rdata')

########## Part IV: Winter 50 and 95 in NP, NA, and EP ##########
########## NP ##########
rm(list = ls())
load("/Volumes/Doris/OCCCI/quantile_regression/winter_ts.Rdata")
year <- 1998:2020
np_tr_50 <- array(data=NA,dim=c(length(lon),length(lat)))
np_pValue_50 <- array(data=NA,dim=c(length(lon),length(lat)))
np_tr_95 <- array(data=NA,dim=c(length(lon),length(lat)))
np_pValue_95 <- array(data=NA,dim=c(length(lon),length(lat)))
naNum<- array(data=0,dim=c(length(lon),length(lat)))

tick_na <- 0
tick <- 0
z1=0.50
z2=0.95

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(np_winter_ts[i,j,]))
    
    if (naNum[i,j] <= 22){
      fit <- rq(np_winter_ts[i,j,] ~ c(1:25),tau=z1,na.action = na.omit)
      np_tr_50[i,j] <- fit$coefficients[[2]]
      np_pValue_50[i,j] <- summary(fit,se='ker')$coefficients[8]#[14]
      
      fit <- rq(np_winter_ts[i,j,] ~ c(1:25),tau=z2,na.action = na.omit)
      np_tr_95[i,j] <- fit$coefficients[[2]]
      np_pValue_95[i,j] <- summary(fit,se='ker')$coefficients[8]#[14]
      tick <- tick + 1
      
    } else {
      tick_na <- tick_na + 1 
    }
  }
}

# center Pacific 50
temp_r <- matrix2raster(np_tr_50[,], x = lon, y = lat) #trend
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
centre_df <- merge(temp1, temp2)
centre_df <- as.data.frame(centre_df,xy=T)
temp_r <- matrix2raster(np_pValue_50[,], x = lon, y = lat) #pvalue
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
temp <- merge(temp1, temp2)
temp <- as.data.frame(temp,xy=T)
centre_df$pValue <- temp[,3]
colnames(centre_df) <- c('lon','lat','trend','pvalue')
np_centre_df_50 <- centre_df

# center Pacific 95
temp_r <- matrix2raster(np_tr_95[,], x = lon, y = lat) #trend
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
centre_df <- merge(temp1, temp2)
centre_df <- as.data.frame(centre_df,xy=T)
temp_r <- matrix2raster(np_pValue_95[,], x = lon, y = lat) #pvalue
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
temp <- merge(temp1, temp2)
temp <- as.data.frame(temp,xy=T)
centre_df$pValue <- temp[,3]
colnames(centre_df) <- c('lon','lat','trend','pvalue')
np_centre_df_95 <- centre_df

########## NA ##########
na_tr_50 <- array(data=NA,dim=c(length(lon),length(lat)))
na_pValue_50 <- array(data=NA,dim=c(length(lon),length(lat)))
na_tr_95 <- array(data=NA,dim=c(length(lon),length(lat)))
na_pValue_95 <- array(data=NA,dim=c(length(lon),length(lat)))
naNum<- array(data=0,dim=c(length(lon),length(lat)))

tick_na <- 0
tick <- 0
z1=0.50
z2=0.95

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(na_winter_ts[i,j,]))
    
    if (naNum[i,j] <= 22){
      fit <- rq(na_winter_ts[i,j,] ~ c(1:25),tau=z1,na.action = na.omit)
      na_tr_50[i,j] <- fit$coefficients[[2]]
      na_pValue_50[i,j] <- summary(fit,se='ker')$coefficients[8]#[14]
      
      fit <- rq(na_winter_ts[i,j,] ~ c(1:25),tau=z2,na.action = na.omit)
      na_tr_95[i,j] <- fit$coefficients[[2]]
      na_pValue_95[i,j] <- summary(fit,se='ker')$coefficients[8]#[14]
      tick <- tick + 1
      
    } else {
      tick_na <- tick_na + 1 
    }
  }
}

# center Pacific 50
temp_r <- matrix2raster(na_tr_50[,], x = lon, y = lat) #trend
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
centre_df <- merge(temp1, temp2)
centre_df <- as.data.frame(centre_df,xy=T)
temp_r <- matrix2raster(na_pValue_50[,], x = lon, y = lat) #pvalue
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
temp <- merge(temp1, temp2)
temp <- as.data.frame(temp,xy=T)
centre_df$pValue <- temp[,3]
colnames(centre_df) <- c('lon','lat','trend','pvalue')
na_centre_df_50 <- centre_df

# center Pacific 95
temp_r <- matrix2raster(na_tr_95[,], x = lon, y = lat) #trend
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
centre_df <- merge(temp1, temp2)
centre_df <- as.data.frame(centre_df,xy=T)
temp_r <- matrix2raster(na_pValue_95[,], x = lon, y = lat) #pvalue
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
temp <- merge(temp1, temp2)
temp <- as.data.frame(temp,xy=T)
centre_df$pValue <- temp[,3]
colnames(centre_df) <- c('lon','lat','trend','pvalue')
na_centre_df_95 <- centre_df

########## EP ##########
ep_tr_50 <- array(data=NA,dim=c(length(lon),length(lat)))
ep_pValue_50 <- array(data=NA,dim=c(length(lon),length(lat)))
ep_tr_95 <- array(data=NA,dim=c(length(lon),length(lat)))
ep_pValue_95 <- array(data=NA,dim=c(length(lon),length(lat)))
naNum<- array(data=0,dim=c(length(lon),length(lat)))

tick_na <- 0
tick <- 0
z1=0.50
z2=0.95

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(ep_winter_ts[i,j,]))
    
    if (naNum[i,j] <= 22){
      fit <- rq(ep_winter_ts[i,j,] ~ c(1:25),tau=z1,na.action = na.omit)
      ep_tr_50[i,j] <- fit$coefficients[[2]]
      ep_pValue_50[i,j] <- summary(fit,se='ker')$coefficients[8]#[14]
      
      fit <- rq(ep_winter_ts[i,j,] ~ c(1:25),tau=z2,na.action = na.omit)
      ep_tr_95[i,j] <- fit$coefficients[[2]]
      ep_pValue_95[i,j] <- summary(fit,se='ker')$coefficients[8]#[14]
      tick <- tick + 1
      
    } else {
      tick_na <- tick_na + 1 
    }
  }
}

# center Pacific 50
temp_r <- matrix2raster(ep_tr_50[,], x = lon, y = lat) #trend
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
centre_df <- merge(temp1, temp2)
centre_df <- as.data.frame(centre_df,xy=T)
temp_r <- matrix2raster(ep_pValue_50[,], x = lon, y = lat) #pvalue
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
temp <- merge(temp1, temp2)
temp <- as.data.frame(temp,xy=T)
centre_df$pValue <- temp[,3]
colnames(centre_df) <- c('lon','lat','trend','pvalue')
ep_centre_df_50 <- centre_df

# center Pacific 95
temp_r <- matrix2raster(ep_tr_95[,], x = lon, y = lat) #trend
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
centre_df <- merge(temp1, temp2)
centre_df <- as.data.frame(centre_df,xy=T)
temp_r <- matrix2raster(ep_pValue_95[,], x = lon, y = lat) #pvalue
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
temp <- merge(temp1, temp2)
temp <- as.data.frame(temp,xy=T)
centre_df$pValue <- temp[,3]
colnames(centre_df) <- c('lon','lat','trend','pvalue')
ep_centre_df_95 <- centre_df

savename <- paste0("winter_regions_quant.Rdata")
save(np_tr_50,np_pValue_50,np_centre_df_50,np_tr_95,np_pValue_95,np_centre_df_95,
     na_tr_50,na_pValue_50,na_centre_df_50,na_tr_95,na_pValue_95,na_centre_df_95,
     ep_tr_50,ep_pValue_50,ep_centre_df_50,ep_tr_95,ep_pValue_95,ep_centre_df_95,
     lon,lat,file = savename)


########## Part V: Summer ts ##########
# NP summer
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/quantile_regression/seasonality_chl_regions.Rdata")
year <- 1998:2021
nyear <- length(year)
np_summer_ts <- array(data=NA,dim=c(length(lon),length(lat),nyear))

indl <- 10
for(i in 1:nyear){
  print(i)
  indr <- indl+2
  summertemp <- npchl[,,c(indl:indr)]
  print(paste0('chl.dim: ',indl,'--',indr))
  summertemp <- apply(summertemp,c(1,2),mean)
  np_summer_ts[,,i] <- summertemp
  indl <- indr+10
}

# NA summer
na_summer_ts <- array(data=NA,dim=c(length(lon),length(lat),nyear))

indl <- 10
for(i in 1:nyear){
  print(i)
  indr <- indl+2
  summertemp <- nachl[,,c(indl:indr)]
  print(paste0('chl.dim: ',indl,'--',indr))
  summertemp <- apply(summertemp,c(1,2),mean)
  na_summer_ts[,,i] <- summertemp
  indl <- indr+10
}

# EP summer
ep_summer_ts <- array(data=NA,dim=c(length(lon),length(lat),nyear))

indl <- 10
for(i in 1:nyear){
  print(i)
  indr <- indl+2
  summertemp <- epchl[,,c(indl:indr)]
  print(paste0('chl.dim: ',indl,'--',indr))
  summertemp <- apply(summertemp,c(1,2),mean)
  ep_summer_ts[,,i] <- summertemp
  indl <- indr+10
}

save(np_summer_ts,na_summer_ts,ep_summer_ts,lon,lat,file='summer_ts.Rdata')

########## Part VI: Summer 50 and 95 in NP, NA, and EP ##########
########## NP ##########
rm(list = ls())
load("/Volumes/Doris/OCCCI/quantile_regression/summer_ts.Rdata")
year <- 1998:2020
np_tr_50 <- array(data=NA,dim=c(length(lon),length(lat)))
np_pValue_50 <- array(data=NA,dim=c(length(lon),length(lat)))
np_tr_95 <- array(data=NA,dim=c(length(lon),length(lat)))
np_pValue_95 <- array(data=NA,dim=c(length(lon),length(lat)))
naNum<- array(data=0,dim=c(length(lon),length(lat)))

tick_na <- 0
tick <- 0
z1=0.50
z2=0.95

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(np_summer_ts[i,j,]))
    
    if (naNum[i,j] <= 22){
      fit <- rq(np_summer_ts[i,j,] ~ c(1:24),tau=z1,na.action = na.omit)
      np_tr_50[i,j] <- fit$coefficients[[2]]
      np_pValue_50[i,j] <- summary(fit,se='ker')$coefficients[8]#[14]
      
      fit <- rq(np_summer_ts[i,j,] ~ c(1:24),tau=z2,na.action = na.omit)
      np_tr_95[i,j] <- fit$coefficients[[2]]
      np_pValue_95[i,j] <- summary(fit,se='ker')$coefficients[8]#[14]
      tick <- tick + 1
      
    } else {
      tick_na <- tick_na + 1 
    }
  }
}

# center Pacific 50
temp_r <- matrix2raster(np_tr_50[,], x = lon, y = lat) #trend
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
centre_df <- merge(temp1, temp2)
centre_df <- as.data.frame(centre_df,xy=T)
temp_r <- matrix2raster(np_pValue_50[,], x = lon, y = lat) #pvalue
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
temp <- merge(temp1, temp2)
temp <- as.data.frame(temp,xy=T)
centre_df$pValue <- temp[,3]
colnames(centre_df) <- c('lon','lat','trend','pvalue')
np_centre_df_50 <- centre_df

# center Pacific 95
temp_r <- matrix2raster(np_tr_95[,], x = lon, y = lat) #trend
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
centre_df <- merge(temp1, temp2)
centre_df <- as.data.frame(centre_df,xy=T)
temp_r <- matrix2raster(np_pValue_95[,], x = lon, y = lat) #pvalue
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
temp <- merge(temp1, temp2)
temp <- as.data.frame(temp,xy=T)
centre_df$pValue <- temp[,3]
colnames(centre_df) <- c('lon','lat','trend','pvalue')
np_centre_df_95 <- centre_df

########## NA ##########
na_tr_50 <- array(data=NA,dim=c(length(lon),length(lat)))
na_pValue_50 <- array(data=NA,dim=c(length(lon),length(lat)))
na_tr_95 <- array(data=NA,dim=c(length(lon),length(lat)))
na_pValue_95 <- array(data=NA,dim=c(length(lon),length(lat)))
naNum<- array(data=0,dim=c(length(lon),length(lat)))

tick_na <- 0
tick <- 0
z1=0.50
z2=0.95

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(na_summer_ts[i,j,]))
    
    if (naNum[i,j] <= 22){
      fit <- rq(na_summer_ts[i,j,] ~ c(1:24),tau=z1,na.action = na.omit)
      na_tr_50[i,j] <- fit$coefficients[[2]]
      na_pValue_50[i,j] <- summary(fit,se='ker')$coefficients[8]#[14]
      
      fit <- rq(na_summer_ts[i,j,] ~ c(1:24),tau=z2,na.action = na.omit)
      na_tr_95[i,j] <- fit$coefficients[[2]]
      na_pValue_95[i,j] <- summary(fit,se='ker')$coefficients[8]#[14]
      tick <- tick + 1
      
    } else {
      tick_na <- tick_na + 1 
    }
  }
}

# center Pacific 50
temp_r <- matrix2raster(na_tr_50[,], x = lon, y = lat) #trend
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
centre_df <- merge(temp1, temp2)
centre_df <- as.data.frame(centre_df,xy=T)
temp_r <- matrix2raster(na_pValue_50[,], x = lon, y = lat) #pvalue
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
temp <- merge(temp1, temp2)
temp <- as.data.frame(temp,xy=T)
centre_df$pValue <- temp[,3]
colnames(centre_df) <- c('lon','lat','trend','pvalue')
na_centre_df_50 <- centre_df

# center Pacific 95
temp_r <- matrix2raster(na_tr_95[,], x = lon, y = lat) #trend
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
centre_df <- merge(temp1, temp2)
centre_df <- as.data.frame(centre_df,xy=T)
temp_r <- matrix2raster(na_pValue_95[,], x = lon, y = lat) #pvalue
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
temp <- merge(temp1, temp2)
temp <- as.data.frame(temp,xy=T)
centre_df$pValue <- temp[,3]
colnames(centre_df) <- c('lon','lat','trend','pvalue')
na_centre_df_95 <- centre_df

########## EP ##########
ep_tr_50 <- array(data=NA,dim=c(length(lon),length(lat)))
ep_pValue_50 <- array(data=NA,dim=c(length(lon),length(lat)))
ep_tr_95 <- array(data=NA,dim=c(length(lon),length(lat)))
ep_pValue_95 <- array(data=NA,dim=c(length(lon),length(lat)))
naNum<- array(data=0,dim=c(length(lon),length(lat)))

tick_na <- 0
tick <- 0
z1=0.50
z2=0.95

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(ep_summer_ts[i,j,]))
    
    if (naNum[i,j] <= 22){
      fit <- rq(ep_summer_ts[i,j,] ~ c(1:24),tau=z1,na.action = na.omit)
      ep_tr_50[i,j] <- fit$coefficients[[2]]
      ep_pValue_50[i,j] <- summary(fit,se='ker')$coefficients[8]#[14]
      
      fit <- rq(ep_summer_ts[i,j,] ~ c(1:24),tau=z2,na.action = na.omit)
      ep_tr_95[i,j] <- fit$coefficients[[2]]
      ep_pValue_95[i,j] <- summary(fit,se='ker')$coefficients[8]#[14]
      tick <- tick + 1
      
    } else {
      tick_na <- tick_na + 1 
    }
  }
}

# center Pacific 50
temp_r <- matrix2raster(ep_tr_50[,], x = lon, y = lat) #trend
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
centre_df <- merge(temp1, temp2)
centre_df <- as.data.frame(centre_df,xy=T)
temp_r <- matrix2raster(ep_pValue_50[,], x = lon, y = lat) #pvalue
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
temp <- merge(temp1, temp2)
temp <- as.data.frame(temp,xy=T)
centre_df$pValue <- temp[,3]
colnames(centre_df) <- c('lon','lat','trend','pvalue')
ep_centre_df_50 <- centre_df

# center Pacific 95
temp_r <- matrix2raster(ep_tr_95[,], x = lon, y = lat) #trend
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
centre_df <- merge(temp1, temp2)
centre_df <- as.data.frame(centre_df,xy=T)
temp_r <- matrix2raster(ep_pValue_95[,], x = lon, y = lat) #pvalue
temp1 <- crop(temp_r, extent(-180, 0, -90, 90))
temp2 <- crop(temp_r, extent(0, 180, -90, 90))
extent(temp1) <- c(180, 360, -90, 90)
temp <- merge(temp1, temp2)
temp <- as.data.frame(temp,xy=T)
centre_df$pValue <- temp[,3]
colnames(centre_df) <- c('lon','lat','trend','pvalue')
ep_centre_df_95 <- centre_df

savename <- paste0("summer_regions_quant.Rdata")
save(np_tr_50,np_pValue_50,np_centre_df_50,np_tr_95,np_pValue_95,np_centre_df_95,
     na_tr_50,na_pValue_50,na_centre_df_50,na_tr_95,na_pValue_95,na_centre_df_95,
     ep_tr_50,ep_pValue_50,ep_centre_df_50,ep_tr_95,ep_pValue_95,ep_centre_df_95,
     lon,lat,file = savename)




########## Part VII: Global seasonal ts and 50th and 95th quantile ##########
# winter
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
year <- 1998:2022
nyear <- length(year)
winter_ts <- array(data=NA,dim=c(length(lon),length(lat),nyear))

indl <- 4
for(i in 1:nyear){
  print(i)
  indr <- indl+2
  wintertemp <- monchl[,,c(indl:indr)]
  print(paste0('chl.dim: ',indl,'--',indr))
  wintertemp <- apply(wintertemp,c(1,2),mean)
  winter_ts[,,i] <- wintertemp
  indl <- indr+10
}

tr_50 <- array(data=NA,dim=c(length(lon),length(lat)))
pValue_50 <- array(data=NA,dim=c(length(lon),length(lat)))
tr_95 <- array(data=NA,dim=c(length(lon),length(lat)))
pValue_95 <- array(data=NA,dim=c(length(lon),length(lat)))
naNum<- array(data=0,dim=c(length(lon),length(lat)))

tick_na <- 0
tick <- 0
z1=0.50
z2=0.95

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(winter_ts[i,j,]))
    
    if (naNum[i,j] <= 22){
      fit <- rq(winter_ts[i,j,] ~ c(1:25),tau=z1,na.action = na.omit)
      tr_50[i,j] <- fit$coefficients[[2]]
       pValue_50[i,j] <- summary(fit,se='ker')$coefficients[8]#[14]
      
      fit <- rq(winter_ts[i,j,] ~ c(1:25),tau=z2,na.action = na.omit)
      tr_95[i,j] <- fit$coefficients[[2]]
      pValue_95[i,j] <- summary(fit,se='ker')$coefficients[8]#[14]
      tick <- tick + 1
      
    } else {
      tick_na <- tick_na + 1 
    }
  }
}

# center Pacific 50
temp_r <- matrix2raster(tr_50[,], x = lon, y = lat) #trend
centre_df <- as.data.frame(temp_r,xy=T)
temp_r <- matrix2raster(pValue_50[,], x = lon, y = lat) #pvalue
temp_r <- as.data.frame(temp_r,xy=T)
centre_df$pValue <- temp_r[,3]
colnames(centre_df) <- c('lon','lat','trend','pvalue')
centre_df_50 <- centre_df

# center Pacific 95
temp_r <- matrix2raster(tr_95[,], x = lon, y = lat) #trend
centre_df <- as.data.frame(temp_r,xy=T)
temp_r <- matrix2raster(pValue_95[,], x = lon, y = lat) #pvalue
temp_r <- as.data.frame(temp_r,xy=T)
centre_df$pValue <- temp_r[,3]
colnames(centre_df) <- c('lon','lat','trend','pvalue')
centre_df_95 <- centre_df

save(winter_ts,tr_50,pValue_50,centre_df_50,tr_95,pValue_95,centre_df_95,
     lon,lat,file='winter_ts_df.Rdata')

# summer
rm(list=ls())
setwd("/Volumes/Doris/OCCCI/quantile_regression")
load("/Volumes/Doris/OCCCI/monthly_chl_rm_coastal_occci.Rdata")
year <- 1998:2021
nyear <- length(year)
summer_ts <- array(data=NA,dim=c(length(lon),length(lat),nyear))

indl <- 10
for(i in 1:nyear){
  print(i)
  indr <- indl+2
  summertemp <- monchl[,,c(indl:indr)]
  print(paste0('chl.dim: ',indl,'--',indr))
  summertemp <- apply(summertemp,c(1,2),mean)
  summer_ts[,,i] <- summertemp
  indl <- indr+10
}

tr_50 <- array(data=NA,dim=c(length(lon),length(lat)))
pValue_50 <- array(data=NA,dim=c(length(lon),length(lat)))
tr_95 <- array(data=NA,dim=c(length(lon),length(lat)))
pValue_95 <- array(data=NA,dim=c(length(lon),length(lat)))
naNum<- array(data=0,dim=c(length(lon),length(lat)))

tick_na <- 0
tick <- 0
z1=0.50
z2=0.95

for(i in 1:length(lon)){
  for(j in 1:length(lat)){
    naNum[i,j] <- sum(is.na(summer_ts[i,j,]))
    
    if (naNum[i,j] <= 21){
      fit <- rq(summer_ts[i,j,] ~ c(1:24),tau=z1,na.action = na.omit)
      tr_50[i,j] <- fit$coefficients[[2]]
      pValue_50[i,j] <- summary(fit,se='ker')$coefficients[8]#[14]
      
      fit <- rq(summer_ts[i,j,] ~ c(1:24),tau=z2,na.action = na.omit)
      tr_95[i,j] <- fit$coefficients[[2]]
      pValue_95[i,j] <- summary(fit,se='ker')$coefficients[8]#[14]
      tick <- tick + 1
      
    } else {
      tick_na <- tick_na + 1 
    }
  }
}

# center Pacific 50
temp_r <- matrix2raster(tr_50[,], x = lon, y = lat) #trend
centre_df <- as.data.frame(temp_r,xy=T)
temp_r <- matrix2raster(pValue_50[,], x = lon, y = lat) #pvalue
temp_r <- as.data.frame(temp_r,xy=T)
centre_df$pValue <- temp_r[,3]
colnames(centre_df) <- c('lon','lat','trend','pvalue')
centre_df_50 <- centre_df

# center Pacific 95
temp_r <- matrix2raster(tr_95[,], x = lon, y = lat) #trend
centre_df <- as.data.frame(temp_r,xy=T)
temp_r <- matrix2raster(pValue_95[,], x = lon, y = lat) #pvalue
temp_r <- as.data.frame(temp_r,xy=T)
centre_df$pValue <- temp_r[,3]
colnames(centre_df) <- c('lon','lat','trend','pvalue')
centre_df_95 <- centre_df

save(summer_ts,tr_50,pValue_50,centre_df_50,tr_95,pValue_95,centre_df_95,
     lon,lat,file='summer_ts_df.Rdata')

