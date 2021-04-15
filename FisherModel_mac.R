## Based on Fisher et al. 2008 model to estimate ET
library(raster)
library(randomForest)
library(solrad)
library(lubridate)
library(sirad)
## Canopy transpiration
# LEc = (1-f_wet)*f_g*f_T*f_M*alpha*(delta/(delta+gamma))*Rnc

## slope aspect
dem <- raster('/Users/NigelKinney/Box/GEE/slope_aspect_otherattributes/LiDAR/DEM30.tif')
latlon <- as.data.frame(dem, xy = TRUE)
latlon$lon <- latlon$x
latlon$lat <- latlon$y
coordinates(latlon) <- ~x + y
gridded(latlon) <- TRUE
lat <- raster(latlon, layer = 'lat')
lon <- raster(latlon, layer = 'lon')

slope <- terrain(dem, opt = 'slope', unit = 'degrees')
slope <- resample(slope, dem, method = 'ngb')
aspect <- terrain(dem, opt = 'aspect', unit = 'degrees')
aspect <- resample(aspect, dem, method = 'ngb')


## Let us do a randomforest to get SW_out and LW_out
RFpd = 1
if(RFpd){
  df <- read.csv('/Users/NigelKinney/Box/DataScience/Ameri/CoreSites/Fluxnet/Para_Import/FULLSET_Daily/FLX_US-NR1_FLUXNET2015_FULLSET_DD_1998-2014_1-3.csv')
  #master_dir <- '/Users/NigelKinney/Box/DataScience/Ameri/CoreSites/Fluxnet/Para_Import/FULLSET_Daily'
  
  Variable_names <- c('TIMESTAMP', 'SW_IN_F', 'SW_OUT', 'LW_IN_F', 'TA_F', 'VPD_F', 'LW_OUT')
  df <- df[, c(which(names(df) %in% Variable_names))]
  df$TIMESTAMP <- as.Date(as.character(df$TIMESTAMP), format = '%Y%m%d')
  df$LW_OUT[which(df$LW_OUT == -9999)] = NA
  df$SW_OUT[which(df$SW_OUT == -9999)] = NA
  df <- na.omit(df)
  df$doy <- yday(df$TIMESTAMP)
  plot(df$LW_OUT)
  lines(df$LW_IN_F, col = 2)
  
  train_idx <- sample(1:dim(df)[1], 0.7 * dim(df)[1])
  rf_train <- df[train_idx, -1]
  rf_test <- df[-train_idx, -1]
  
  RF_LWOUT <- randomForest(x = rf_train[, c('TA_F', 'LW_IN_F', 'VPD_F', 'doy')], y = rf_train$LW_OUT, importance = TRUE)
  train_pd <- predict(RF_LWOUT, rf_train[, c('TA_F', 'LW_IN_F', 'VPD_F', 'doy')])
  test_pd <- predict(RF_LWOUT, rf_test[, c('TA_F', 'LW_IN_F', 'VPD_F', 'doy')])
  
  plot(train_pd, rf_train$LW_OUT)
  plot(test_pd, rf_test$LW_OUT)
  
  RF_SWOUT <- randomForest(x = rf_train[, c('TA_F', 'LW_IN_F', 'VPD_F', 'doy')], y = rf_train$SW_OUT, importance = TRUE)
  train_pd <- predict(RF_SWOUT, rf_train[, c('TA_F',  'LW_IN_F', 'VPD_F', 'doy')])
  test_pd <- predict(RF_SWOUT, rf_test[, c('TA_F', 'LW_IN_F', 'VPD_F', 'doy')])
  
  plot(train_pd, rf_train$SW_OUT)
  plot(test_pd, rf_test$SW_OUT)
  
}
load('/Users/NigelKinney/Box/GEE/slope_aspect_otherattributes/BB.RData')

## set the extent to the pumphouse site.
doyouwantPHsite = 0
if(doyouwantPHsite){
  ex <- extent(-106.9555, -106.9434, 38.91701, 38.9274)
  dem <- crop(dem, ex)
  slope <- crop(slope, ex)
  aspect <- crop(aspect, ex)
}

## ONLY CONSIDER SMALL DOMAIN
ONLYSMALL = 1
if(ONLYSMALL == 1){
   ex <- extent(-106.9555, -106.9434, 38.91701, 38.9274)
  dem <- crop(dem, ex)
  slope <- crop(slope, ex)
  aspect <- crop(aspect, ex)
  lat <- crop(lat, ex)
  lon <- crop(lon, ex)
}

## Fisher model begins here
filelist <- list.files('/Users/NigelKinney/Box/GEE/Landsat8_Geotiff_SR/RF_Daily/NDVI/ras')
length(filelist)
ranges <- which(substr(filelist, 5, 6) == '07')
for(idx in ranges){
  VP <- raster(paste('/Users/NigelKinney/Box/GEE/DAYMET/RF_clean/vp/', filelist[idx], sep = ''))
  TMEAN <- raster(paste('/Users/NigelKinney/Box/GEE/DAYMET/RF_clean/tmean/', filelist[idx], sep = ''))
  TMAX <- raster(paste('/Users/NigelKinney/Box/GEE/DAYMET/RF_clean/tmax/', filelist[idx], sep = ''))
  TMIN <- raster(paste('/Users/NigelKinney/Box/GEE/DAYMET/RF_clean/tmin/', filelist[idx], sep = ''))
  
  crs(TMAX) <- crs(TMIN) <- crs(VP) <- crs(TMEAN) <- '+proj=utm +zone=13 +datum=WGS84 +units=m +ellps=WGS84'
  r_proj <- "+proj=longlat +datum=WGS84"
  VP <- resample(projectRaster(from = VP, crs = r_proj, method = 'ngb'), dem, method = 'ngb')
  TMEAN <- resample(projectRaster(from = TMEAN, crs = r_proj, method = 'ngb'), dem, method = 'ngb')
  TMAX <- resample(projectRaster(from = TMAX, crs = r_proj, method = 'ngb'), dem, method = 'ngb') 
  TMIN <- resample(projectRaster(from = TMIN, crs = r_proj, method = 'ngb'), dem, method = 'ngb')  

  
  SVP <- 610.7*10^(7.5*TMEAN / (237.3+TMEAN)) / 100 # saturated vapor pressure, hPa
  SVP <- 0.61094 * exp(17.625*TMEAN/(TMEAN + 243.04)) * 1000
  RH <- (VP/SVP)
  f_wet <- RH^4
  
  # f_g: green canopy fraction, f_g = fapar/fipar
  # fapar = m1SAVI + b1; fipar = m2NDVI +b2
  NDVI <- raster(paste('/Users/NigelKinney/Box/GEE/Landsat8_Geotiff_SR/RF_Daily/NDVI/ras/', filelist[idx], sep = ''))
  NDVI <- resample(NDVI, dem, method = 'ngb')
  SAVI <- raster(paste('/Users/NigelKinney/Box/GEE/Landsat8_Geotiff_SR/RF_Daily/SAVI/ras/', filelist[idx], sep = ''))
  SAVI <- resample(SAVI, dem, method = 'ngb')
  m1 = 1.2*1.136
  b1 = -0.04*1.2
  m2 = 1
  b2 = -0.05  #-0.05
  f_apar <- m1 * SAVI + b1
  f_apar[which(f_apar[]<0)] = 0
  f_ipar <- m2 * NDVI + b2
  f_g = f_apar/f_ipar
  f_g[which(f_g[]<0)] <- 0
  
  ## f_M: plant moisture constraint. f_m = f_apar/f_aparmax
  f_aparmax = 1
  f_M = f_apar/f_aparmax
  
  # alpha, delta and gamma
  alpha  = 1.26 
  # delta is the slope of saturation-to-vapor pressurre curve, http://www.fao.org/3/X0490E/x0490e07.htm#:~:text=The%20slope%20of%20the%20saturation%20vapour%20pressure%20curve%2C%20D%2C%20is,ETo%20from%20climatic%20data.&text=The%20actual%20vapour%20pressure%20(e,the%20water%20in%20the%20air.
  #delta = 4098*0.6108*exp(12.27*TMEAN/(TMEAN + 237.3))/(TMEAN + 237.3)^2
  delta = deltaVP(TMAX, TMIN)
  # gamma is the psychrometric constat = 0.066
  gamma = 0.066
  
  ## Rnc
  DirSW <- 0
  DifSW <- 0
  dateinput <- as.Date(substr(filelist[idx], 1, 8), format = '%Y%m%d')
  for(hour in 0:23){
    temp_DIR <- DirectRadiation(yday(dateinput)+hour/24, lat, lon, lon, 0, dem, slope, aspect+180)
    temp_DIF <- DiffuseRadiation(yday(dateinput)+hour/24, lat, lon, lon, 0, dem, slope)
    DirSW <- DirSW + temp_DIR/24
    DifSW <- DifSW + temp_DIF/24
    print(hour)
  }
  scale <- BB$Ratio[which(BB$Date == as.Date(substr(filelist[idx], 1, 8), format = '%Y%m%d'))]
  DirSW <- DirSW * scale
  SW <- DirSW + DifSW
  ## Try longwave radiation
  bozman = 5.67*10^(-8)
  
  #Herrero-Polo_cs ea = a + bRH + cTa, a = -1.17; b = 0.0016; c = 0.0062
  ea <- -1.17 + 0.0016 * RH + 0.0062 * (TMEAN + 273.15)
  LW <- ea * bozman * (TMEAN + 273.15) ^ 4
  
  rr <- brick(TMEAN, SW, LW, (SVP - RH)/100)
  names(rr) <- c('TA_F', 'SW_IN_F', 'LW_IN_F', 'VPD_F')
  rdf <- na.omit(as.data.frame(rr, xy = TRUE))
  rdf$doy <- yday(as.Date(substr(filelist[idx], 1, 8), format = '%Y%m%d'))
  rdf$LWout <- predict(RF_LWOUT, rdf[, c('TA_F', 'LW_IN_F', 'VPD_F', 'doy')])
  rdf$SWout <- predict(RF_SWOUT, rdf[, c('TA_F', 'LW_IN_F', 'VPD_F', 'doy')])
  
  coordinates(rdf) <- ~x+y
  gridded(rdf) <- TRUE
  LW_out <- raster(rdf, layer = 'LWout')
  SW_out <- raster(rdf, layer = 'SWout')
  SW_out <- resample(SW_out, SW)
  LW_out <- resample(LW_out, LW)
  Rn = SW - SW_out + LW - LW_out
  Rn[which(Rn[]<0)] = 0
  plot(Rn)

  # f_T: Plant temperature constraint f_T = exp(-((T_max - T_opt)/lambda)^2)
  # As we do not have PAR etc..let us just assume a T_opt = 20
  T_opt = 10 # optimum plant growth tempreature
  lambda = T_opt
  f_T = exp(-((TMEAN - T_opt)/lambda)^2)
  
  # net radiation to the soil
  k_rn = 0.6
  k_par = 0.5
  f_c = f_ipar
  LAI = -log(1-f_c)/k_par
  Rns = Rn * exp(-k_rn*LAI)
  Rnc = Rn - Rns
  
  # Canopy transpiration
  LEc = (1 - f_wet) * f_g * f_T * f_M * alpha * delta/(delta + gamma) * Rnc
  
  ## interception evaporation
  LEi = f_wet*alpha*(delta/(delta + gamma)) * Rnc
  
  ## LEs soil evaporation
  VPD = SVP - VP
  f_SM = RH^(VPD/1000)
  LEs = (f_wet + f_SM * (1 - f_wet)) * alpha * (delta/(delta + gamma)) * Rns
  plot(LEs)
  plot(LEi)
  plot(LEc)
  LE = LEs + LEi + LEc
  
  # Now let us convert LE to ET
  LE_fun <- function(T){
    L = 2500.8 - 2.36*T + 0.0016*T^2 - 0.00006*T^3
    return(L)
  }
  ET <- (60*60*24)/(LE_fun(TMEAN))/10^3 * LE
  ETc <- (60*60*24)/(LE_fun(TMEAN))/10^3 * LEc
  ETi <- (60*60*24)/(LE_fun(TMEAN))/10^3 * LEi
  ETs <- (60*60*24)/(LE_fun(TMEAN))/10^3 * LEs
  
  par(mfrow = c(1, 1))
  plot(ET, main = as.Date(substr(filelist[idx], 1, 8), format = '%Y%m%d'))
  
  figname <- paste('/Users/NigelKinney/Box/GEE/ET_FisherModel/pic/', substr(filelist[idx], 1, 8), '.png', sep = '')
  png(figname, width = 600, height = 600)
  par(mfrow = c(2, 2),     # 2x2 layout
      oma = c(2, 0, 0, 2),
       # two rows of text at the outer left and bottom margin
      mar = c(3, 3, 3, 3)) # space for one row of text at ticks and to separate plot) 
  plot(ET, main = paste('ET:', substr(filelist[idx], 1, 8)))
  plot(ETc, main = 'Transpiration')
  plot(ETi, main = 'Canopy evaporation')
  plot(ETs, main = 'Soil evaporation')
  dev.off()
  
  ETcname <- paste('/Users/NigelKinney/Box/GEE/ET_FisherModel/ETc/', filelist[idx], sep = '')
  writeRaster(ETc, ETcname, overwrite = TRUE)
  ETiname <- paste('/Users/NigelKinney/Box/GEE/ET_FisherModel/ETi/', filelist[idx], sep = '')
  writeRaster(ETi, ETiname, overwrite = TRUE)
  ETsname <- paste('/Users/NigelKinney/Box/GEE/ET_FisherModel/ETs/', filelist[idx], sep = '')
  writeRaster(ETs, ETsname, overwrite = TRUE)
  ETname <- paste('/Users/NigelKinney/Box/GEE/ET_FisherModel/ET/', filelist[idx], sep = '')
  writeRaster(ET, ETname, overwrite = TRUE)
  
  
  rET <- brick(dem, slope, aspect, Rn, SW, LW, SW_out, LW_out, NDVI, SAVI, TMEAN, TMAX, TMIN, VP, ETi, ETs, ETc, ET)
  rET <- na.omit(as.data.frame(rET, xy = TRUE))
  names(rET) <- c('x', 'y', 'elevation', 'slope', 'aspect', 'Rn', 'SW_in', 'LW_in', 'SW_out', 'LW_out',
                  'NDVI', 'SAVI', 'Tmean', 'TMAX', 'TMIN', 'VP', 'ETi', 'ETs', 'ETc', 'ET')
  rETname <- paste('/Users/NigelKinney/Box/GEE/ET_FisherModel/df/', substr(filelist[idx], 1, 8), '.csv', sep = '')
  write.csv(rET, rETname)  
}

