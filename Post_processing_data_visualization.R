## Large domain small domain comparison
## Fisher results post analysis
rm(list = ls())
library(ggplot2)
library(lubridate)
library(wesanderson)
library(raster)
library(TTR)
pal <- wes_palette("Zissou1", 100, type = "continuous")
catergroups = function(df){
  df$slopegroup <- '<10'
  df$slopegroup[which(df$slope >10 & df$slope <= 20)] = '10-20'
  df$slopegroup[which(df$slope >20 & df$slope <= 30)] = '20-30'
  df$slopegroup[which(df$slope >30 )] = '>30'
  df$slopegroup <- factor(df$slopegroup, levels = c('<10', '10-20', '20-30', '>30'))
  
  df$aspectgroup <- 'N'
  df$aspectgroup[which(df$aspect > 22.5 & df$aspect <= 67.5)] <- 'NE'
  df$aspectgroup[which(df$aspect > 67.5 & df$aspect <= 112.5)] <- 'E'
  df$aspectgroup[which(df$aspect > 112.5 & df$aspect <= 157.5)] <- 'SE'
  df$aspectgroup[which(df$aspect > 157.5 & df$aspect <= 202.5)] <- 'S'
  df$aspectgroup[which(df$aspect > 202.5 & df$aspect <= 247.5)] <- 'SW'
  df$aspectgroup[which(df$aspect > 247.5 & df$aspect <= 292.5)] <- 'W'
  df$aspectgroup[which(df$aspect > 292.5 & df$aspect <= 337.5)] <- 'NW'
  df$aspectgroup <- factor(df$aspectgroup, levels = c('N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'))
  
  df$elevgroup <- 'Low: < 3000'
  df$elevgroup[which(df$elevation > 3000 & df$elevation <= 3500)] <- 'Medium: 3000 - 3500'
  df$elevgroup[which(df$elevation > 3500)] <- 'High: >3500'
  df$elevgroup <- factor(df$elevgroup, levels = c('Low: < 3000', 'Medium: 3000 - 3500', 'High: > 3500'))
  return(df)
}

## Snow
Snow_18 <- raster('C:/Users/Nigel/Box/GEE/slope_aspect_otherattributes/SWE/20180331.tif')
Snow_18 <- addLayer(Snow_18, raster('C:/Users/Nigel/Box/GEE/slope_aspect_otherattributes/SWE/20180524.tif')) 
Snow_19 <- raster('C:/Users/Nigel/Box/GEE/slope_aspect_otherattributes/SWE/20190407.tif')
Snow_19 <- addLayer(Snow_19, raster('C:/Users/Nigel/Box/GEE/slope_aspect_otherattributes/SWE/20190610.tif')) 

ex <- extent(-107, -106.94, 38.92, 38.98)
#ex <- extent(-106.9555, -106.9434, 38.91701, 38.9274)
sm_18 <- crop(Snow_18, ex)
sm_19 <- crop(Snow_19, ex)
plot(sm_18$X20180331, xlab = 'lon', ylab = 'lat', main = '2018-03-31 SWE')
plot(sm_19$X20190407, xlab = 'lon', ylab = 'lat', main = '2019-04-07 SWE')
plot(sm_18$X20180524, xlab = 'lon', ylab = 'lat', main = '2018-05-24 SWE')
plot(sm_19$X20190610)

sm_18 <- as.data.frame(sm_18, xy = TRUE)
names(sm_18)[3:4] <- c('SWE1', 'SWE2')
sm_19 <- as.data.frame(sm_19, xy = TRUE)
names(sm_19)[3:4] <- c('SWE1', 'SWE2')

## ER
ERWR <- data.frame(lon = c(-106.9511161, -106.9855073, -106.9835507, -106.9911594),
                   lat = c(38.92230895, 38.92949275, 38.96065217, 38.96021739))

setwd('C:/Users/Nigel/Box/GEE/ET_FisherModel/True_RS_albedo_correct')
files <- list.files()
pos <- read.csv('C:/Users/Nigel/Box/GEE/Probe_SurfaceTemp/sensor_loc.csv')
pos <- pos[which(pos$Status == 'OK'),]
load('LS_df.RData')
load('S2_df.RData')
LS_df$sate <- 'L8'
S2_df$sate <- 'S2'
df <- rbind(LS_df, S2_df)
rm(list = c('LS_df', 'S2_df'))
df <- catergroups(df)
#df_small <- df[which(df$x >= -106.9555 & df$x <= -106.9434 & df$y >= 38.91701 & df$y <= 38.9274),]
df_small <- df
df_small$compgr <- paste(df_small$aspectgroup, ' + ', df_small$slopegroup, sep = '')
df_small$site <- 'Watershed'
df_small$site[which(df_small$x >= -106.9555 & df_small$x <= -106.9434 & df_small$y >= 38.91701 & df_small$y <= 38.9274)] = 'Pumphouse'

## Need to calculate some constraints following Fisher 2008
# plant moisture constraint
#df_small$SVP <- 610.7*10^(7.5*df_small$Tmean / (237.3+df_small$Tmean)) / 100 # saturated vapor pressure, hPa
df_small$SVP <- 0.61094 * exp(17.625*df_small$Tmean/(df_small$Tmean + 243.04)) * 1000
df_small$RH <- (df_small$VP/df_small$SVP)
df_small$f_SM <- df_small$RH ^((df_small$SVP - df_small$VP)/1000)

T_opt = 15 # optimum plant growth tempreature
lambda = T_opt
#plant temperature constraint
df_small$f_T <- exp(-((df_small$Tmean - T_opt)/lambda)^2)
sm <- df_small[which(df_small$site == 'Pumphouse'),]
##
load('C:/Users/Nigel/Box/GEE/ET_FisherModel/NR.RData')
#baddays <- as.Date(c("2017-06-15", "2018-04-22", "2018-05-10", "2018-06-18", "2018-07-04"))
u <- aggregate(df_small, by = list(df_small$Date), var)
u$Group.1[which(u$ET > 0.5)] # this returns significantly large variacne days
baddays <- as.Date(c("2017-06-15", "2018-04-22", "2018-05-10", "2018-06-18", "2018-07-04"))
'%ni%' <- Negate('%in%')
df_small <- df_small[which(df_small$Date %ni% baddays),]
library(dplyr)
df_stat <- df_small
m <- df_stat%>% 
  group_by(Date) %>% 
  summarize(mET = mean(ET),
            sdET = sd(ET),
            q75 = quantile(ET, 0.75),
            q25 = quantile(ET, 0.25))

ggplot() + geom_point(data = NR, aes(x = Date, y = ET, color = 'Niwot Ridge')) +
  geom_point(data = m, aes(x = Date, y = mET, color = 'East River')) + 
  geom_errorbar(data = m, aes(x = Date, ymin = q25, ymax = q75)) + geom_line() +
  ggtitle('East River versus Niwot Ridge, ET') + xlab('Date') + ylab('ET') + ylim(0, 4) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

mm <- aggregate(df_small[, c('Date', 'ET')], by = list(df_small$Date), mean)
mm <- aggregate(mm, by = list(month(mm$Date), year(mm$Date)), mean)
nn <- aggregate(NR[, c('Date', 'ET')], by = list(month(NR$Date), year(NR$Date)), mean)
cc <- merge(mm, nn, by = c('Group.1', 'Group.2'))

#names(mm)[1] <- 'Date'

aa <- merge(mm, NR[, c('Date', 'ET')], by = 'Date')
plot(aa$Date, aa$ET.x)
points(aa$Date, aa$ET.y, col = 2, pch = 16)
plot(aa$ET.x, aa$ET.y, xlab = 'Mean ET of East River', ylab = 'ET of US-NR1', pch = 16,
     main = 'Day to day comparison, use with Caustion')

NR_18 <- NR[which(year(NR$Date)==2018),]
NR_18$ET <- SMA(NR_18$ET, 7)
NR_19 <- NR[which(year(NR$Date)==2019),] 
NR_19$ET <- SMA(NR_19$ET, 7)
ggplot() + geom_line(data = NR_18, aes(x = yday(Date), y = ET, color = '2018')) +
  geom_line(data = NR_19, aes(x = yday(Date), y = ET, color = '2019')) +
  ggtitle('Niwot Ridge ET, smoothed') + xlab('Date') + ylab('ET') + ylim(0, 4) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

plot(yday(NR_18$Date), SMA(NR_18$ET, 7), type = 'l')
lines(yday(NR_19$Date), SMA(NR_19$ET, 7), col = 2)
## Time series and radiation
## Let us loop over the albedos
uniq_da <- unique(df_small$Date)

uniq_da18 <- uniq_da[which(year(uniq_da)==2018)][c(3, 5, 7, 9, 11, 15, 17, 20)]
uniq_da19 <- uniq_da[which(year(uniq_da)==2019)][c(3, 5, 7, 8, 10, 12, 16, 19)]
uniq_da18 <- as.Date(c('2018-03-07', '2018-04-15', '2018-05-17', '2018-06-02', '2018-06-27', 
                       '2018-07-29', '2018-08-14', '2018-09-15', '2018-11-09', '2018-12-11'))
uniq_da19 <- as.Date(c('2019-03-17', '2019-04-18', '2019-05-13', '2019-06-05', '2019-06-30', 
                       '2019-07-21', '2019-08-15', '2019-09-14', '2019-11-08', '2019-12-07'))

sub_18 <- df_small[which(df_small$Date %in% uniq_da18),]
sub_19 <- df_small[which(df_small$Date %in% uniq_da19),]
sub_18 <- sub_18[which(sub_18$aspectgroup == 'N' | sub_18$aspectgroup == 'S'),]
sub_19 <- sub_19[which(sub_19$aspectgroup == 'N' | sub_19$aspectgroup == 'S'),]




## Let us select some snowy dates.
snowmeltperiods = 1
if(snowmeltperiods){
  bg_180330 <- df[which(df$Date == as.Date('2018-03-30')),]
  sm_180330 <- df_small[which(df_small$Date == as.Date('2018-03-30')),]
  bg_190326 <- df[which(df$Date == as.Date('2019-03-26')),]
  sm_190326 <- df_small[which(df_small$Date == as.Date('2019-03-26')),]
  
  sm_180330NS <- sm_180330[which(sm_180330$aspectgroup == 'N' | sm_180330$aspectgroup == 'S'),]
  bg_180330NS <- bg_180330[which(bg_180330$aspectgroup == 'N' | bg_180330$aspectgroup == 'S'),]
  sm_190326NS <- sm_190326[which(sm_190326$aspectgroup == 'N' | sm_190326$aspectgroup == 'S'),]
  bg_190326NS <- bg_190326[which(bg_190326$aspectgroup == 'N' | bg_190326$aspectgroup == 'S'),]
  
  ## Energy balances and albedo
  
  ggplot() + geom_raster(data = sm_180330, aes(x = x, y = y, fill = Rn), interpolate = TRUE) + 
    coord_quickmap() + ggtitle('Rn [W/m2]') + scale_fill_gradientn(colours = pal) + 
    geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2018-03-30 Rn') + 
    geom_point(data = ERWR, aes(x = lon, y = lat)) +
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
  
  ggplot() + geom_raster(data = sm_180330, aes(x = x, y = y, fill = SW_out/SW_in), interpolate = TRUE) + 
    coord_quickmap() + ggtitle('Albedo') + scale_fill_gradientn(colours = pal) + 
    geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2018-03-30 albedo') + 
    geom_point(data = ERWR, aes(x = lon, y = lat)) +
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
  
  ggplot() + geom_raster(data = sm_180330, aes(x = x, y = y, fill = SW_in), interpolate = TRUE) + 
    coord_quickmap() + ggtitle('SW_in') + scale_fill_gradientn(colours = pal) + 
    geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2018-03-30 SW_in') + 
    geom_point(data = ERWR, aes(x = lon, y = lat)) +
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
  
  ggplot() + geom_raster(data = sm_180330, aes(x = x, y = y, fill = NDVI), interpolate = TRUE) + 
    coord_quickmap() + ggtitle('NDVI') + scale_fill_gradientn(colours = pal) + 
    geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2018-03-30 NDVI') + 
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
  
  ## Now let us go to 2019 
  ggplot() + geom_raster(data = sm_190326, aes(x = x, y = y, fill = Rn), interpolate = TRUE) + 
    coord_quickmap() + ggtitle('Rn') + scale_fill_gradientn(colours = pal) + 
    geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2019-03-26 Rn') + 
    geom_point(data = ERWR, aes(x = lon, y = lat)) +
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
  
  ggplot() + geom_raster(data = sm_190326, aes(x = x, y = y, fill = SW_in), interpolate = TRUE) + 
    coord_quickmap() + ggtitle('SW_in') + scale_fill_gradientn(colours = pal) + 
    geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2019-03-26 SW_in') + 
    geom_point(data = ERWR, aes(x = lon, y = lat)) +
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
  
  ggplot() + geom_raster(data = sm_190326, aes(x = x, y = y, fill = SW_out/SW_in), interpolate = TRUE) + 
    coord_quickmap() + ggtitle('Albedo') + scale_fill_gradientn(colours = pal) + 
    geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2019-03-26 albedo') + 
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
  
  ggplot() + geom_raster(data = sm_190326, aes(x = x, y = y, fill = NDVI), interpolate = TRUE) + 
    coord_quickmap() + ggtitle('NDVI') + scale_fill_gradientn(colours = pal) + 
    geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2019-03-26 NDVI') + 
    geom_point(data = ERWR, aes(x = lon, y = lat)) +
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
  
}

## Snowmeltperiods
bg_180517 <- df[which(df$Date == as.Date('2018-05-17')),]
sm_180517 <- df_small[which(df_small$Date == as.Date('2018-05-17')),]
bg_190605 <- df[which(df$Date == as.Date('2019-06-05')),]
sm_190605 <- df_small[which(df_small$Date == as.Date('2019-06-05')),]

sm_180517NS <- sm_180517[which(sm_180517$aspectgroup == 'N' | sm_180517$aspectgroup == 'S'),]
bg_180517NS <- bg_180517[which(bg_180517$aspectgroup == 'N' | bg_180517$aspectgroup == 'S'),]
sm_190605NS <- sm_190605[which(sm_190605$aspectgroup == 'N' | sm_190605$aspectgroup == 'S'),]
bg_190605NS <- bg_190605[which(bg_190605$aspectgroup == 'N' | bg_190605$aspectgroup == 'S'),]

## Energy balances and albedo

ggplot() + geom_raster(data = sm_180517, aes(x = x, y = y, fill = Rn), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('Rn [W/m2]') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2018-05-17 Rn') + 
  geom_point(data = ERWR, aes(x = lon, y = lat)) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

ggplot() + geom_raster(data = sm_180517, aes(x = x, y = y, fill = SW_out/SW_in), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('Albedo') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2018-05-17 albedo') + 
  geom_point(data = ERWR, aes(x = lon, y = lat)) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

ggplot() + geom_raster(data = sm_180517, aes(x = x, y = y, fill = SW_in), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('SW_in') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2018-05-17 SW_in') + 
  geom_point(data = ERWR, aes(x = lon, y = lat)) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

ggplot() + geom_raster(data = sm_180517, aes(x = x, y = y, fill = NDVI), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('NDVI') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2018-05-17 NDVI') + 
  geom_point(data = ERWR, aes(x = lon, y = lat)) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

## Now let us go to 2019 
ggplot() + geom_raster(data = sm_190605, aes(x = x, y = y, fill = Rn), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('Rn') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2019-06-05 Rn') + 
  geom_point(data = ERWR, aes(x = lon, y = lat)) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

ggplot() + geom_raster(data = sm_190605, aes(x = x, y = y, fill = SW_in), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('SW_in') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2019-06-05 SW_in') + 
  geom_point(data = ERWR, aes(x = lon, y = lat)) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

ggplot() + geom_raster(data = sm_190605, aes(x = x, y = y, fill = SW_out/SW_in), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('Albedo') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2019-06-05 albedo') + 
  geom_point(data = ERWR, aes(x = lon, y = lat)) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

ggplot() + geom_raster(data = sm_190605, aes(x = x, y = y, fill = NDVI), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('NDVI') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2019-06-05 NDVI') + 
  geom_point(data = ERWR, aes(x = lon, y = lat)) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

#11-02-2020 compare small large domain

ggplot() + geom_raster(data = sm_190605, aes(x = x, y = y, fill = ET), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('NDVI') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2019-06-05 ET') + 
  geom_point(data = ERWR, aes(x = lon, y = lat)) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

ggplot() + geom_raster(data = sm_190605, aes(x = x, y = y, fill = ETs), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('NDVI') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2019-06-05 Evaporation') + 
  geom_point(data = ERWR, aes(x = lon, y = lat)) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 


ggplot() + geom_raster(data = sm_190605, aes(x = x, y = y, fill = ETc), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('NDVI') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2019-06-05 Transpiration') + 
  geom_point(data = ERWR, aes(x = lon, y = lat)) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

ggplot() + geom_raster(data = sm_190605, aes(x = x, y = y, fill = NDVI), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('NDVI') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2019-06-05 NDVI') + 
  geom_point(data = ERWR, aes(x = lon, y = lat)) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

## small domain
ggplot() + geom_raster(data = sm_190605[which(sm_190605$site == 'Pumphouse'),], aes(x = x, y = y, fill = ET), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('NDVI') + scale_fill_gradientn(colours = pal) +   
  geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2019-06-05 ET') + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

ggplot() + geom_raster(data = sm_190605[which(sm_190605$site == 'Pumphouse'),], aes(x = x, y = y, fill = ETs), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('NDVI') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2019-06-05 Evaporation') + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 


ggplot() + geom_raster(data = sm_190605[which(sm_190605$site == 'Pumphouse'),], aes(x = x, y = y, fill = ETc), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('NDVI') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2019-06-05 Transpiration') + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

ggplot() + geom_raster(data = sm_190605[which(sm_190605$site == 'Pumphouse'),], aes(x = x, y = y, fill = NDVI), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('NDVI') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2019-06-05 NDVI') + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 



## Let us loop over the albedos
uniq_da <- unique(df$Date)
uniq_da18 <- uniq_da[which(year(uniq_da)==2018 & month(uniq_da) >=3 & month(uniq_da) <= 10)][c(1, 3, 5, 7, 10, 13, 15)]
uniq_da19 <- uniq_da[which(year(uniq_da)==2019)][c(3, 5, 7, 9, 12, 16, 19)]

load('C:/Users/Nigel/Box/GEE/Ameri/NR.RData')
for(idx in 1:length(uniq_da18)){
  tt <- df_small[which(df_small$Date==uniq_da18[idx]),]
  g <- ggplot() + geom_raster(data = tt, aes(x = x, y = y, fill = SW_out/SW_in), interpolate = TRUE) + 
    coord_quickmap() + ggtitle('Albedo') + scale_fill_gradientn(colours = pal) + 
    geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle(paste('Albedo', uniq_da18[idx])) + 
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
  print(g)
}

for(idx in 1:length(uniq_da19)){
  tt <- df_small[which(df_small$Date==uniq_da19[idx]),]
  g <- ggplot() + geom_raster(data = tt, aes(x = x, y = y, fill = SW_out/SW_in), interpolate = TRUE) + 
    coord_quickmap() + ggtitle('Albedo') + scale_fill_gradientn(colours = pal) + 
    geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle(paste('Albedo', uniq_da19[idx])) + 
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
  print(g)
}


## ET spatial domain
for(idx in 1:length(uniq_da18)){
  tt <- df_small[which(df_small$Date==uniq_da18[idx]),]
  g <- ggplot() + geom_raster(data = tt, aes(x = x, y = y, fill = ET), interpolate = TRUE) + 
    coord_quickmap() + ggtitle('ET') + scale_fill_gradientn(colours = pal) + 
    geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle(paste('ET', uniq_da18[idx])) + 
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
  print(g)
}

for(idx in 1:length(uniq_da19)){
  tt <- df_small[which(df_small$Date==uniq_da19[idx]),]
  g <- ggplot() + geom_raster(data = tt, aes(x = x, y = y, fill = ET), interpolate = TRUE) + 
    coord_quickmap() + ggtitle('ET') + scale_fill_gradientn(colours = pal) + 
    geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle(paste('ET', uniq_da19[idx])) + 
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
  print(g)
}


## Transpiration
for(idx in 1:length(uniq_da18)){
  tt <- df_small[which(df_small$Date==uniq_da18[idx]),]
  g <- ggplot() + geom_raster(data = tt, aes(x = x, y = y, fill = ETc), interpolate = TRUE) + 
    coord_quickmap() + ggtitle('Transpiration') + scale_fill_gradientn(colours = pal) + 
    geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle(paste('Transpiration', uniq_da18[idx])) + 
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
  print(g)
}

for(idx in 1:length(uniq_da19)){
  tt <- df_small[which(df_small$Date==uniq_da19[idx]),]
  g <- ggplot() + geom_raster(data = tt, aes(x = x, y = y, fill = ETc), interpolate = TRUE) + 
    coord_quickmap() + ggtitle('Transpiration') + scale_fill_gradientn(colours = pal) + 
    geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle(paste('Transpiration', uniq_da19[idx])) + 
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
  print(g)
}

## soil evaporation
for(idx in 1:length(uniq_da18)){
  tt <- df_small[which(df_small$Date==uniq_da18[idx]),]
  g <- ggplot() + geom_raster(data = tt, aes(x = x, y = y, fill = ETs), interpolate = TRUE) + 
    coord_quickmap() + ggtitle('Soil evaporation') + scale_fill_gradientn(colours = pal) + 
    geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle(paste('Soil evaporation', uniq_da18[idx])) + 
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
  print(g)
}

for(idx in 1:length(uniq_da19)){
  tt <- df_small[which(df_small$Date==uniq_da19[idx]),]
  g <- ggplot() + geom_raster(data = tt, aes(x = x, y = y, fill = ETs), interpolate = TRUE) + 
    coord_quickmap() + ggtitle('Soil evaporation') + scale_fill_gradientn(colours = pal) + 
    geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle(paste('Soil evaporation', uniq_da19[idx])) + 
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
  print(g)
}

## ET spatial domain
for(idx in 1:length(uniq_da18)){
  tt <- df_small[which(df_small$Date==uniq_da18[idx]),]
  g <- ggplot() + geom_raster(data = tt, aes(x = x, y = y, fill = NDVI), interpolate = TRUE) + 
    coord_quickmap() + ggtitle('NDVI') + scale_fill_gradientn(colours = pal) + 
    geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle(paste('NDVI', uniq_da18[idx])) + 
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
  print(g)
}

for(idx in 1:length(uniq_da19)){
  tt <- df_small[which(df_small$Date==uniq_da19[idx]),]
  g <- ggplot() + geom_raster(data = tt, aes(x = x, y = y, fill = NDVI), interpolate = TRUE) + 
    coord_quickmap() + ggtitle('NDVI') + scale_fill_gradientn(colours = pal) + 
    geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle(paste('NDVI', uniq_da19[idx])) + 
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
  print(g)
}

## Peak growing season periods
## Snowmeltperiods
bg_180729 <- df[which(df$Date == as.Date('2018-07-29')),]
sm_180729 <- df_small[which(df_small$Date == as.Date('2018-07-29')),]
bg_190723 <- df[which(df$Date == as.Date('2019-07-23')),]
sm_190723 <- df_small[which(df_small$Date == as.Date('2019-07-23')),]

sm_180729NS <- sm_180729[which(sm_180729$aspectgroup == 'N' | sm_180729$aspectgroup == 'S'),]
bg_180729NS <- bg_180729[which(bg_180729$aspectgroup == 'N' | bg_180729$aspectgroup == 'S'),]
sm_190723NS <- sm_190723[which(sm_190723$aspectgroup == 'N' | sm_190723$aspectgroup == 'S'),]
bg_190723NS <- bg_190723[which(bg_190723$aspectgroup == 'N' | bg_190723$aspectgroup == 'S'),]

## energy
ggplot() + geom_raster(data = sm_180729, aes(x = x, y = y, fill = Rn), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('Rn') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2018-07-29 Rn') + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

ggplot() + geom_raster(data = sm_190723, aes(x = x, y = y, fill = Rn), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('Rn') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + ggtitle('2019-07-23 Rn') + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 




# ## Time series and radiation
# uniq_da18 <- uniq_da[which(year(uniq_da)==2018)][c(3, 5, 7, 9, 11, 15, 17, 20)]
# uniq_da19 <- uniq_da[which(year(uniq_da)==2019)][c(3, 5, 7, 8, 10, 12, 16, 19)]
uniq_da18 <- as.Date(c('2018-01-02', '2018-02-26', '2018-03-07', '2018-04-15', '2018-05-17', 
                       '2018-06-02', '2018-06-27', '2018-07-21', '2018-08-14', '2018-09-22', 
                       '2018-10-19', '2018-11-09', '2018-12-11'))
uniq_da19 <- as.Date(c('2019-01-05', '2019-02-01', '2019-03-17', '2019-04-18', '2019-05-13', 
                       '2019-06-05', '2019-06-30', '2019-07-21', '2019-08-15', '2019-09-14',
                       '2019-10-14', '2019-11-08', '2019-12-07'))

sub_18 <- df_small[which(df_small$Date %in% uniq_da18),]
sub_19 <- df_small[which(df_small$Date %in% uniq_da19),]
sub_18 <- sub_18[which(sub_18$aspectgroup == 'N' | sub_18$aspectgroup == 'S'),]
sub_19 <- sub_19[which(sub_19$aspectgroup == 'N' | sub_19$aspectgroup == 'S'),]


small_18 <- sm[which(sm$Date %in% uniq_da18),]
small_18 <- small_18[-which(small_18$Date == as.Date('2018-05-17') & small_18$sate == 'S2'),]
small_19 <- sm[which(sm$Date %in% uniq_da19),]

### focus on the pumphouse domain
small_18 <- small_18[which(small_18$aspectgroup %in% c('N', 'S', 'NE', 'SW')),]
small_19 <- small_19[which(small_19$aspectgroup %in% c('N', 'S', 'NE', 'SW')),]

## 01-04-2021 Raster plots to show slope and aspect together

## 01-04-2021 Raster plots to show slope and aspect together

## 01-04-2021 Raster plots to show slope and aspect together

## 01-04-2021 Raster plots to show slope and aspect together

## 01-04-2021 Raster plots to show slope and aspect together

## 01-04-2021 Raster plots to show slope and aspect together

## Dominance analysis
library(dominanceanalysis)
listdates <- c(uniq_da18[c(-1, -13)], uniq_da19[c(-1, -13)])
ETc_contr_S <- data.frame(Date = listdates, Rn = 0, SAVI = 0, aspectgroup = 'S')
ETc_contr_N <- data.frame(Date = listdates, Rn = 0, SAVI = 0, aspectgroup = 'N')
#ET_contr_S <- data.frame(Date = listdates, Rn = 0, SAVI = 0, f_SM = 0, VP = 0, aspectgroup = 'S')
#ET_contr_N <- data.frame(Date = listdates, Rn = 0, SAVI = 0, f_SM = 0, VP = 0, aspectgroup = 'N')

for(udx in 1:length(listdates)){
  temp_S <- df_small[which((df_small$Date %in% listdates[udx]) & df_small$aspectgroup %in% c('SW','S')),]
  lm.temp_S <- lm(ETc~Rn + SAVI, temp_S)
  da.temp_S <- dominanceAnalysis(lm.temp_S)
  ETc_contr_S[udx, 2:3] <- da.temp_S$contribution.average$r2

  temp_N <- df_small[which((df_small$Date %in% listdates[udx]) & df_small$aspectgroup %in% c('N','NE')),]
  lm.temp_N <- lm(ETc~Rn + SAVI, temp_N)
  da.temp_N <- dominanceAnalysis(lm.temp_N)
  ETc_contr_N[udx, 2:3] <- da.temp_N$contribution.average$r2
  print(listdates[udx])
}
ETc_contr <- rbind(ETc_contr_S, ETc_contr_N)
ETc_contr$Rn <- ETc_contr$Rn/(ETc_contr$Rn + ETc_contr$SAVI)
ETc_contr$SAVI <- 1 - ETc_contr$Rn
ETc_contr$yr <- year(ETc_contr$Date)
library(reshape2)
ETc_contr_m <- melt(ETc_contr, id = c('Date', 'aspectgroup', 'yr'))
ETc_contr_18 <- ETc_contr_m[which(ETc_contr_m$yr==2018),]
ETc_contr_19 <- ETc_contr_m[which(ETc_contr_m$yr==2019),]
ggplot(ETc_contr_18, aes(x = as.factor(Date), y = value)) + 
  geom_bar(aes(fill = variable), stat = 'identity', width = 0.8) +
  facet_grid(aspectgroup~.) + ggtitle('Transpiration dominance in year 2018') + 
  xlab('Date') + ylab('Contribution') + labs(fill = 'Attributes') + 
  theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold"),
        axis.text.x = element_text(angle = 90, size = 10),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20))

ggplot(ETc_contr_19, aes(x = as.factor(Date), y = value)) + 
  geom_bar(aes(fill = variable), stat = 'identity', width = 0.8) +
  facet_grid(aspectgroup~.) + ggtitle('Transpiration dominance in year 2019') + 
  xlab('Date') + ylab('Contribution') + labs(fill = 'Attributes') + 
  theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold"),
        axis.text.x = element_text(angle = 90, size = 10),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20))

### ETs and 
ETs_contr_S <- data.frame(Date = listdates, Rn = 0, SAVI = 0, aspectgroup = 'S')
ETs_contr_N <- data.frame(Date = listdates, Rn = 0, SAVI = 0, aspectgroup = 'N')
#ET_contr_S <- data.frame(Date = listdates, Rn = 0, SAVI = 0, f_SM = 0, VP = 0, aspectgroup = 'S')
#ET_contr_N <- data.frame(Date = listdates, Rn = 0, SAVI = 0, f_SM = 0, VP = 0, aspectgroup = 'N')

for(udx in 1:length(listdates)){
  temp_S <- df_small[which((df_small$Date %in% listdates[udx]) & df_small$aspectgroup %in% c('SW','S')),]
  lm.temp_S <- lm(ETs~Rn + SAVI, temp_S)
  da.temp_S <- dominanceAnalysis(lm.temp_S)
  ETs_contr_S[udx, 2:3] <- da.temp_S$contribution.average$r2
  
  temp_N <- df_small[which((df_small$Date %in% listdates[udx]) & df_small$aspectgroup %in% c('N','NE')),]
  lm.temp_N <- lm(ETs~Rn + SAVI, temp_N)
  da.temp_N <- dominanceAnalysis(lm.temp_N)
  ETs_contr_N[udx, 2:3] <- da.temp_N$contribution.average$r2
  print(listdates[udx])
}
ETs_contr <- rbind(ETs_contr_S, ETs_contr_N)
ETs_contr$Rn <- ETs_contr$Rn/(ETs_contr$Rn + ETs_contr$SAVI)
ETs_contr$SAVI <- 1 - ETs_contr$Rn
ETs_contr$yr <- year(ETs_contr$Date)
library(reshape2)
ETs_contr_m <- melt(ETs_contr, id = c('Date', 'aspectgroup', 'yr'))
ETs_contr_18 <- ETs_contr_m[which(ETs_contr_m$yr==2018),]
ETs_contr_19 <- ETs_contr_m[which(ETs_contr_m$yr==2019),]
ggplot(ETs_contr_18, aes(x = as.factor(Date), y = value)) + 
  geom_bar(aes(fill = variable), stat = 'identity', width = 0.8) +
  facet_grid(aspectgroup~.) + ggtitle('Evaporation dominance in year 2018') + 
  xlab('Date') + ylab('Contribution') + labs(fill = 'Attributes') + 
  theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold"),
        axis.text.x = element_text(angle = 90, size = 10),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20))

ggplot(ETs_contr_19, aes(x = as.factor(Date), y = value)) + 
  geom_bar(aes(fill = variable), stat = 'identity', width = 0.8) +
  facet_grid(aspectgroup~.) + ggtitle('Evaporation dominance in year 2019') + 
  xlab('Date') + ylab('Contribution') + labs(fill = 'Attributes') + 
  theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold"),
        axis.text.x = element_text(angle = 90, size = 10),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20))

### ET components
ET_contr_S <- data.frame(Date = listdates, Rn = 0, SAVI = 0, aspectgroup = 'S')
ET_contr_N <- data.frame(Date = listdates, Rn = 0, SAVI = 0, aspectgroup = 'N')
#ET_contr_S <- data.frame(Date = listdates, Rn = 0, SAVI = 0, f_SM = 0, VP = 0, aspectgroup = 'S')
#ET_contr_N <- data.frame(Date = listdates, Rn = 0, SAVI = 0, f_SM = 0, VP = 0, aspectgroup = 'N')

for(udx in 1:length(listdates)){
  temp_S <- df_small[which((df_small$Date %in% listdates[udx]) & df_small$aspectgroup %in% c('SW','S')),]
  lm.temp_S <- lm(ET~Rn + SAVI, temp_S)
  da.temp_S <- dominanceAnalysis(lm.temp_S)
  ET_contr_S[udx, 2:3] <- da.temp_S$contribution.average$r2
  
  temp_N <- df_small[which((df_small$Date %in% listdates[udx]) & df_small$aspectgroup %in% c('N','NE')),]
  lm.temp_N <- lm(ET~Rn + SAVI, temp_N)
  da.temp_N <- dominanceAnalysis(lm.temp_N)
  ET_contr_N[udx, 2:3] <- da.temp_N$contribution.average$r2
  print(listdates[udx])
}
ET_contr <- rbind(ET_contr_S, ET_contr_N)
ET_contr$Rn <- ET_contr$Rn/(ET_contr$Rn + ET_contr$SAVI)
ET_contr$SAVI <- 1 - ET_contr$Rn
ET_contr$yr <- year(ET_contr$Date)
library(reshape2)
ET_contr_m <- melt(ET_contr, id = c('Date', 'aspectgroup', 'yr'))
ET_contr_18 <- ET_contr_m[which(ET_contr_m$yr==2018),]
ET_contr_19 <- ET_contr_m[which(ET_contr_m$yr==2019),]
ggplot(ET_contr_18, aes(x = as.factor(Date), y = value)) + 
  geom_bar(aes(fill = variable), stat = 'identity', width = 0.8) +
  facet_grid(aspectgroup~.) + ggtitle('ET dominance in year 2018') + 
  xlab('Date') + ylab('Contribution') + labs(fill = 'Attributes') + 
  theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold"),
        axis.text.x = element_text(angle = 90, size = 10),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20))

ggplot(ET_contr_19, aes(x = as.factor(Date), y = value)) + 
  geom_bar(aes(fill = variable), stat = 'identity', width = 0.8) +
  facet_grid(aspectgroup~.) + ggtitle('ET dominance in year 2019') + 
  xlab('Date') + ylab('Contribution') + labs(fill = 'Attributes') + 
  theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold"),
        axis.text.x = element_text(angle = 90, size = 10),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20))

## Correlation coefficient between ET and other variables

library(plyr)
#not including Jan and Dec data
Sou_coeff <- df_small[which((df_small$Date %in% c(uniq_da18[2:12], uniq_da19[2:12])) & (df_small$aspectgroup %in% c('S', 'SW'))),]
Nor_coeff <- df_small[which((df_small$Date %in% c(uniq_da18[2:12], uniq_da19[2:12])) & (df_small$aspectgroup %in% c('N', 'NE'))),]

S.ET_SAVI <- ddply(Sou_coeff, 'Date', function(x)cor(x$ET, x$SAVI))
S.ET_Rn <- ddply(Sou_coeff, 'Date', function(x)cor(x$ET, x$Rn))
N.ET_SAVI <- ddply(Nor_coeff, 'Date', function(x)cor(x$ET, x$SAVI))
N.ET_Rn <- ddply(Nor_coeff, 'Date', function(x)cor(x$ET, x$Rn))

Coeff <- data.frame(Date = c(S.ET_SAVI$Date, N.ET_SAVI$Date),
                    SAVI = c(S.ET_SAVI$V1, N.ET_SAVI$V1),
                    Rn = c(S.ET_Rn$V1, N.ET_Rn$V1),
                    aspect = c(rep('S', 11), rep('N', 11)))
Coeff$doy <- yday(Coeff$Date)
Coeff$yr <- year(Coeff$Date)
ggplot(Coeff, aes(x = doy, y = Rn, group = aspect, shape = aspect, colour = factor(yr))) + 
  geom_line() + geom_point() + 
  facet_grid(facets = factor(yr)~.) + theme_bw()

## ET rasters. also look at slopes
ggplot() + geom_boxplot(data = small_18, aes(x = as.factor(Date), y = NDVI, fill = aspectgroup)) +
  ggtitle('Temporal distribution of NDVI in 2018') + xlab('Date') + ylab('NDVI') + ylim(0, 0.8) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90))
ggplot() + geom_boxplot(data = small_19, aes(x = as.factor(Date), y = NDVI, fill = aspectgroup)) +
  ggtitle('Temporal distribution of NDVI in 2019') + xlab('Date') + ylab('NDVI') + ylim(0, 0.8) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90))

ggplot() + geom_boxplot(data = small_18, aes(x = as.factor(Date), y = Rn, fill = aspectgroup)) +
  ggtitle('Temporal distribution of Rn in 2018') + xlab('Date') + ylab('Rn')  + ylim(0, 250) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90))
ggplot() + geom_boxplot(data = small_19, aes(x = as.factor(Date), y = Rn, fill = aspectgroup)) +
  ggtitle('Temporal distribution of Rn in 2019') + xlab('Date') + ylab('Rn') + ylim(0, 250) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90))

ggplot() + geom_boxplot(data = small_18, aes(x = as.factor(Date), y = ET, fill = aspectgroup)) +
  ggtitle('Temporal distribution of ET in 2018') + xlab('Date') + ylab('ET') + ylim(0, 3.5) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90))
ggplot() + geom_boxplot(data = small_19, aes(x = as.factor(Date), y = ET, fill = aspectgroup)) +
  ggtitle('Temporal distribution of ET in 2019') + xlab('Date') + ylab('ET') + ylim(0, 3.5) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90))

## ETc

ggplot() + geom_boxplot(data = small_18, aes(x = as.factor(Date), y = ETc, fill = compgr)) +
  ggtitle('Temporal distribution of transpiration in 2018') + xlab('Date') + ylab('ETc') + ylim(0, 3.5) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90))
ggplot() + geom_boxplot(data = small_19, aes(x = as.factor(Date), y = ETc, fill = compgr)) +
  ggtitle('Temporal distribution of transpiration in 2019') + xlab('Date') + ylab('ETc') + ylim(0, 3.5) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90))

ggplot() + geom_boxplot(data = small_18, aes(x = as.factor(Date), y = ETs, fill = compgr)) +
  ggtitle('Temporal distribution of evaporation in 2018') + xlab('Date') + ylab('ETc') + ylim(0, 3) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90))
ggplot() + geom_boxplot(data = small_19, aes(x = as.factor(Date), y = ETs, fill = compgr)) +
  ggtitle('Temporal distribution of evaporation in 2019') + xlab('Date') + ylab('ETc') + ylim(0, 3) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90))

## North and south only
ggplot() + geom_boxplot(data = small_18[which(small_18$aspectgroup=='N'),], aes(x = as.factor(Date), y = ET, fill = slopegroup)) +
  ggtitle('2018 ET on North slopes') + xlab('Date') + ylab('ET') + ylim(0, 3.5) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90))
ggplot() + geom_boxplot(data = small_18[which(small_18$aspectgroup=='S'),], aes(x = as.factor(Date), y = ET, fill = slopegroup)) +
  ggtitle('2018 ET on South slopes') + xlab('Date') + ylab('ET') + ylim(0, 3.5) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90))

ggplot() + geom_boxplot(data = small_19[which(small_19$aspectgroup=='N'),], aes(x = as.factor(Date), y = ET, fill = slopegroup)) +
  ggtitle('2019 ET on North slopes') + xlab('Date') + ylab('ET') + ylim(0, 3.5) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90))
ggplot() + geom_boxplot(data = small_19[which(small_19$aspectgroup=='S'),], aes(x = as.factor(Date), y = ET, fill = slopegroup)) +
  ggtitle('2019 ET on South slopes') + xlab('Date') + ylab('ET') + ylim(0, 3.5) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90))



## Radiation versus ET
slo_18 <- aggregate(small_18[, c('Rn', 'SAVI', 'ETi', 'ETs', 'ETc', 'ET')], by = list(small_18$Date, small_18$slopegroup, small_18$aspectgroup), mean)
names(slo_18)[1:3] <- c('Date', 'slopegroup', 'aspectgroup')
slo_19 <- aggregate(small_19[, c('Rn', 'SAVI', 'ETi', 'ETs', 'ETc', 'ET')], by = list(small_19$Date, small_19$slopegroup, small_19$aspectgroup), mean)
names(slo_19)[1:3] <- c('Date', 'slopegroup', 'aspectgroup')

ggplot() + geom_point(data = small_18[which(small_18$aspectgroup=='N'),],
                      aes(x = Rn, y = ET, color = as.factor(Date), shape = slopegroup, size = 0.5)) +
  geom_point(data = slo_18[which(slo_18$aspectgroup == 'N'),],
             aes(x = Rn, y = ET, color = as.factor(Date), shape = slopegroup, size = 1.3)) + 
  ggtitle('2018 ET versus Radiation, North-slopes') + xlab('Rn') + ylab('ET') + ylim(0, 3.5) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90)) +
  guides(alpha = FALSE, size = FALSE) + 
  labs(color = 'Date')
ggplot() + geom_point(data = na.omit(small_19[which(small_19$aspectgroup=='N'),]),
                      aes(x = Rn, y = ET, color = as.factor(Date), shape = slopegroup, size = 0.5)) +
  geom_point(data = slo_19[which(slo_19$aspectgroup == 'N'),],
             aes(x = Rn, y = ET, color = as.factor(Date), shape = slopegroup, size = 1.3)) + 
  ggtitle('2019 ET versus Radiation, North-slopes') + xlab('Rn') + ylab('ET') + ylim(0, 3.5) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90)) +
  guides(alpha = FALSE, size = FALSE) + 
  labs(color = 'Date')

ggplot() + geom_point(data = na.omit(small_18[which(small_18$aspectgroup=='S'),]),
                      aes(x = Rn, y = ET, color = as.factor(Date), shape = slopegroup, size = 0.5)) +
  geom_point(data = slo_18[which(slo_18$aspectgroup == 'S'),],
             aes(x = Rn, y = ET, color = as.factor(Date), shape = slopegroup, size = 1.3)) + 
  ggtitle('2018 ET versus Radiation, South-slopes') + xlab('Rn') + ylab('ET') + ylim(0, 3.5) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90)) +
  guides(alpha = FALSE, size = FALSE) + 
  labs(color = 'Date')
ggplot() + geom_point(data = na.omit(small_19[which(small_18$aspectgroup=='S'),]),
                      aes(x = Rn, y = ET, color = as.factor(Date), shape = slopegroup, size = 0.5)) +
  geom_point(data = slo_19[which(slo_19$aspectgroup == 'S'),],
             aes(x = Rn, y = ET, color = as.factor(Date), shape = slopegroup, size = 1.3)) + 
  ggtitle('2019 ET versus Radiation, South-slopes') + xlab('Rn') + ylab('ET') + ylim(0, 3.5) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90)) +
  guides(alpha = FALSE, size = FALSE) + 
  labs(color = 'Date')


## ETs
ggplot() + geom_point(data = small_18[which(small_18$aspectgroup=='N'),],
                      aes(x = Rn, y = ETs, color = as.factor(Date), shape = slopegroup, size = 0.5)) +
  geom_point(data = slo_18[which(slo_18$aspectgroup == 'N'),],
             aes(x = Rn, y = ETs, color = as.factor(Date), shape = slopegroup, size = 1.3)) + 
  ggtitle('2018 evaporation versus Radiation, North-slopes') + xlab('Rn') + ylab('ET') + ylim(0, 3) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90)) +
  guides(alpha = FALSE, size = FALSE) + 
  labs(color = 'Date')
ggplot() + geom_point(data = na.omit(small_19[which(small_19$aspectgroup=='N'),]),
                      aes(x = Rn, y = ETs, color = as.factor(Date), shape = slopegroup, size = 0.5)) +
  geom_point(data = slo_19[which(slo_19$aspectgroup == 'N'),],
             aes(x = Rn, y = ETs, color = as.factor(Date), shape = slopegroup, size = 1.3)) + 
  ggtitle('2019 evaporation versus Radiation, North-slopes') + xlab('Rn') + ylab('ET') + ylim(0, 3) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90)) +
  guides(alpha = FALSE, size = FALSE) + 
  labs(color = 'Date')

ggplot() + geom_point(data = na.omit(small_18[which(small_18$aspectgroup=='S'),]),
                      aes(x = Rn, y = ETs, color = as.factor(Date), shape = slopegroup, size = 0.5)) +
  geom_point(data = slo_18[which(slo_18$aspectgroup == 'S'),],
             aes(x = Rn, y = ETs, color = as.factor(Date), shape = slopegroup, size = 1.3)) + 
  ggtitle('2018 evaporation versus Radiation, South-slopes') + xlab('Rn') + ylab('ET') + ylim(0, 3) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90)) +
  guides(alpha = FALSE, size = FALSE) + 
  labs(color = 'Date')
ggplot() + geom_point(data = na.omit(small_19[which(small_18$aspectgroup=='S'),]),
                      aes(x = Rn, y = ETs, color = as.factor(Date), shape = slopegroup, size = 0.5)) +
  geom_point(data = slo_19[which(slo_19$aspectgroup == 'S'),],
             aes(x = Rn, y = ETs, color = as.factor(Date), shape = slopegroup, size = 1.3)) + 
  ggtitle('2019 evaporation versus Radiation, South-slopes') + xlab('Rn') + ylab('ET') + ylim(0, 3) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90)) +
  guides(alpha = FALSE, size = FALSE) + 
  labs(color = 'Date')


## ETc

ggplot() + geom_point(data = small_18[which(small_18$aspectgroup=='N'),],
                      aes(x = Rn, y = ETc, color = as.factor(Date), shape = slopegroup, size = 0.5)) +
  geom_point(data = slo_18[which(slo_18$aspectgroup == 'N'),],
             aes(x = Rn, y = ETc, color = as.factor(Date), shape = slopegroup, size = 1.3)) + 
  ggtitle('2018 transpiration versus Radiation, North-slopes') + xlab('Rn') + ylab('ET') + ylim(0, 3) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90)) +
  guides(alpha = FALSE, size = FALSE) + 
  labs(color = 'Date')
ggplot() + geom_point(data = na.omit(small_19[which(small_19$aspectgroup=='N'),]),
                      aes(x = Rn, y = ETc, color = as.factor(Date), shape = slopegroup, size = 0.5)) +
  geom_point(data = slo_19[which(slo_19$aspectgroup == 'N'),],
             aes(x = Rn, y = ETc, color = as.factor(Date), shape = slopegroup, size = 1.3)) + 
  ggtitle('2019 transpiration versus Radiation, North-slopes') + xlab('Rn') + ylab('ET') + ylim(0, 3) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90)) +
  guides(alpha = FALSE, size = FALSE) + 
  labs(color = 'Date')

ggplot() + geom_point(data = na.omit(small_18[which(small_18$aspectgroup=='S'),]),
                      aes(x = Rn, y = ETc, color = as.factor(Date), shape = slopegroup, size = 0.5)) +
  geom_point(data = slo_18[which(slo_18$aspectgroup == 'S'),],
             aes(x = Rn, y = ETc, color = as.factor(Date), shape = slopegroup, size = 1.3)) + 
  ggtitle('2018 transpiration versus Radiation, South-slopes') + xlab('Rn') + ylab('ET') + ylim(0, 3) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90)) +
  guides(alpha = FALSE, size = FALSE) + 
  labs(color = 'Date')
ggplot() + geom_point(data = na.omit(small_19[which(small_18$aspectgroup=='S'),]),
                      aes(x = Rn, y = ETc, color = as.factor(Date), shape = slopegroup, size = 0.5)) +
  geom_point(data = slo_19[which(slo_19$aspectgroup == 'S'),],
             aes(x = Rn, y = ETc, color = as.factor(Date), shape = slopegroup, size = 1.3)) + 
  ggtitle('2019 transpiration versus Radiation, South-slopes') + xlab('Rn') + ylab('ET') + ylim(0, 3) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90)) +
  guides(alpha = FALSE, size = FALSE) + 
  labs(color = 'Date')

## TRANSPIRATION MOSTLY DEPENDS ON VEGETATION STATUS

ggplot() + geom_point(data = na.omit(small_18[which(small_18$aspectgroup=='S'),]),
                      aes(x = SAVI, y = ETc, color = as.factor(Date), shape = slopegroup, size = 0.5)) +
  geom_point(data = slo_18[which(slo_18$aspectgroup == 'S'),],
             aes(x = SAVI, y = ETc, color = as.factor(Date), shape = slopegroup, size = 1.3)) + 
  ggtitle('2018 transpiration versus SAVI, South-slopes') + xlab('SAVI') + ylab('ET') + ylim(0, 3) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90)) +
  guides(alpha = FALSE, size = FALSE) + 
  labs(color = 'Date')

ggplot() + geom_point(data = na.omit(small_19[which(small_19$aspectgroup=='S'),]),
                      aes(x = SAVI, y = ETc, color = as.factor(Date), shape = slopegroup, size = 0.5)) +
  geom_point(data = slo_19[which(slo_19$aspectgroup == 'S'),],
             aes(x = SAVI, y = ETc, color = as.factor(Date), shape = slopegroup, size = 1.3)) + 
  ggtitle('2019 transpiration versus SAVI, South-slopes') + xlab('SAVI') + ylab('ET') + ylim(0, 3) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90)) +
  guides(alpha = FALSE, size = FALSE) + 
  labs(color = 'Date')

ggplot() + geom_point(data = na.omit(small_18[which(small_18$aspectgroup=='N'),]),
                      aes(x = SAVI, y = ETc, color = as.factor(Date), shape = slopegroup, size = 0.5)) +
  geom_point(data = slo_18[which(slo_18$aspectgroup == 'N'),],
             aes(x = SAVI, y = ETc, color = as.factor(Date), shape = slopegroup, size = 1.3)) + 
  ggtitle('2018 transpiration versus SAVI, North-slopes') + xlab('SAVI') + ylab('ET') + ylim(0, 3) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90)) +
  guides(alpha = FALSE, size = FALSE) + 
  labs(color = 'Date')

ggplot() + geom_point(data = na.omit(small_19[which(small_19$aspectgroup=='N'),]),
                      aes(x = SAVI, y = ETc, color = as.factor(Date), shape = slopegroup, size = 0.5)) +
  geom_point(data = slo_19[which(slo_19$aspectgroup == 'N'),],
             aes(x = SAVI, y = ETc, color = as.factor(Date), shape = slopegroup, size = 1.3)) + 
  ggtitle('2019 transpiration versus SAVI, North-slopes') + xlab('SAVI') + ylab('ET') + ylim(0, 3) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90)) +
  guides(alpha = FALSE, size = FALSE) + 
  labs(color = 'Date')
#########################################
#########################################
#########################################
#########################################

## Plot the ET rasters, larger domain

ggplot() + geom_boxplot(data = sub_18, aes(x = as.factor(Date), y = NDVI, fill = aspectgroup)) +
  ggtitle('Temporal distribution of NDVI') + xlab('Date') + ylab('NDVI') + ylim(0, 0.8) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggplot() + geom_boxplot(data = sub_19, aes(x = as.factor(Date), y = NDVI, fill = aspectgroup))  +
  ggtitle('Temporal distribution of NDVI') + xlab('Date') + ylab('NDVI') + ylim(0, 0.8) +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

## ET components

ggplot() + geom_boxplot(data = sub_18, aes(x = as.factor(Date), y = ETc, fill = aspectgroup)) +
  ggtitle('Temporal distribution of transpiration') + xlab('Date') + ylab('ETc') + ylim(0, 4) +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold"),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90))


ggplot() + geom_boxplot(data = sub_19, aes(x = as.factor(Date), y = ETc, fill = aspectgroup))  +
  ggtitle('Temporal distribution of transpiration') + xlab('Date') + ylab('ETc') + ylim(0, 4) +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold"),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20))

ggplot() + geom_boxplot(data = sub_18, aes(x = as.factor(Date), y = ETs, fill = aspectgroup)) +
  ggtitle('Temporal distribution of evaporation') + xlab('Date') + ylab('ETs') + ylim(0, 4) +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold"),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20))

ggplot() + geom_boxplot(data = sub_19, aes(x = as.factor(Date), y = ETs, fill = aspectgroup))  +
  ggtitle('Temporal distribution of evaporation') + xlab('Date') + ylab('ETs') + ylim(0, 4) +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold"),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20))

ggplot() + geom_boxplot(data = sub_18, aes(x = as.factor(Date), y = 100 * ETc/ET, fill = aspectgroup))  +
  ggtitle('Temporal distribution of Tran%') + xlab('Date') + ylab('Tran%') +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold"),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20))

ggplot() + geom_boxplot(data = sub_19, aes(x = as.factor(Date), y = 100 * ETc/ET, fill = aspectgroup))  +
  ggtitle('Temporal distribution of Tran%') + xlab('Date') + ylab('Tran%') +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold"),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20))

## Moisture and temperature constraints
ggplot() + geom_boxplot(data = sub_18, aes(x = as.factor(Date), y = f_T, fill = aspectgroup))  +
  ggtitle('plant temperature constraint') + xlab('Date') + ylab('f_T') + ylim(0, 1) + 
  theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold"),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20))

ggplot() + geom_boxplot(data = sub_19, aes(x = as.factor(Date), y = f_T, fill = aspectgroup))  +
  ggtitle('plant temperature constraint') + xlab('Date') + ylab('f_T') + ylim(0, 1) + 
  theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold"),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20))


ggplot() + geom_boxplot(data = sub_18, aes(x = as.factor(Date), y = f_SM, fill = aspectgroup))  +
  ggtitle('soil moisture constraint') + xlab('Date') + ylab('f_sm') + ylim(0, 1) + 
  theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold"),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20))

ggplot() + geom_boxplot(data = sub_19, aes(x = as.factor(Date), y = f_SM, fill = aspectgroup))  +
  ggtitle('soil moisture constraint') + xlab('Date') + ylab('f_sm') + ylim(0, 1) + 
  theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold"),
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20))


## savi/ndvi changes over time (rate of change)
## spatial dependent temporal increase rate of NDVI
df_small <- df[which(df$x >= -106.9555 & df$x <= -106.9434 & df$y >= 38.91701 & df$y <= 38.9274),]

Date_cd_18 <- c(as.Date('2018-04-15'), as.Date('2018-05-17'), as.Date('2018-06-02'), as.Date('2018-06-11'))
Date_cd_19 <- c(as.Date('2019-04-18'), as.Date('2019-05-13'), as.Date('2019-06-06'), as.Date('2019-06-30'))

a307_18 <- df_small[which(df_small$Date == as.Date('2018-03-07')), ]
a415_18 <- df_small[which(df_small$Date == as.Date('2018-04-15')), ]
# There's 05-17-2018 from both Landsat and Sentinel2
a517_18 <- df_small[which(df_small$Date == as.Date('2018-05-17') & df_small$sate == 'L8'), ]
 
a611_18 <- df_small[which(df_small$Date == as.Date('2018-06-11')), ]
a704_18 <- df_small[which(df_small$Date == as.Date('2018-07-04')), ]


a310_19 <- df_small[which(df_small$Date == as.Date('2019-03-10')), ]
a418_18 <- df_small[which(df_small$Date == as.Date('2019-04-18')), ]
a513_19 <- df_small[which(df_small$Date == as.Date('2019-05-13')), ]
a605_19 <- df_small[which(df_small$Date == as.Date('2019-06-05')), ]
a716_19 <- df_small[which(df_small$Date == as.Date('2019-07-16') & df_small$sate == 'L8'), ]


cd_18 <- cbind(a307_18[, c('x', 'y', 'slopegroup', 'aspectgroup', 'SAVI')], 
               a415_18[, 'SAVI'], a517_18[, 'SAVI'], a611_18[, 'SAVI'],  a704_18[, 'SAVI'])
names(cd_18) <- c('x', 'y', 'slopegroup', 'aspectgroup', 'VI180307', 'VI180415',
                  'VI180517', 'VI180611', 'VI180704')
cd_19 <- cbind(a310_19[, c('x', 'y', 'slopegroup', 'aspectgroup', 'SAVI')], 
               a418_18[, 'SAVI'], a513_19[, 'SAVI'], a605_19[, 'SAVI'], a716_19[, 'SAVI'])
names(cd_19) <- c('x', 'y', 'slopegroup', 'aspectgroup', 'VI190310', 'VI190418', 'VI190513',
                  'VI190605', 'VI190716')

cd <- cbind(cd_18, cd_19[, c(-1, -2, -3, -4)])
ggplot() + geom_point(data = cd, aes(x = VI180307, y = VI190310, color = '2018-03-07 vs. 2019-03-10', shape = aspectgroup)) + 
  geom_point(data = cd, aes(x = VI180415, y = VI190418, color = '2018-04-15 vs. 2019-04-18', shape = aspectgroup)) +
  geom_point(data = cd, aes(x = VI180517, y = VI190513, color = '2018-05-17 vs. 2019-05-13', shape = aspectgroup)) +
  geom_point(data = cd, aes(x = VI180611, y = VI190605, color = '2018-06-11 vs. 2019-06-05', shape = aspectgroup)) +
  geom_point(data = cd, aes(x = VI180704, y = VI190716, color = '2018-07-04 vs. 2019-07-16', shape = aspectgroup)) +
  labs(x = '2018', y = '2019', color = 'time pairs') + ggtitle('SAVI in 2018 and 2019') + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

cd_ <- cd[which(cd$aspectgroup %in% c('S', 'N')),]

ggplot() + geom_point(data = cd_, aes(x = VI180307, y = VI190310, color = '2018-03-07 vs. 2019-03-10', shape = aspectgroup)) + 
  geom_point(data = cd_, aes(x = VI180415, y = VI190418, color = '2018-04-15 vs. 2019-04-18', shape = aspectgroup)) +
  geom_point(data = cd_, aes(x = VI180517, y = VI190513, color = '2018-05-17 vs. 2019-05-13', shape = aspectgroup)) +
  geom_point(data = cd_, aes(x = VI180611, y = VI190605, color = '2018-06-11 vs. 2019-06-05', shape = aspectgroup)) +
  geom_point(data = cd_, aes(x = VI180704, y = VI190716, color = '2018-07-04 vs. 2019-07-16', shape = aspectgroup)) +
  labs(x = '2018', y = '2019', color = 'time pairs', shape = 'aspect group') + 
  ggtitle('SAVI in 2018 and 2019') + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggplot() + geom_point(data = cd_, aes(x = VI180307, y = VI190418, color = '2018-03-07 vs. 2019-04-18', shape = aspectgroup)) + 
  geom_point(data = cd_, aes(x = VI180415, y = VI190513, color = '2018-04-15 vs. 2019-05-13', shape = aspectgroup)) +
  geom_point(data = cd_, aes(x = VI180517, y = VI190605, color = '2018-05-17 vs. 2019-06-05', shape = aspectgroup)) +
  geom_point(data = cd_, aes(x = VI180611, y = VI190716, color = '2018-06-11 vs. 2019-07-16', shape = aspectgroup)) +
  xlim(-0.2, 0.75) + ylim(-0.2, 0.75) + 
  geom_abline(slope = 1, intercept = 0) +
  labs(x = '2018', y = '2019', color = 'time pairs', shape = 'aspect group') + ggtitle('SAVI in 2018 and 2019, time lagged') + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

## For ET
et_18 <- cbind(a307_18[, c('x', 'y', 'slopegroup', 'aspectgroup', 'ET')], 
               a415_18[, 'ET'], a517_18[, 'ET'], a611_18[, 'ET'],  a704_18[, 'ET'])
names(et_18) <- c('x', 'y', 'slopegroup', 'aspectgroup', 'VI180307', 'VI180415',
                  'VI180517', 'VI180611', 'VI180704')
et_19 <- cbind(a310_19[, c('x', 'y', 'slopegroup', 'aspectgroup', 'ET')], 
               a418_18[, 'ET'], a513_19[, 'ET'], a605_19[, 'ET'], a716_19[, 'ET'])
names(et_19) <- c('x', 'y', 'slopegroup', 'aspectgroup', 'VI190310', 'VI190418', 'VI190513',
                  'VI190605', 'VI190716')
et_ <- cbind(et_18, et_19[, c(-1, -2, -3, -4)])
ggplot() + geom_point(data = et_, aes(x = VI180307, y = VI190418, color = '2018-03-07 vs. 2019-04-18', shape = aspectgroup)) + 
  geom_point(data = et_, aes(x = VI180415, y = VI190513, color = '2018-04-15 vs. 2019-05-13', shape = aspectgroup)) +
  geom_point(data = et_, aes(x = VI180517, y = VI190605, color = '2018-05-17 vs. 2019-06-05', shape = aspectgroup)) +
  geom_point(data = et_, aes(x = VI180611, y = VI190716, color = '2018-06-11 vs. 2019-07-16', shape = aspectgroup)) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = '2018', y = '2019', color = 'aspectgroup', shape = 'pairs') + ggtitle('ET in 2018 and 2019, time lagged') + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))


ggplot() + geom_boxplot(data = cd_, aes(x = '2019-07-16 - 2019-06-05', y = (VI190716 - VI190605)/41, fill = aspectgroup)) +
  geom_boxplot(data = cd_, aes(x = '2018-06-11 - 2018-05-17', y = (VI180611 - VI180517)/25, fill = aspectgroup)) +
  labs(x = 'Year', y = expression(paste(delta, 'SAVI'))) +
  ggtitle('SAVI increase rate') + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))




bg_170608 <- df[which(df$Date == as.Date('2017-06-08')),]
sm_170608 <- df_small[which(df_small$Date == as.Date('2017-06-08')),]
bg_180611 <- df[which(df$Date == as.Date('2018-06-11')),]
sm_180611 <- df_small[which(df_small$Date == as.Date('2018-06-11')),]


sm_170608NS <- sm_170608[which(sm_170608$aspectgroup == 'N' | sm_170608$aspectgroup == 'S'),]
bg_170608NS <- bg_170608[which(bg_170608$aspectgroup == 'N' | bg_170608$aspectgroup == 'S'),]
sm_180611NS <- sm_180611[which(sm_180611$aspectgroup == 'N' | sm_180611$aspectgroup == 'S'),]
bg_180611NS <- bg_180611[which(bg_180611$aspectgroup == 'N' | bg_180611$aspectgroup == 'S'),]


### Topography features
## small domain
ggplot() + geom_raster(data = sm_170608, aes(x = x, y = y, fill = elevation), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('Elevation') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

ggplot() + geom_raster(data = sm_170608, aes(x = x, y = y, fill = slope), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('Slope') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

ggplot() + geom_raster(data = sm_170608, aes(x = x, y = y, fill = aspect), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('Aspect') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

## big domain
ggplot() + geom_raster(data = bg_170608, aes(x = x, y = y, fill = elevation), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('Elevation') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

ggplot() + geom_raster(data = bg_170608, aes(x = x, y = y, fill = slope), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('Slope') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

ggplot() + geom_raster(data = bg_170608, aes(x = x, y = y, fill = aspect), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('Aspect') + scale_fill_gradientn(colours = pal) + 
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

##do we observe differences in statistical distributions
ggplot() + geom_density(data = sm_170608, aes(x = elevation, color = 'PH')) +
  geom_density(data = bg_170608, aes(x = elevation, color = 'Large Domain')) +
  ggtitle('Elevation distribution') + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

ggplot() + geom_density(data = sm_170608, aes(x = slope, color = 'PH')) +
  geom_density(data = bg_170608, aes(x = slope, color = 'Large Domain')) +
  ggtitle('Slope distribution') + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

ggplot() + geom_density(data = sm_170608, aes(x = aspect, color = 'PH')) +
  geom_density(data = bg_170608, aes(x = aspect, color = 'Large Domain')) +
  ggtitle('Aspect distribution') + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

## Energy balance
ggplot() + geom_raster(data = sm_170608, aes(x = x, y = y, fill = Rn), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('Net radiation, 2017-06-08') + scale_fill_gradientn(colours = pal) +
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
ggplot() + geom_raster(data = sm_180611, aes(x = x, y = y, fill = Rn), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('Net radiation, 2018-06-11') + scale_fill_gradientn(colours = pal) +
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggplot() + geom_raster(data = bg_170608, aes(x = x, y = y, fill = Rn), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('Net radiation, 2017-06-08') + scale_fill_gradientn(colours = pal) +
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
ggplot() + geom_raster(data = bg_180611, aes(x = x, y = y, fill = Rn), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('Net radiation, 2018-06-11') + scale_fill_gradientn(colours = pal) +
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggplot() + geom_density(data = sm_170608, aes(x = Rn, color = 'PH, 2017-06-08')) +
  geom_density(data = bg_170608, aes(x = Rn, color = 'L, 2017-06-08')) +
  geom_density(data = sm_180611, aes(x = Rn, color = 'PH, 2018-06-11'))+
  geom_density(data = bg_180611, aes(x = Rn, color = 'L, 2018-06-11')) +
  ggtitle('Rn distribution') + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

ggplot() + geom_boxplot(data = sm_170608NS, aes(x = slopegroup, y = Rn, fill = aspectgroup)) +
  ggtitle('PH: Rn 2017-06-08') + xlab('Slopegroup') + ylab('Rn') +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggplot() + geom_boxplot(data = bg_170608NS, aes(x = slopegroup, y = Rn, fill = aspectgroup)) +
  ggtitle('L: Rn 2017-06-08') + xlab('Slopegroup') + ylab('Rn') +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))


ggplot() + geom_boxplot(data = sm_180611NS, aes(x = slopegroup, y = Rn, fill = aspectgroup)) +
  ggtitle('PH: Rn 2018-06-11') + xlab('Slopegroup') + ylab('Rn') +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggplot() + geom_boxplot(data = bg_180611NS, aes(x = slopegroup, y = Rn, fill = aspectgroup)) +
  ggtitle('L: Rn 2018-06-11') + xlab('Slopegroup') + ylab('Rn') +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))


## Vegetation index
ggplot() + geom_raster(data = sm_170608, aes(x = x, y = y, fill = NDVI), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('NDVI, 2017-06-08') + scale_fill_gradientn(colours = pal) +
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggplot() + geom_raster(data = sm_180611, aes(x = x, y = y, fill = NDVI), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('NDVI, 2018-06-11') + scale_fill_gradientn(colours = pal) +
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggplot() + geom_raster(data = bg_170608, aes(x = x, y = y, fill = NDVI), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('NDVI, 2017-06-08') + scale_fill_gradientn(colours = pal) +
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggplot() + geom_raster(data = bg_180611, aes(x = x, y = y, fill = NDVI), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('NDVI, 2018-06-11') + scale_fill_gradientn(colours = pal) +
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))


ggplot() + geom_boxplot(data = sm_180611NS, aes(x = slopegroup, y = NDVI, fill = aspectgroup)) +
  ggtitle('PH: NDVI 2018-06-11') + xlab('Slopegroup') + ylab('NDVI') +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggplot() + geom_boxplot(data = sm_170608NS, aes(x = slopegroup, y = NDVI, fill = aspectgroup)) +
  ggtitle('PH: NDVI 2017-06-08') + xlab('Slopegroup') + ylab('NDVI') +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))


ggplot() + geom_boxplot(data = bg_180611NS, aes(x = slopegroup, y = NDVI, fill = aspectgroup)) +
  ggtitle('L: NDVI 2018-06-11') + xlab('Slopegroup') + ylab('NDVI') +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggplot() + geom_boxplot(data = bg_170608NS, aes(x = slopegroup, y = NDVI, fill = aspectgroup)) +
  ggtitle('L: NDVI 2017-06-08') + xlab('Slopegroup') + ylab('NDVI') +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))


## ET and components
ggplot() + geom_raster(data = sm_170608, aes(x = x, y = y, fill = ET), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('ET, 2017-06-08') + scale_fill_gradientn(colours = pal) +
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggplot() + geom_raster(data = sm_180611, aes(x = x, y = y, fill = ET), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('ET, 2018-06-11') + scale_fill_gradientn(colours = pal) +
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

## ETs
ggplot() + geom_raster(data = sm_170608, aes(x = x, y = y, fill = ETs), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('Evaporation, 2017-06-08') + scale_fill_gradientn(colours = pal) +
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggplot() + geom_raster(data = sm_180611, aes(x = x, y = y, fill = ETs), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('Evaporation, 2018-06-11') + scale_fill_gradientn(colours = pal) +
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))


## ETc
ggplot() + geom_raster(data = sm_170608, aes(x = x, y = y, fill = ETc), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('Transpiration, 2017-06-08') + scale_fill_gradientn(colours = pal) +
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggplot() + geom_raster(data = sm_180611, aes(x = x, y = y, fill = ETc), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('Transpiration, 2018-06-11') + scale_fill_gradientn(colours = pal) +
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

sm_170608NS <- sm_170608[which(sm_170608$aspectgroup == 'N' | sm_170608$aspectgroup == 'S'),]
bg_170608NS <- bg_170608[which(bg_170608$aspectgroup == 'N' | bg_170608$aspectgroup == 'S'),]
sm_180611NS <- sm_180611[which(sm_180611$aspectgroup == 'N' | sm_180611$aspectgroup == 'S'),]
bg_180611NS <- bg_180611[which(bg_180611$aspectgroup == 'N' | bg_180611$aspectgroup == 'S'),]


ggplot() + geom_boxplot(data = sm_170608NS, aes(x = slopegroup, y = ET, fill = aspectgroup)) +
  ggtitle('PH: ET 2017-06-08') + xlab('Slopegroup') + ylab('ET') +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
ggplot() + geom_boxplot(data = sm_180611NS, aes(x = slopegroup, y = ET, fill = aspectgroup)) +
  ggtitle('PH: ET 2018-06-11') + xlab('Slopegroup') + ylab('ET') +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggplot() + geom_boxplot(data = sm_170608NS, aes(x = slopegroup, y = ET, fill = aspectgroup)) +
  ggtitle('PH: ET 2017-06-08') + xlab('Slopegroup') + ylab('ET') +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))


## Checking the troublesome 2019-04-18
sm_190418 <- df_small[which(df_small$Date == as.Date('2019-04-18')),]
ggplot() + geom_raster(data = sm_190418, aes(x = x, y = y, fill = ET), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('ET, 2019-04-18') + scale_fill_gradientn(colours = pal) +
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggplot() + geom_raster(data = sm_190418, aes(x = x, y = y, fill = NDVI), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('NDVI, 2019-04-18') + scale_fill_gradientn(colours = pal) +
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggplot() + geom_raster(data = sm_190418, aes(x = x, y = y, fill = ETc), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('ETc, 2019-04-18') + scale_fill_gradientn(colours = pal) +
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggplot() + geom_raster(data = sm_190418, aes(x = x, y = y, fill = ETs), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('ETs, 2019-04-18') + scale_fill_gradientn(colours = pal) +
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggplot() + geom_raster(data = sm_190418, aes(x = x, y = y, fill = VP), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('ETs, 2019-04-18') + scale_fill_gradientn(colours = pal) +
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggplot() + geom_raster(data = sm_190418, aes(x = x, y = y, fill = RH^((SVP-VP)/1000)), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('ETs, 2019-04-18') + scale_fill_gradientn(colours = pal) +
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

ggplot() + geom_raster(data = sm_190418, aes(x = x, y = y, fill = Tmean), interpolate = TRUE) + 
  coord_quickmap() + ggtitle('ETs, 2019-04-18') + scale_fill_gradientn(colours = pal) +
  geom_point(data = pos, aes(x = lon, y = lat)) + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))


## savi/ndvi changes over time (rate of change)







