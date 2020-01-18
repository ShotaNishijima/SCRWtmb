
devtools::install_github("ShotaNishijima/SCRWtmb")
library(SCRWtmb)

setwd("C:/Users/00007802/Dropbox/switching_state-space_model")

use_scrw_tmb() # activate CPP file

dat = read.csv("ringo.csv", header = T)
# dat = read.csv("momosiro.csv", header = T)

colnames(dat)
lonlat = dat[,c(6,5)]
lonlat = na.omit(lonlat)

library(rgeos)
library(sp)
library(rgdal)
# GPSデータ（緯度経度形式）をまず sp 形式にする

d_longlat<-SpatialPoints(lonlat,proj4string = CRS("+proj=longlat +ellps=WGS84"))
# UTM に変換
d_utm<-spTransform(d_longlat,CRS("+proj=utm +north +zone=53 +ellps=WGS84"))
d_utm = d_utm@coords
d_utm = data.frame(d_utm)
colnames(d_utm) = c("Easting", "Northing")

# install.packages("tidyverse")
library(tidyverse)
lonlat = bind_cols(lonlat, d_utm)
dat2 = full_join(dat, lonlat)
loc = dat2[,c("Easting","Northing")]
loc2 = data.frame(Easting = loc$Easting - loc$Easting[1],Northing = loc$Northing - loc$Northing[1]) #初期値を0に
head(loc2)
y3 = loc2/1000

model0 = scrw(y3,mode = FALSE,obs_error_type = "normal", sigma_key = c(0,1))
model0$AIC
model0$AICc
model0$gamma
model0$alpha
model0$theta
model0$sigma

model1 = scrw(y3,mode = FALSE,obs_error_type = "normal",
              sigma_key = c(0,1), theta_fix = 0)
model1$AIC
model1$AICc
model1$gamma
model1$alpha
model1$theta
model1$rho
model1$sigma

start <- proc.time()
model2 = scrw(y3,mode =TRUE,obs_error_type = "normal",sigma_key = c(0,1),
              gamma_fix = c(1,0),theta_fix = c(0,0))
(time <- proc.time()-start)

model2$AIC
model2$AICc
model2$gamma
model2$alpha
model2$theta
model2$rho
model2$sigma
model2$det_p
model2$rep

?scrw  # help
names(model2)

plot(y3,xlim=range(y3[,1], na.rm=TRUE),ylim=range(y3[,2], na.rm=TRUE),type = "b")
par(new=T)
plot(model2$x,xlim=range(y3[,1], na.rm=TRUE),ylim=range(y3[,2], na.rm=TRUE),type = "l",col = "red")

model3 = scrw(y3,mode =TRUE,obs_error_type = "normal",sigma_key = c(0,1),
              gamma_fix = c(1,0),theta_fix = c(-1,0))
model3$AIC
model3$AICc
model3$gamma
model3$alpha
model3$theta #小さい
model3$rho
model3$sigma
model3$rep
model3$det_p

model4 = scrw(y3,mode =TRUE,obs_error_type = "normal",sigma_key = c(0,0),
              gamma_fix = c(1,0),theta_fix = c(0,0))
model4$AIC
model4$AICc
model4$gamma
model4$alpha
model4$theta #小さい
model4$rho
model4$sigma
model4$det_p
