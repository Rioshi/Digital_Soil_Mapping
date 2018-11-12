##Lectura de datos##
library(raster)
library(caret)
setwd("H:/TESIS/2018")

#Lectura de covariables RELIEVE
startdir <- getwd()
setwd(paste(getwd(), "/RASTER/Relieve/SAGA", sep=""))
files <- list.files(pattern="sdat$")
stack1 <- list()
for(i in 1:length(files)) {
  stack1[[i]] <- raster(files[i])}
relieve <- do.call(stack, stack1) ### JO!
setwd(startdir)
relieve

#Lectura de covariables Clima
load("H:/TESIS/2018/RDATA/clima.RData")
HS_max <- max(HS)
HS_min <- min(HS)
HS_acu <- sum(HS)
ppt_max <- max(prec)
ppt_min <- min(prec)
ppt_acu <- sum(prec)
tmax <- mean(t_max)
tmin <- mean(t_min)
tmean <- mean(T_mean)
clima <- stack(HS_max,HS_min,HS_acu,ppt_max,ppt_min,ppt_acu,tmax,tmin,tmean)
nam <-c("HS_max","HS_min","HS_acu","ppt_max","ppt_min","ppt_acu","tmax","tmin","tmean")
names(clima) <- nam
rm(files,HS_max,HS_min,HS_acu,ppt_max,ppt_min,ppt_acu,tmax,tmin,tmean,i,HS,prec,
   stack1,startdir,t_max,T_mean,t_min,BH,nam)
clima <-projectRaster(clima, crs='+proj=utm +zone=18 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs')

#Lectura de covariables organismos
load("H:/TESIS/2018/RDATA/organismos.RData")
organismos <- stack(ndvi,ci)
nam <- c("ndvi","ci")
names(organismos) <- nam
rm(ci,nam,ndvi)
organismos <-projectRaster(organismos, crs='+proj=utm +zone=18 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs')

#Lectura de covariables litologia
load("H:/TESIS/2018/RDATA/aster.RData")
qri <- (aster[[11]]*aster[[11]])/(aster[[10]]*aster[[12]]) #quartz index
carb <- aster[[13]]/aster[[12]] #carbonate index
mafic <- (aster[[12]])*(aster[[14]]^3)/(aster[[13]]^4) #mafic index
lito <- stack(qri,carb,mafic)
nam <- c("qri","carb","mafic")
names(lito) <- nam
rm(qri,carb,mafic,aster)
lito <- projectRaster(lito, crs='+proj=utm +zone=18 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs')

#lectura unidades homogeneas de terreno
uhomo <- shapefile("H:/TESIS/2018/Unidades_homogeneas/uhomo.shp")

#Lectura de calicatas
cali <- read.csv("H:/TESIS/2018/calic.csv",header = TRUE)
coordinates(cali) = ~X+Y
proj4string(cali) <- CRS("+proj=utm +zone=18 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
