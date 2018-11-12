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
relieve <- dropLayer(relieve, i=c(2,5,6,9,11,12,13,16,18,20))

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
clima <- dropLayer(clima, i=c(1,2,4,5,7,8))


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
cali.sp <- cali
coordinates(cali.sp) = ~X+Y
proj4string(cali.sp) <- CRS("+proj=utm +zone=18 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

#Extraccion de datos de covariables
cali <- cbind(cali, extract(clima, cali.sp))
cali <- cbind(cali, extract(lito, cali.sp))
cali <- cbind(cali, extract(organismos, cali.sp))
cali <- cbind(cali, extract(relieve, cali.sp))

#Corte de covariables en base a la zona de estudio
masc <- rasterize(uhomo,relieve$Aspect)
masc2 <- rasterize(uhomo,clima$HS_acu)
masc3 <- rasterize(uhomo,organismos$ndvi)
masc4 <- rasterize(uhomo,lito$qri)


relieve <- mask(relieve,masc)
clima <- mask(clima,masc2)
organismos <- mask(organismos,masc3)
lito <- mask(lito,masc4)

#######################################################
## ANALISIS EXPLORATORIO DE DATOS
######################################################

#Correlacion entre las covariables
rel.corr <- layerStats(relieve,stat="pearson",na.rm = TRUE) #hacer antes y despues de droplayers
clim.corr <- layerStats(clima,stat="pearson",na.rm = TRUE)
org.corr <- layerStats(organismos,stat="pearson",na.rm = TRUE)
lito.corr <- layerStats(lito,stat="pearson",na.rm = TRUE)

library(agricolae)
abs(correlation(cali[,5:22])$`correlation`)>0.5 #reobtener las covariables luego de droplayers


#Kolgomorov-Smirnoff
ks.test(x=na.omit(relieve$MDE),y=cali$MDE)
ks.test(x=na.omit(relieve$Slope),y=cali$Slope)
ks.test(x=na.omit(relieve$Aspect),y=cali$Aspect)
ks.test(x=na.omit(relieve$Convergence_Index),y=cali$Convergence_Index)
ks.test(x=na.omit(relieve$Cross.Sectional_Curvature),y=cali$Cross.Sectional_Curvature)


#Histograma / density plots POBLACION vs MUESTRA
#AGREGAR PRUEBA KOLGOMOROV-SMIRNOFF para ver si provienen de la misma distribucion
options("scipen"=100, "digits"=5)
hist(na.omit(relieve$MDE),maxpixels=1500000,xlim=c(0,5000),main="Elevación",col="gray",
     prob=TRUE,ylim=c(0,0.0005),xlab="ks:  D = 0.222   p-valor = 0.0036",ylab="Densidad",breaks=20)
lines(density(cali$MDE),col="red",lwd = 3)

hist(na.omit(relieve$Slope),maxpixels=1500000,xlim=c(0,200),main="Pendiente",col="gray",
     prob=TRUE,ylim=c(0,0.02),xlab="ks:  D = 0.347   p-valor < 0.05",ylab="Densidad",breaks=50)
lines(density(cali$Slope),col="red",lwd = 3)

hist(na.omit(relieve$Aspect),maxpixels=1500000,xlim=c(0,380),main="Orientación",col="gray",
     prob=TRUE,ylim=c(0,0.006),xlab="ks:  D = 0.0511   p-valor > 0.05",ylab="Densidad")
lines(density(cali$Aspect),col="red",lwd = 3)

hist(na.omit(relieve$Convergence_Index),maxpixels=1500000,xlim=c(-20,20),main="Orientación",col="gray",
     prob=TRUE,ylim=c(0,0.15),xlab="ks:  D = 0.106   p-valor = 0.47",ylab="Densidad",breaks=100)
lines(density(cali$Convergence_Index,bw = "sj"),col="red",lwd = 3)

hist(na.omit(relieve$Cross.Sectional_Curvature),maxpixels=1500000,xlim=c(-0.1,0.1),main="Orientación",col="gray",
     prob=TRUE,ylim=c(0,30),xlab="ks:  D = 0.124   p-valor = 0.28",ylab="Densidad",breaks=100)
lines(density(cali$Cross.Sectional_Curvature,bw = "sj"),col="red",lwd = 3)
