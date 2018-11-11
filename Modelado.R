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

#Lectura de covariables organismos
load("H:/TESIS/2018/RDATA/organismos.RData")
organismos <- stack(ndvi,ci)
nam <- c("ndvi","ci")
names(organismos) <- nam
rm(ci,nam,ndvi)
