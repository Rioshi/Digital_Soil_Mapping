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

#relieve$Topographic_Wetness_Index <- log2(relieve$Topographic_Wetness_Index)

#Lectura de covariables Clima
load("RDATA/clima.RData")
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
load("RDATA/organismos.RData")
organismos <- stack(ndvi,ci)
nam <- c("ndvi","ci")
names(organismos) <- nam
rm(ci,nam,ndvi)
organismos <-projectRaster(organismos, crs='+proj=utm +zone=18 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs')

#Lectura de covariables litologia
load("RDATA/aster.RData")
qri <- (aster[[11]]*aster[[11]])/(aster[[10]]*aster[[12]]) #quartz index
carb <- aster[[13]]/aster[[12]] #carbonate index
mafic <- (aster[[12]])*(aster[[14]]^3)/(aster[[13]]^4) #mafic index
lito <- stack(qri,carb,mafic)
nam <- c("qri","carb","mafic")
names(lito) <- nam
rm(qri,carb,mafic,aster)
lito <- projectRaster(lito, crs='+proj=utm +zone=18 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs')

#lectura unidades homogeneas de terreno
uhomo <- shapefile("Unidades_homogeneas/uhomo.shp")

#Lectura de calicatas
cali <- read.csv("calic.csv",header = TRUE)
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

#Correlacion entre las covariables por factor
rel.corr <- layerStats(relieve,stat="pearson",na.rm = TRUE) #hacer antes y despues de droplayers
clim.corr <- layerStats(clima,stat="pearson",na.rm = TRUE)
org.corr <- layerStats(organismos,stat="pearson",na.rm = TRUE)
lito.corr <- layerStats(lito,stat="pearson",na.rm = TRUE)

#Prueba variables
nearZeroVar(cali[,5:23],saveMetrics = TRUE)

###Covariables remuestreo
relieve2 <- resample(relieve,organismos$ndvi, method="bilinear")
clima2 <- resample(clima,organismos$ndvi, method="bilinear")
lito2 <- resample(lito,organismos$ndvi, method="bilinear")
covariables <- stack(organismos,relieve2,clima2,lito2)
names(covariables) <- c("NDVI","CI","Or","Icv","Cpa","Cpe","LS","Elv","Psp","Pdt","Ith","Pv","ET","Pt","Tm","Icu","Ica","Ima")

#Estadisticas basicas


media <- cellStats(x=covariables,stat='mean',na.rm=TRUE,asSample=FALSE)
sd <- cellStats(x=covariables,stat='sd',na.rm=TRUE,asSample=FALSE)
CV <- sd*100/abs(media)
library(psych)
ske <- skew(as.data.frame(covariables),na.rm=TRUE) #asimetria
sts <- cbind(media,sd,CV,ske)
write.csv(x=sts,file="estadisticas_basicas.csv")


#Correlacion entre covariables
options("scipen"=100, "digits"=5)
covar.corr <- layerStats(covariables,stat="pearson",na.rm = TRUE,asSample = FALSE)

library(corrplot)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
#col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
#                           "cyan", "#007FFF", "blue","#00007F"))
#col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
#                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
#                           "#4393C3", "#2166AC", "#053061"))
#col3 <- colorRampPalette(c("red", "white", "blue"))
#col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
#                           "cyan", "#007FFF", "blue", "#00007F"))

#Para obtener los p-valor
#covar.matrix <- as.data.frame(covariables,na.rm=TRUE) 
#library(agricolae)
#correl.matrix <- correlation(covar.matrix)

#Guardarlo en 1500 width ,,, heght 1300 
corrplot(covar.corr$`pearson correlation coefficient`, method="color", col=col(10),  
         type="lower", order="original", number.digits = 2, number.cex = 1.2, tl.cex = 1.5,cl.cex=1.5,
         addCoef.col = "black",
         tl.col="black", tl.srt=45,
         diag=FALSE 
)





#Kolgomorov-Smirnoff
ks.test(x=scale(na.omit(relieve$MDE)),y=scale(cali$MDE))
ks.test(x=scale(na.omit(relieve$Slope)),y=scale(cali$Slope))
ks.test(x=scale(na.omit(relieve$Aspect)),y=scale(cali$Aspect))
ks.test(x=scale(na.omit(relieve$Convergence_Index)),y=scale(cali$Convergence_Index))
ks.test(x=scale(na.omit(relieve$Cross.Sectional_Curvature)),y=scale(cali$Cross.Sectional_Curvature))
ks.test(x=scale(na.omit(relieve$Longitudinal_Curvature)),y=scale(cali$Longitudinal_Curvature))
ks.test(x=scale(na.omit(relieve$LS_Factor)),y=scale(cali$LS_Factor))
ks.test(x=scale(na.omit(relieve$Relative_Slope_Position)),y=scale(cali$Relative_Slope_Position))
ks.test(x=scale(na.omit(relieve$Topographic_Wetness_Index)),y=scale(cali$Topographic_Wetness_Index))
ks.test(x=scale(na.omit(relieve$Valley_Depth)),y=scale(cali$Valley_Depth))

ks.test(x=scale(na.omit(clima$HS_acu)),y=scale(cali$HS_acu))
ks.test(x=scale(na.omit(clima$ppt_acu)),y=scale(cali$ppt_acu))
ks.test(x=scale(na.omit(clima$tmean)),y=scale(cali$tmean))

ks.test(x=scale(na.omit(lito$qri)),y=scale(cali$qri))
ks.test(x=scale(na.omit(lito$carb)),y=scale(cali$carb))
ks.test(x=scale(na.omit(lito$mafic)),y=scale(cali$mafic))

ks.test(x=scale(na.omit(organismos$ndvi)),y=scale(cali$ndvi))
ks.test(x=scale(na.omit(organismos$ci)),y=scale(cali$ci))


#Histograma / density plots POBLACION vs MUESTRA
#AGREGAR PRUEBA KOLGOMOROV-SMIRNOFF para ver si provienen de la misma distribucion
options("scipen"=100, "digits"=5)
hist(na.omit(relieve$MDE),maxpixels=1500000,xlim=c(0,5000),main="Elevación",col="gray",
     prob=TRUE,ylim=c(0,0.0005),xlab="ks:  D = 0.222   p-valor = 0.0036",ylab="Densidad",breaks=20, cex.lab=1.3,cex.main=1.3)
lines(density(cali$MDE),col="red",lwd = 3)

hist(na.omit(relieve$Slope),maxpixels=1500000,xlim=c(0,200),main="Pendiente",col="gray",
     prob=TRUE,ylim=c(0,0.02),xlab="ks:  D = 0.347   p-valor < 0.05",ylab="Densidad",breaks=50, cex.lab=1.3,cex.main=1.3)
lines(density(cali$Slope),col="red",lwd = 3)

hist(na.omit(relieve$Aspect),maxpixels=1500000,xlim=c(0,380),main="Orientación",col="gray",
     prob=TRUE,ylim=c(0,0.006),xlab="ks:  D = 0.0511   p-valor > 0.05",ylab="Densidad", cex.lab=1.3,cex.main=1.3)
lines(density(cali$Aspect),col="red",lwd = 3)

hist(na.omit(relieve$Convergence_Index),maxpixels=1500000,xlim=c(-20,20),main="Índice de Convergencia",col="gray",
     prob=TRUE,ylim=c(0,0.15),xlab="ks:  D = 0.106   p-valor = 0.47",ylab="Densidad",breaks=100, cex.lab=1.3,cex.main=1.3)
lines(density(cali$Convergence_Index,bw = "sj"),col="red",lwd = 3)

hist(na.omit(relieve$Cross.Sectional_Curvature),maxpixels=1500000,xlim=c(-0.1,0.1),main="Curvatura plana",col="gray",
     prob=TRUE,ylim=c(0,30),xlab="ks:  D = 0.124   p-valor = 0.28",ylab="Densidad",breaks=100, cex.lab=1.3,cex.main=1.3)
lines(density(cali$Cross.Sectional_Curvature,bw = "sj"),col="red",lwd = 3)

hist(na.omit(relieve$Longitudinal_Curvature),maxpixels=1500000,xlim=c(-0.15,0.15),main="Curvatura de perfil",col="gray",
     prob=TRUE,ylim=c(0,25),xlab="ks:  D = 0.125   p-valor = 0.273",ylab="Densidad",breaks=100, cex.lab=1.3,cex.main=1.3)
lines(density(cali$Longitudinal_Curvature,bw = "sj"),col="red",lwd = 3)

hist(na.omit(relieve$LS_Factor),maxpixels=1500000,xlim=c(0,40),main="Factor LS",col="gray",
     prob=TRUE,ylim=c(0,0.10),xlab="ks:  D = 0.125   p-valor = 0.273",ylab="Densidad",breaks=100, cex.lab=1.3,cex.main=1.3)
lines(density(cali$LS_Factor,bw = "sj"),col="red",lwd = 3)

hist(na.omit(relieve$Relative_Slope_Position),maxpixels=1500000,xlim=c(0,1),main="Posición relativa a la pendiente",col="gray",
     prob=TRUE,ylim=c(0,10),xlab="ks:  D = 0.278   p-valor < 0.05",ylab="Densidad",breaks=50, cex.lab=1.3,cex.main=1.3)
lines(density(cali$Relative_Slope_Position,bw = "sj"),col="red",lwd = 3)

hist(na.omit(relieve$Topographic_Wetness_Index),maxpixels=1500000,xlim=c(0,20),main="Índice topográfico de humedad",col="gray",
     prob=TRUE,ylim=c(0,0.4),xlab="ks:  D = 0.192   p-valor = 0.018",ylab="Densidad",breaks=50, cex.lab=1.3,cex.main=1.3)
lines(density(cali$Topographic_Wetness_Index,bw = "sj"),col="red",lwd = 3)

hist(na.omit(relieve$Valley_Depth),maxpixels=1500000,xlim=c(0,700),main="Profundidad de los valles",col="gray",
     prob=TRUE,ylim=c(0,0.008),xlab="ks:  D = 0.189   p-valor = 0.019",ylab="Densidad",breaks=50, cex.lab=1.3,cex.main=1.3)
lines(density(cali$Valley_Depth,bw = "sj"),col="red",lwd = 3)

###########
hist(na.omit(clima$tmean),maxpixels=1500000,xlim=c(5,18),main="Temperatura media anual",col="gray",
     prob=TRUE,ylim=c(0,0.25),xlab="ks:  D = 0.130   p-valor = 0.191",ylab="Densidad",breaks=50, cex.lab=1.3,cex.main=1.3)
lines(density(cali$tmean,bw = "sj"),col="red",lwd = 3)

hist(na.omit(clima$ppt_acu),maxpixels=1500000,xlim=c(100,700),main="Precipitación anual acumulada",col="gray",
     prob=TRUE,ylim=c(0,0.01),xlab="ks:  D = 0.161   p-valor = 0.282",ylab="Densidad",breaks=50, cex.lab=1.3,cex.main=1.3)
lines(density(cali$ppt_acu,bw = "sj"),col="red",lwd = 3)

hist(na.omit(clima$HS_acu),maxpixels=1500000,xlim=c(1050,1500),main="Evapotranspiración anual acumulada",col="gray",
     prob=TRUE,ylim=c(0,0.01),xlab="ks:  D = 0.213   p-valor = 0.067",ylab="Densidad",breaks=50, cex.lab=1.3,cex.main=1.3)
lines(density(cali$HS_acu,bw = "sj"),col="red",lwd = 3)
##########

hist(na.omit(lito$qri),maxpixels=1500000,xlim=c(0.995,1.010),main="Índice de Cuarzo",col="gray",
     prob=TRUE,ylim=c(0,350),xlab="ks:  D = 0.165   p-valor = 0.061",ylab="Densidad",breaks=100, cex.lab=1.3,cex.main=1.3)
lines(density(cali$qri,bw = "sj"),col="red",lwd = 3)

hist(na.omit(lito$carb),maxpixels=1500000,xlim=c(1,1.010),main="Índice de carbonatos",col="gray",
     prob=TRUE,ylim=c(0,300),xlab="ks:  D = 0.171   p-valor = 0.047",ylab="Densidad",breaks=100, cex.lab=1.3,cex.main=1.3)
lines(density(cali$carb,bw = "sj"),col="red",lwd = 3)

hist(na.omit(lito$mafic),maxpixels=1500000,xlim=c(0.98,1.01),main="Índice máfico",col="gray",
     prob=TRUE,ylim=c(0,200),xlab="ks:  D = 0.152   p-valor = 0.106",ylab="Densidad",breaks=100, cex.lab=1.3,cex.main=1.3)
lines(density(cali$mafic,bw = "sj"),col="red",lwd = 3)
############

hist(na.omit(organismos$ndvi),maxpixels=1500000,xlim=c(0,1),main="Índice Diferenciado de Vegetación Normalizado",col="gray",
     prob=TRUE,ylim=c(0,3),xlab="ks:  D = 0.201   p-valor = 0.012",ylab="Densidad",breaks=100, cex.lab=1.3,cex.main=1.3)
lines(density(cali$ndvi,bw = "sj"),col="red",lwd = 3)

hist(na.omit(organismos$ci),maxpixels=1500000,xlim=c(0.5,2),main="Costra biológica",col="gray",
     prob=TRUE,ylim=c(0,6),xlab="ks:  D = 0.266   p-valor < 0.05",ylab="Densidad",breaks=100, cex.lab=1.3,cex.main=1.3)
lines(density(cali$ci,bw = "sj"),col="red",lwd = 3)


######################################
#MODELADO#
library(caret)
library("klaR")
library(e1071)


levels(cali$ST) <- make.names(levels(factor(cali$ST)))


#Cross Validation

set.seed(25)
setControl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 1,
  verboseIter = FALSE,
  classProbs=TRUE, 
  sampling = "up"
)

#Naive Bayes Parameters
grid <- data.frame(fL=1000, usekernel=TRUE,adjust=1)
#Naive Bayes Model
model <- train(ST~.,data=cali[,4:22],'nb',
               metric=c("Accuracy","Kappa"),
               tuneGrid=grid,
               preProcess=c("center", "scale", "YeoJohnson"),
               trControl=setControl,
               na.action = na.omit
               )

#Prediction of classes
mm <- predict(object=covariables, model=model, fun=predict, type="raw") #type raw = probability, prob = class
mm@data@attributes

mm2 <- stack(predict(object=covariables, model=model, fun=predict.train, type="prob"))

#Importancia de la variable
plot(varImp(model))



###Resultados
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), 
       widths=c(1,1), heights=c(1,1))

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE), 
       widths=c(1,1), heights=c(1,1))




