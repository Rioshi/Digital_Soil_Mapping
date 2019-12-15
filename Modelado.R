#####################
##Lectura de datos##
#####################
library(raster)
library(caret)
setwd("F:/TESIS/2018")

rm(list = ls())

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
HS_acu <- sum(HS)
ppt_acu <- sum(prec)
tmean <- mean(T_mean)
wind.mean <- mean(Wind)
Hd.mean <- mean(Hd)
Rd.mean <- mean(Rd)
clima <- stack(HS_acu,ppt_acu,tmean,wind.mean,Hd.mean,Rd.mean)
nam <-c("HS_acu","ppt_acu","tmean","wind","Hd","Rd")
names(clima) <- nam
rm(files,HS_acu,ppt_acu,tmean,i,HS,prec,
   stack1,startdir,t_max,T_mean,t_min,BH,nam,Hd,Hd.mean,Rd,Rd.mean,Wind,wind.mean)
clima <-projectRaster(clima, crs='+proj=utm +zone=18 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs')


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
cali <- cali[,1:4]
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

###Covariables remuestreo
relieve2 <- resample(relieve,organismos$ndvi, method="bilinear")
clima2 <- resample(clima,organismos$ndvi, method="bilinear")
lito2 <- resample(lito,organismos$ndvi, method="bilinear")
covariables <- stack(organismos,relieve2,clima2,lito2)
saveRDS(covariables, file = "F:/TESIS/2018/RDATA/covariables.rds")

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

#Estadisticas basicas
media <- cellStats(x=covariables,stat='mean',na.rm=TRUE,asSample=FALSE)
sd <- cellStats(x=covariables,stat='sd',na.rm=TRUE,asSample=FALSE)
CV <- sd*100/abs(media)
max <- cellStats(x=covariables,stat='max',na.rm=TRUE,asSample=FALSE)
min <- cellStats(x=covariables,stat='min',na.rm=TRUE,asSample=FALSE)
library(psych)
ske <- skew(as.data.frame(covariables),na.rm=TRUE) #asimetria
sts <- cbind(media,max,min,sd,CV,ske)
write.csv(x=sts,file="estadisticas_basicas.csv")
as.data.frame()
names(covariables) <- c("NDVI","CB","OR","ICO","CPP","CPL","LS","MDE","RSP","PD","TWI","VD","ET",
                        "PPT","TM","WS","WVP","RS","QI","CI","MI")

#Correlacion entre covariables
options("scipen"=100, "digits"=5)
covar.corr <- layerStats(covariables,stat="pearson",na.rm = TRUE,asSample = FALSE)
write.csv(covar.corr$`pearson correlation coefficient`,file="cormatrix.csv")
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
corrplot(covar.corr$`pearson correlation coefficient`, method="pie", col=col(10),  
         type="lower", order="original", number.digits = 2, number.cex = 1.2, tl.cex = 1.2,cl.cex=1.5,
#         addCoef.col = "black",
         tl.col="black", tl.srt=45,
         diag=TRUE 
)

#Ploteado de covariables
r1 <- plot(covariables$QI, xlim=c(324705,340000),ylim=c(8682745,8695000))
r2 <- plot(covariables$CI)
r3 <- plot(covariables$MI)

library(rasterVis) #magmaTheme RdBuTheme viridisTheme
#indices litologicos
rasterVis::levelplot(covariables, layers = c(19:21), contour=FALSE,
          par.settings = magmaTheme, xlim=c(324705,340000),ylim=c(8682745,8695000))
#indices de relieve
rasterVis::levelplot(covariables, layers = c(5,6), contour=FALSE,
                     par.settings = magmaTheme, xlim=c(324705,340000),ylim=c(8682745,8695000))
r2<-rasterVis::levelplot(covariables, layers = c(13,14), contour=FALSE,margin=FALSE,main="?ndice de Cuarzos",
                     par.settings = RdBuTheme, xlim=c(324705,340000),ylim=c(8682745,8695000),
                     colorkey=list(
                       space='bottom',                   
                       labels=list(font=4,cex=1),
                       axis.line=list(col='black'),
                       width=1
                     )) 
gridExtra::grid.arrange(r1,r2,ncol=2)
#######################################################
## DISTRIBUCION DE LOS DATOS
######################################################

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
ks.test(x=scale(na.omit(clima$wind)),y=scale(cali$wind))
ks.test(x=scale(na.omit(clima$Hd)),y=scale(cali$Hd))
ks.test(x=scale(na.omit(clima$Rd)),y=scale(cali$Rd))

ks.test(x=scale(na.omit(lito$qri)),y=scale(cali$qri))
ks.test(x=scale(na.omit(lito$carb)),y=scale(cali$carb))
ks.test(x=scale(na.omit(lito$mafic)),y=scale(cali$mafic))

ks.test(x=scale(na.omit(organismos$ndvi)),y=scale(cali$ndvi))
ks.test(x=scale(na.omit(organismos$ci)),y=scale(cali$ci))


#Histograma / density plots POBLACION vs MUESTRA
#AGREGAR PRUEBA KOLGOMOROV-SMIRNOFF para ver si provienen de la misma distribucion
options("scipen"=100, "digits"=5)
par(mfrow=c(2,2))

jpeg("Imagenes/explora1.jpeg", width = 30, height = 20, units = 'cm', res = 600, pointsize = 20)
par(mfrow=c(2,2))
hist(na.omit(relieve$MDE),maxpixels=1500000,xlim=c(0,5000),main="Elevaci?n",col="gray",
     prob=TRUE,ylim=c(0,0.0005),xlab="ks:  D = 0.08   p-valor = 0.8635",ylab="Densidad",breaks=20, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$MDE),col="red",lwd = 3)

hist(na.omit(relieve$Slope),maxpixels=1500000,xlim=c(0,200),main="Pendiente",col="gray",
     prob=TRUE,ylim=c(0,0.02),xlab="ks:  D = 0.07   p-valor = 0.90",ylab="Densidad",breaks=50, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$Slope),col="red",lwd = 3)

hist(na.omit(relieve$Aspect),maxpixels=1500000,xlim=c(0,380),main="Orientaci?n",col="gray",
     prob=TRUE,ylim=c(0,0.006),xlab="ks:  D = 0.07   p-valor = 0.92",ylab="Densidad", cex.lab=2, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$Aspect),col="red",lwd = 3)

hist(na.omit(relieve$Convergence_Index),maxpixels=1500000,xlim=c(-20,20),main="?ndice de Convergencia",col="gray",
     prob=TRUE,ylim=c(0,0.15),xlab="ks:  D = 0.23   p-valor < 0.01",ylab="Densidad",breaks=100, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$Convergence_Index,bw = "sj"),col="red",lwd = 3)
dev.off()

jpeg("Imagenes/explora2.jpeg", width = 30, height = 20, units = 'cm', res = 600, pointsize = 20)
par(mfrow=c(2,2))
hist(na.omit(relieve$Cross.Sectional_Curvature),maxpixels=1500000,xlim=c(-0.1,0.1),main="Curvatura plana",col="gray",
     prob=TRUE,ylim=c(0,30),xlab="ks:  D = 0.10   p-valor = 0.51",ylab="Densidad",breaks=100, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$Cross.Sectional_Curvature,bw = "sj"),col="red",lwd = 3)

hist(na.omit(relieve$Longitudinal_Curvature),maxpixels=1500000,xlim=c(-0.15,0.15),main="Curvatura de perfil",col="gray",
     prob=TRUE,ylim=c(0,25),xlab="ks:  D = 0.09   p-valor = 0.71",ylab="Densidad",breaks=100, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$Longitudinal_Curvature,bw = "sj"),col="red",lwd = 3)

hist(na.omit(relieve$LS_Factor),maxpixels=1500000,xlim=c(0,40),main="Factor LS",col="gray",
     prob=TRUE,ylim=c(0,0.10),xlab="ks:  D = 0.08   p-valor = 0.78",ylab="Densidad",breaks=100, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$LS_Factor,bw = "sj"),col="red",lwd = 3)

hist(na.omit(relieve$Relative_Slope_Position),maxpixels=1500000,xlim=c(0,1),main="Posici?n relativa a la pendiente",col="gray",
     prob=TRUE,ylim=c(0,10),xlab="ks:  D = 0.36   p-valor < 0.01",ylab="Densidad",breaks=50, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$Relative_Slope_Position,bw = "sj"),col="red",lwd = 3)
dev.off()

jpeg("Imagenes/explora3.jpeg", width = 30, height = 20, units = 'cm', res = 600, pointsize = 20)
par(mfrow=c(2,2))
hist(na.omit(relieve$Topographic_Wetness_Index),maxpixels=1500000,xlim=c(0,20),main="?ndice topogr?fico de humedad",col="gray",
     prob=TRUE,ylim=c(0,0.4),xlab="ks:  D = 0.15   p-valor = 0.14",ylab="Densidad",breaks=50, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$Topographic_Wetness_Index,bw = "sj"),col="red",lwd = 3)

hist(na.omit(relieve$Valley_Depth),maxpixels=1500000,xlim=c(0,700),main="Profundidad de los valles",col="gray",
     prob=TRUE,ylim=c(0,0.008),xlab="ks:  D = 0.16   p-valor = 0.09",ylab="Densidad",breaks=50, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$Valley_Depth,bw = "sj"),col="red",lwd = 3)

###########
hist(na.omit(clima$tmean),maxpixels=1500000,xlim=c(5,18),main="Temperatura media anual",col="gray",
     prob=TRUE,ylim=c(0,0.25),xlab="ks:  D = 0.08   p-valor = 0.97",ylab="Densidad",breaks=50, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$tmean,bw = "sj"),col="red",lwd = 3)

hist(na.omit(clima$ppt_acu),maxpixels=1500000,xlim=c(100,700),main="Precipitaci?n anual acumulada",col="gray",
     prob=TRUE,ylim=c(0,0.01),xlab="ks:  D = 0.28   p-valor < 0.01",ylab="Densidad",breaks=50, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$ppt_acu,bw = "sj"),col="red",lwd = 3)
dev.off()

jpeg("Imagenes/explora4.jpeg", width = 30, height = 20, units = 'cm', res = 600, pointsize = 20)
par(mfrow=c(2,2))
hist(na.omit(clima$HS_acu),maxpixels=1500000,xlim=c(1050,1500),main="Evapotranspiraci?n anual acumulada",col="gray",
     prob=TRUE,ylim=c(0,0.01),xlab="ks:  D = 0.13   p-valor = 0.58",ylab="Densidad",breaks=50, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$HS_acu,bw = "sj"),col="red",lwd = 3)

hist(na.omit(clima$wind),maxpixels=1500000,xlim=c(2.2,3.6),main="Velocidad del viento",col="gray",
     prob=TRUE,ylim=c(0,4),xlab="ks:  D = 0.11   p-valor = 0.72",ylab="Densidad",breaks=50, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$wind,bw = "sj"),col="red",lwd = 3)

hist(na.omit(clima$Hd),maxpixels=1500000,xlim=c(0.5,1.6),main="Presi?n devapor de agua",col="gray",
     prob=TRUE,ylim=c(0,3),xlab="ks:  D = 0.10   p-valor = 0.82",ylab="Densidad",breaks=50, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$Hd,bw = "sj"),col="red",lwd = 3)

hist(na.omit(clima$Rd),maxpixels=1500000,xlim=c(16200,18500),main="Radiaci?n solar total anual",col="gray",
     prob=TRUE,ylim=c(0,0.003),xlab="ks:  D = 0.24   p-valor = 0.02",ylab="Densidad",breaks=50, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$Rd,bw = "sj"),col="red",lwd = 3)
dev.off()
##########

jpeg("Imagenes/explora5.jpeg", width = 30, height = 20, units = 'cm', res = 600, pointsize = 20)
par(mfrow=c(2,2))
hist(na.omit(lito$qri),maxpixels=1500000,xlim=c(0.995,1.010),main="?ndice de Cuarzo",col="gray",
     prob=TRUE,ylim=c(0,350),xlab="ks:  D = 0.08   p-valor = 0.76",ylab="Densidad",breaks=100, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$qri,bw = "sj"),col="red",lwd = 3)

hist(na.omit(lito$carb),maxpixels=1500000,xlim=c(1,1.010),main="?ndice de carbonatos",col="gray",
     prob=TRUE,ylim=c(0,300),xlab="ks:  D = 0.09   p-valor = 0.72",ylab="Densidad",breaks=100, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$carb,bw = "sj"),col="red",lwd = 3)

hist(na.omit(lito$mafic),maxpixels=1500000,xlim=c(0.98,1.01),main="?ndice m?fico",col="gray",
     prob=TRUE,ylim=c(0,200),xlab="ks:  D = 0.10   p-valor = 0.54",ylab="Densidad",breaks=100, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$mafic,bw = "sj"),col="red",lwd = 3)
############

hist(na.omit(organismos$ndvi),maxpixels=1500000,xlim=c(0,1),main="?ndice Diferenciado de Vegetaci?n Normalizado",col="gray",
     prob=TRUE,ylim=c(0,3),xlab="ks:  D = 0.15   p-valor = 0.10",ylab="Densidad",breaks=100, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$ndvi,bw = "sj"),col="red",lwd = 3)
dev.off()

jpeg("Imagenes/explora6.jpeg", width = 30, height = 20, units = 'cm', res = 600, pointsize = 20)
par(mfrow=c(1,1))
hist(na.omit(organismos$ci),maxpixels=1500000,xlim=c(0.5,2),main="Costra biol?gica",col="gray",
     prob=TRUE,ylim=c(0,6),xlab="ks:  D = 0.10   p-valor = 0.52",ylab="Densidad",breaks=100, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$ci,bw = "sj"),col="red",lwd = 3)
dev.off()

jpeg("Imagenes/cali_cov.jpeg", width = 30, height = 20, units = 'cm', res = 600, pointsize = 20)
par(mfrow=c(2,2))
hist(na.omit(relieve$Convergence_Index),maxpixels=1500000,xlim=c(-20,20),main="?ndice de Convergencia",col="gray",
     prob=TRUE,ylim=c(0,0.15),xlab="ks:  D = 0.23   p-valor < 0.01",ylab="Densidad",breaks=100, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$Convergence_Index,bw = "sj"),col="red",lwd = 3)
hist(na.omit(relieve$Relative_Slope_Position),maxpixels=1500000,xlim=c(0,1),main="Posici?n relativa a la pendiente",col="gray",
     prob=TRUE,ylim=c(0,10),xlab="ks:  D = 0.36   p-valor < 0.01",ylab="Densidad",breaks=50, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$Relative_Slope_Position,bw = "sj"),col="red",lwd = 3)
hist(na.omit(clima$ppt_acu),maxpixels=1500000,xlim=c(100,700),main="Precipitaci?n anual acumulada",col="gray",
     prob=TRUE,ylim=c(0,0.01),xlab="ks:  D = 0.28   p-valor < 0.01",ylab="Densidad",breaks=50, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$ppt_acu,bw = "sj"),col="red",lwd = 3)
hist(na.omit(clima$Rd),maxpixels=1500000,xlim=c(16200,18500),main="Radiaci?n solar total anual",col="gray",
     prob=TRUE,ylim=c(0,0.003),xlab="ks:  D = 0.24   p-valor = 0.02",ylab="Densidad",breaks=50, cex.lab=1.5,cex.main=1.5,cex.axis=1.0)
lines(density(cali$Rd,bw = "sj"),col="red",lwd = 3)
dev.off()

######################################
#MODELADO#
######################################
library(caret)
library("klaR")
library(e1071)
levels(cali$ST) <- make.names(levels(factor(cali$ST)))
#Balanceo de clases
set.seed(9560)
up_train <- upSample(x = cali[,5:25],
                     y = cali$ST)

###Metodo con Validacion Cruzada
set.seed(25)
setControl <- trainControl(
  method = "repeatedcv",
  number = 10, #usar 64/n = numero de taxas
  repeats=100,
  verboseIter = FALSE,
  classProbs=TRUE # ,
  # sampling = "up" 
)
#Naive Bayes Parameters
grid <- data.frame(fL=100, usekernel=TRUE,adjust=1)
set.seed(11)
#Naive Bayes Model
modelcv <- train(Class~.,data=up_train,'nb',
               metric=c("Accuracy","Kappa"),
               tuneGrid=grid,
               preProcess=c("center", "scale", "YeoJohnson"),
               trControl=setControl,
               na.action = na.omit
)
confusionMatrix(data=predict(modelcv,up_train),reference=up_train$Class)





###Metodo con Data splitting
set.seed(3456)
trainIndex <- createDataPartition(up_train$Class, p = 0.7, 
                                  list = FALSE, 
                                  times = 1)
Train <- up_train[ trainIndex,]
Test  <- up_train[-trainIndex,]
#Validacion Cruzada
set.seed(25)
setControl <- trainControl(
  method = "none")
#Naive Bayes Parameters
grid <- data.frame(fL=100, usekernel=TRUE,adjust=1)
set.seed(11)
#Naive Bayes Model
model <- train(Class~.,data=Train,'nb',
               metric=c("Accuracy","Kappa"),
               tuneGrid=grid,
               preProcess=c("center", "scale", "YeoJohnson"),
               trControl=setControl,
               na.action = na.omit
)
#Matriz de confusion
confusionMatrix(data=predict(model,Test),reference=Test$Class)



###############################
#Importancia de la variable
##############################
vi <- varImp(modelcv,scale=TRUE,sort=TRUE)
colnames(vi$importance) <- c("Fluventic Haplocambids", "Lithic Haplocambids","Lithic Torriorthents",
                             "Sodic Haplocambids", "Typic Haplocambids", "Typic Torriorthents")
rownames(vi$importance) <- c("ET","PPT","TM","WS","WVP","RS","QI","CI","MI","NDVI","CB","OR","ICO",
                             "CPP","CPL","LS","MDE","RSP","PD","TWI","VD")
plot(vi)
jpeg("Imagenes/varimptot.jpeg", width = 15, height = 15, units = 'cm', res = 400)
ggplot(vi) + xlab("Variables Predictoras") +
  ylab("Importancia de la Variable (0 - 100") +
  theme(axis.text = element_text(size=11.5),
        strip.text = element_text(size = 12, face="bold"), 
        axis.title.x = element_text(size=12, face="bold"),
        axis.title.y = element_text(size=12, face="bold"))
dev.off()
#para mas info https://github.com/topepo/caret/blob/master/pkg/caret/R/plot.varImp.train.R
plotObj <- sortImp(vi, dim(vi$importance)[1])


write.csv(vi$importance,file="importancivariable.csv")

#Prediction of classes
mm <- predict(object=covariables, model=modelcv, fun=predict, type="raw") #type raw = probability, prob = class
mm@data@attributes

#Prediccion de clases como dataframe
covariables.DF <- as.data.frame(covariables,xy=TRUE,na.rm=TRUE)
mm2 <- predict(object=modelcv,newdata=covariables.DF,type="prob")

#
rnam <- row.names(mm2)
rnam2 <- row.names(covariables.DF)
mm2["ID"] <- rnam
covariables.DF["ID"] <- rnam2
ST <- merge(covariables.DF,mm2,by="ID")
save(ST, file = "H:/TESIS/2018/ST.RData")
load("H:/TESIS/2018/ST.RData")

FH <- ST[,c(2,3,25)] #fluventic Haplocambids
coordinates(FH) <- ~ x + y
gridded(FH) <- TRUE
rasterFH <- raster(FH)
rm(FH)

LH <- ST[,c(2,3,26)] #Lithic Haplocambids
coordinates(LH) <- ~ x + y
gridded(LH) <- TRUE
rasterLH <- raster(LH)
rm(LH)

LT <- ST[,c(2,3,27)] #Lithic Torriorthents
coordinates(LT) <- ~ x + y
gridded(LT) <- TRUE
rasterLT <- raster(LT)
rm(LT)

SH <- ST[,c(2,3,28)] #Sodic Haplocambids
coordinates(SH) <- ~ x + y
gridded(SH) <- TRUE
rasterSH <- raster(SH)
rm(SH)

TH <- ST[,c(2,3,29)] #Typic Haplocambids
coordinates(TH) <- ~ x + y
gridded(TH) <- TRUE
rasterTH <- raster(TH)
rm(TH)

TT <- ST[,c(2,3,30)] #Typic Torriorthents
coordinates(TT) <- ~ x + y
gridded(TT) <- TRUE
rasterTT <- raster(TT)
rm(TT)

Probs <- stack(rasterFH,rasterLH,rasterLT,rasterSH,rasterTH,rasterTT)
rm(rasterFH,rasterLH,rasterLT,rasterSH,rasterTH,rasterTT)
#


mm3 <- predict(object = modelcv,newdata=covariables.DF,type = "raw",threshold = 0.5)
covariables.DF["ST"] <- as.numeric(mm3)
class.ST <- covariables.DF[,c(1,2,25)]
coordinates(class.ST) <- ~ x + y
gridded(class.ST) <- TRUE
rasterclass.ST <- raster(class.ST)
plot(rasterclass.ST)
levels(mm3)
head(mm3,71)
as.numeric(mm3)

writeRaster(x = rasterclass.ST,filename = "RASTER/classrast.tif",overwrite=TRUE)
rasterclass.ST <- raster("RASTER/classrast.tif")

################################################################
#Graficar el mapa conjunto de taxones de suelo
#################################################################
r <- ratify(rasterclass.ST)
rat <- levels(r)[[1]]
#rat$Soil <- c('Fluventic Haplocambids','Lithic Haplocambids','Lithic Torriorthents',
#               'Sodic Haplocambids','Typic Haplocambids','Typic Torriorthents')
rat$Soil <- c('FH','LH','LT',
               'SH','TH','TT')
levels(r) <- rat
jpeg("Imagenes/taxaplot.jpeg", width = 30, height = 25, units = 'cm', res = 400)
rasterVis::levelplot(r, col.regions=c("#00b144", "#ffff00", "#86b1ec","#0000ff", "#D2691E", "#BB0000"),
                     colorkey=list(
                       space='bottom'
                     ))
dev.off()

################################################################
#Graficar cada taxon de suelo con su importancia de la variable
#################################################################
#Cambiar taxon por taxon
library(rasterVis)
library(RColorBrewer)
p1 <- ggplot(data=plotObj,aes(x=reorder(row.names(plotObj), plotObj$`Typic Torriorthents`),y=plotObj$`Typic Torriorthents`)) +
  geom_bar(stat="identity") + xlab("") +
  ylab("Importancia de la variable") + theme(axis.text.x=element_text(face="bold",size=12,angle = 90),
                                             axis.title.y=element_text(face="bold",size=12))
p2 <-rasterVis::levelplot(Probs,par.settings = RdBuTheme, layers=6,
          xlab=NULL, ylab=NULL, scales=list(draw=FALSE),
          margin=FALSE, maxplixels =1e40, main = "Typic Torriorthents"
)
jpeg("Imagenes/TT.jpeg", width = 25, height = 30, units = 'cm', res = 400)
grid.arrange(p2,p1,nrow=2)
dev.off()

################################################################
#Comparacion TAXA vs UHOMO y UCARTO
#################################################################
comp <- read.table("Unidades_homogeneas/INTERSEC_TAXA_UHOMO.txt",header = TRUE,sep=";")
str(comp); comp$classrast <- as.factor(comp$classrast) ; comp$Clave <- as.factor(comp$Clave)
comp1 <- read.table("Unidades_homogeneas/INTERSEC_TAXA_UCARTO.txt",header = TRUE,sep=";")
str(comp1); comp1$classrast <- as.factor(comp1$classrast) ; comp1$Clave <- as.factor(comp1$Clave)
levels(comp1$UCarto)[levels(comp1$UCarto)=="Asociacion Cumb-Chacalla"] <- "Asociacion Salpin-Chacalla"
levels(comp1$SimbUcar)[levels(comp1$SimbUcar)=="Cu-Ch"] <- "Sa-Ch"
levels(comp1$classrast)[levels(comp1$classrast)=="1"] <- "FH"
levels(comp1$classrast)[levels(comp1$classrast)=="2"] <- "LH"
levels(comp1$classrast)[levels(comp1$classrast)=="3"] <- "LT"
levels(comp1$classrast)[levels(comp1$classrast)=="4"] <- "SH"
levels(comp1$classrast)[levels(comp1$classrast)=="5"] <- "TH"
levels(comp1$classrast)[levels(comp1$classrast)=="6"] <- "TT"


nam <- c("Fluventic Haplocambids","Lithic Haplocambids","Lithic Torriorthents","Sodic Haplocambids",
         "Typic Haplocambids","Typic Torriorthents")


#Datos frecuenciales con zonas de vida
table1 <- table(comp$classrast,comp$Abrev_)
prop.table(table1,margin=1)*100
#Datos frecuenciales con geologia
table2 <- table(comp$classrast,comp$NAME)
prop.table(table2,margin=1)*100
#Datos frecuenciales con uhomo
table3 <- table(comp1$classrast,comp$uhomo)
prop.table(table3,margin=1)*100
#Datos frecuenciales con Unidades Cartograficas
table4 <- table(comp1$UCarto,comp1$classrast)
prop.table(table4,margin=1)*100

#Funcion para prueba de Pawlik o Chi2 corregigo
library(DescTools)
paw <- function(x,y){
  tabla <- table(x,y)
  CP <-ContCoef(tabla)
  q <- min(dim(tabla)[1],dim(tabla)[2])
  Cmax <- sqrt((q-1)/q)
  pawlik <- CP/Cmax
  return(pawlik)
}

#Prueba de independencia UHomogenea vs UCartografica
table5 <- table(comp1$UCarto,comp1$uhomo)
chisq.test(table5,simulate.p.value = TRUE)
paw(comp1$UCarto,comp1$uhomo)

#Prueba de independencia UTaxonomica vs UCartografica
table6 <- table(comp1$UCarto,comp1$classrast)
chisq.test(table6,simulate.p.value = TRUE)
paw(comp1$UCarto,comp1$classrast)

#Prueba de independencia UTaxonomica vs UHomogenea
table7 <- table(comp1$uhomo,comp1$classrast)
chisq.test(table7,simulate.p.value = TRUE)
paw(comp1$uhomo,comp1$classrast)

#Analisis de correspondencias simples
library("FactoMineR")
library("factoextra")
library("gplots")

table8 <- round(prop.table(table(comp1$SimbUcar,comp1$classrast))*100,2)
dt <- as.table(as.matrix(table8))
balloonplot(t(dt), main ="housetasks", xlab ="", ylab="",label=FALSE,
            show.margins = FALSE,dotsize=3,scale.method="volume")


res.ca = CA(X=as.data.frame.matrix(table(comp1$SimbUcar,comp1$classrast)))
summary(res.ca)
jpeg("Imagenes/Correspondencias.jpeg", width = 20, height = 20, units = 'cm', res = 400, pointsize = 20)
fviz_ca_biplot(res.ca, repel = TRUE,title="An?lisis de Correspondencias")
dev.off()

res.ca = CA(as.data.frame.matrix(table(comp1$uhomo,comp1$classrast)))
summary(res.ca)
ellipseCA(res.ca,method="boot")
#################################################################
##PLOT rasters
#################################################################
library(rasterVis)
library(RColorBrewer)
colr <- colorRampPalette(brewer.pal(11, 'RdYlBu'))
nam <- c("Fluventic Haplocambids","Lithic Haplocambids","Lithic Torriorthents","Sodic Haplocambids",
         "Typic Haplocambids","Typic Torriorthents")
leg.col <- list(space='right',labels=list(at=seq(0,1,by=0.2), font=4),axis.line=list(col='black'),width=1)

jpeg("Imagenes/Taxa_compilado.jpeg", width = 25, height = 17, units = 'cm', res = 600, pointsize = 600)
rasterVis::levelplot(Probs,layout=c(3,2),col.regions=colr,
                     scales=list(draw=FALSE),names.attr=nam,colorkey=leg.col,par.strip.text=list(
                       cex=1.5
                     ))
dev.off()

jpeg("Imagenes/Vimp_compilado.jpeg", width = 25, height = 17, units = 'cm', res = 600, pointsize = 20)
ggplot(vi) + xlab("Covariables - Factores de formaci?n de suelos") +
  ylab("Importancia de la Variable") +
  theme(axis.text = element_text(size=11.5),
        strip.text = element_text(size = 12, face="bold"), 
        axis.title.x = element_text(size=12, face="bold"),
        axis.title.y = element_text(size=12, face="bold"))
dev.off()