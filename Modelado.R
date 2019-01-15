#####################
##Lectura de datos##
#####################
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
rasterVis::levelplot(covariables, layers = 10, contour=FALSE,margin=FALSE,
                     par.settings = RdBuTheme, xlim=c(324705,340000),ylim=c(8682745,8695000),
                     colorkey=list(
                       space='bottom',                   
                       labels=list(font=4,cex=1.5),
                       axis.line=list(col='black'),
                       width=1
                     ))

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

hist(na.omit(clima$wind),maxpixels=1500000,xlim=c(2.2,3.6),main="Velocidad del viento",col="gray",
     prob=TRUE,ylim=c(0,4),xlab="ks:  D = 0.114   p-valor = 0.720",ylab="Densidad",breaks=50, cex.lab=1.3,cex.main=1.3)
lines(density(cali$wind,bw = "sj"),col="red",lwd = 3)

hist(na.omit(clima$Hd),maxpixels=1500000,xlim=c(0.5,1.6),main="Presión devapor de agua",col="gray",
     prob=TRUE,ylim=c(0,3),xlab="ks:  D = 0.103   p-valor = 0.820",ylab="Densidad",breaks=50, cex.lab=1.3,cex.main=1.3)
lines(density(cali$Hd,bw = "sj"),col="red",lwd = 3)

hist(na.omit(clima$Rd),maxpixels=1500000,xlim=c(16200,18500),main="Radiación solar total anual",col="gray",
     prob=TRUE,ylim=c(0,0.003),xlab="ks:  D = 0.243   p-valor = 0.023",ylab="Densidad",breaks=50, cex.lab=1.3,cex.main=1.3)
lines(density(cali$Rd,bw = "sj"),col="red",lwd = 3)

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
######################################
library(caret)
library("klaR")
library(e1071)


levels(cali$ST) <- make.names(levels(factor(cali$ST)))

#Balanceo de clases
set.seed(9560)
up_train <- upSample(x = cali[,5:25],
                     y = cali$ST)
#Data splitting
set.seed(3456)
trainIndex <- createDataPartition(up_train$Class, p = .7, 
                                  list = FALSE, 
                                  times = 1)
Train <- up_train[ trainIndex,]
Test  <- up_train[-trainIndex,]

#Validacion Cruzada
set.seed(25)
setControl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 100,
  verboseIter = FALSE,
  classProbs=TRUE # ,
# sampling = "up" 
)

#Naive Bayes Parameters
grid <- data.frame(fL=1000, usekernel=TRUE,adjust=1)
set.seed(11)
#Naive Bayes Model
model <- train(Class~.,data=Train,'nb',
               metric=c("Accuracy","Kappa"),
               tuneGrid=grid,
               preProcess=c("center", "scale", "YeoJohnson"),
               trControl=setControl,
               na.action = na.omit
)

cm1 <- confusionMatrix(data=predict(model,Train),reference=Train$Class)

#Matriz de confusion
pred.test <- predict(model,Test)
confusionMatrix(data=pred.test,reference=Test$Class)

###############################
#Importancia de la variable
##############################
vi <- varImp(model,scale=TRUE,sort=TRUE)
colnames(vi$importance) <- c("Aridic Haplustolls", "Pachic Haplustolls","Sodic Haplocambids",
                             "Torriorthentic Haplustolls", "Typic Haplocambids", "Typic Torriorthents")
rownames(vi$importance) <- c("ET","PPT","TM","WS","WVP","RS","QI","CI","MI","NDVI","CB","OR","ICO",
                             "CPP","CPL","LS","MDE","RSP","PD","TWI","VD")
plot(vi)
ggplot(vi) + xlab("Covariables - Factores de formación de suelos") +
  ylab("Importancia de la Variable") +
  theme(axis.text = element_text(size=11.5),
        strip.text = element_text(size = 12, face="bold"), 
        axis.title.x = element_text(size=12, face="bold"),
        axis.title.y = element_text(size=12, face="bold"))

#para mas info https://github.com/topepo/caret/blob/master/pkg/caret/R/plot.varImp.train.R
plotObj <- sortImp(vi, dim(vi$importance)[1])


write.csv(vi$importance,file="importancivariable.csv")

#boxplot variables comunes a los7 taxones
p1 <- ggplot(cali, aes(x=ST,y=Topographic_Wetness_Index))+
  geom_boxplot(fill="gray")  + scale_fill_grey()  +
  theme(axis.text.y=element_text(face="bold"))+
  xlab("") + ylab("Índice Topográfico de Humedad")

p2 <- ggplot(cali, aes(x=ST,y=LS_Factor))+
  geom_boxplot(fill="gray") + coord_flip() + scale_fill_grey()  +
  theme(axis.text.y = element_blank())+
  xlab("") + ylab("Factor LS")

p3 <- ggplot(cali, aes(x=ST,y=mafic))+
  geom_boxplot(fill="gray") + coord_flip() + scale_fill_grey()  + 
  theme(axis.text.y = element_blank())+
  xlab("") + ylab("Índice Máfico")

library("gridExtra")
grid.arrange(p1,p2,p3,ncol=3)
grid.arrange(p1,                             # First row with one plot spaning over 2 columns
             arrangeGrob(p2, p3, ncol = 2), # Second row with 2 plots in 2 different columns
             nrow = 2)                       # Number of rows

#Prediction of classes
mm <- predict(object=covariables, model=model, fun=predict, type="raw") #type raw = probability, prob = class
mm@data@attributes

#Prediccion de clases como dataframe
covariables.DF <- as.data.frame(covariables,xy=TRUE,na.rm=TRUE)
mm2 <- predict(object=model,newdata=covariables.DF,type="prob")


rnam <- row.names(mm2)
rnam2 <- row.names(covariables.DF)
mm2["ID"] <- rnam
covariables.DF["ID"] <- rnam2

ST <- merge(covariables.DF,mm2,by="ID")
save(ST, file = "H:/TESIS/2018/ST.RData")
load("H:/TESIS/2018/ST.RData")

AH <- ST[,c(2,3,25)]
coordinates(AH) <- ~ x + y
gridded(AH) <- TRUE
rasterAH <- raster(AH)


PH <- ST[,c(2,3,26)]
coordinates(PH) <- ~ x + y
gridded(PH) <- TRUE
rasterPH <- raster(PH)


SH <- ST[,c(2,3,27)]
coordinates(SH) <- ~ x + y
gridded(SH) <- TRUE
rasterSH <- raster(SH)


ToH <- ST[,c(2,3,28)]
coordinates(ToH) <- ~ x + y
gridded(ToH) <- TRUE
rasterToH <- raster(ToH)


TyH <- ST[,c(2,3,29)]
coordinates(TyH) <- ~ x + y
gridded(TyH) <- TRUE
rasterTyH <- raster(TyH)


TyT <- ST[,c(2,3,30)]
coordinates(TyT) <- ~ x + y
gridded(TyT) <- TRUE
rasterTyT <- raster(TyT)

rm("TyT","TyH","ToH","SH","PH","SH","AH")


plot(rasterAH,main="Aridic Haplustolls")
plot(rasterPH,main="Pachic Haplustolls")
plot(rasterSH,main="Sodic Haplocambids")
plot(rasterToH,main="Torriorthentic Haplustolls")
plot(rasterTyH,main="Typic Haplocambids")
plot(rasterTyT,main="Typic Torriorthents")

Probs <- stack(rasterAH,rasterPH,rasterSH,rasterToH,rasterTyH,rasterTyT)

mm3 <- predict(object = model,newdata=covariables.DF,type = "raw",threshold = 0.5)
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
#Graficar cada taxon de suelo con su importancia de la variable
#################################################################
#Aridic Haplustolls
p1 <- ggplot(data=plotObj,aes(x=reorder(row.names(plotObj), plotObj$`Aridic Haplustolls`),y=plotObj$`Aridic Haplustolls`)) +
  geom_bar(stat="identity") + xlab("") +
  ylab("Importancia de la variable") + theme(axis.text.x=element_text(face="bold",size=12),
                                             axis.title.y=element_text(face="bold",size=12))




p2 <-levelplot(Probs,par.settings = RdBuTheme, layers=1,
          xlab=NULL, ylab=NULL, scales=list(draw=FALSE),
          margin=FALSE, maxplixels =1e40, main = "Aridic Haplustolls"
)
grid.arrange(p2,p1,nrow=2)
#Sodic Haplocambids
p1 <- ggplot(data=plotObj,aes(x=reorder(row.names(plotObj), plotObj$`Sodic Haplocambids`),y=plotObj$`Sodic Haplocambids`)) +
  geom_bar(stat="identity") + xlab("") +
  ylab("Importancia de la variable") + theme(axis.text.x=element_text(face="bold",size=12),
                                             axis.title.y=element_text(face="bold",size=12))
p2 <-levelplot(Probs,par.settings = RdBuTheme, layers=3,
               xlab=NULL, ylab=NULL, scales=list(draw=FALSE),
               margin=FALSE, maxplixels =1e40, main = "Sodic Haplocambids"
)
grid.arrange(p2,p1,nrow=2)
#Pachic Haplustolls
p1 <- ggplot(data=plotObj,aes(x=reorder(row.names(plotObj), plotObj$`Pachic Haplustolls`),y=plotObj$`Pachic Haplustolls`)) +
  geom_bar(stat="identity") + xlab("") +
  ylab("Importancia de la variable") + theme(axis.text.x=element_text(face="bold",size=12),
                                             axis.title.y=element_text(face="bold",size=12))
p2 <-levelplot(Probs,par.settings = RdBuTheme, layers=2,
               xlab=NULL, ylab=NULL, scales=list(draw=FALSE),
               margin=FALSE, maxplixels =1e40, main = "Pachic Haplustolls"
)
grid.arrange(p2,p1,nrow=2)
#Torriorthentic Haplustolls
p1 <- ggplot(data=plotObj,aes(x=reorder(row.names(plotObj), plotObj$`Torriorthentic Haplustolls`),y=plotObj$`Torriorthentic Haplustolls`)) +
  geom_bar(stat="identity") + xlab("") +
  ylab("Importancia de la variable") + theme(axis.text.x=element_text(face="bold",size=12),
                                             axis.title.y=element_text(face="bold",size=12))
p2 <-levelplot(Probs,par.settings = RdBuTheme, layers=4,
               xlab=NULL, ylab=NULL, scales=list(draw=FALSE),
               margin=FALSE, maxplixels =1e40, main = "Torriorthentic Haplustolls"
)
grid.arrange(p2,p1,nrow=2)
#Typic Haplocambids
p1 <- ggplot(data=plotObj,aes(x=reorder(row.names(plotObj), plotObj$`Typic Haplocambids`),y=plotObj$`Typic Haplocambids`)) +
  geom_bar(stat="identity") + xlab("") +
  ylab("Importancia de la variable") + theme(axis.text.x=element_text(face="bold",size=12),
                                             axis.title.y=element_text(face="bold",size=12))
p2 <-levelplot(Probs,par.settings = RdBuTheme, layers=5,
               xlab=NULL, ylab=NULL, scales=list(draw=FALSE),
               margin=FALSE, maxplixels =1e40, main = "Typic Haplocambids"
)
grid.arrange(p2,p1,nrow=2)
#Typic Torriorthents
p1 <- ggplot(data=plotObj,aes(x=reorder(row.names(plotObj), plotObj$`Typic Torriorthents`),y=plotObj$`Typic Torriorthents`)) +
  geom_bar(stat="identity") + xlab("") +
  ylab("Importancia de la variable") + theme(axis.text.x=element_text(face="bold",size=12),
                                             axis.title.y=element_text(face="bold",size=12))
p2 <-levelplot(Probs,par.settings = RdBuTheme, layers=6,
               xlab=NULL, ylab=NULL, scales=list(draw=FALSE),
               margin=FALSE, maxplixels =1e40, main = "Typic Torriorthents"
)
grid.arrange(p2,p1,nrow=2)



#################################################################
#Extraer datos de unidades homogeneas
uhomo2 <- raster(as.numeric(uhomo$))

datos <- cbind(datos, extract(covariables, datos.sp))


##PLOT rasters
library(rasterVis)
library(RColorBrewer)
colr <- colorRampPalette(brewer.pal(11, 'RdYlBu'))
nam <- c("Aridic Haplustolls","Pachic Haplustolls","Sodic Haplocambids","Torriorthentic Haplustolls",
         "Typic Haplocambids","Typic Torriorthents")


levelplot(Probs, 
          margin=TRUE,   
          layers=c(1,3,4,5),
          colorkey=list(
            space='right',                   
            labels=list(at=seq(0,1,by=0.2), font=4),
            axis.line=list(col='black'),
            width=1
          ),  
          par.settings=list(
            strip.border=list(col='black'),
            strip.background=list(col='transparent'),
            axis.line=list(col='transparent')
          ),
          scales=list(draw=FALSE),  
          par.strip.text=list(
            cex=1.5
          ),
          col.regions=colr,                   
          at=seq(0,1, len=101),
          names.attr=nam[c(1,3,4,5)]) #+           
  layer(sp.polygons(uhomo, lwd=0.5))

#Plot categorical data in levelplot#
r <- ratify(rasterclass.ST)
rat <- levels(r)[[1]]
rat$Soil <- c('Aridic Haplustolls','Pachic Haplustolls','Sodic Haplocambids',
         'Torriorthentic Haplustolls','Typic Haplocambids','Typic Torriorthents')
levels(r) <- rat
levelplot(r, col.regions=c("#00b144", "#ffff00", "#86b1ec","#0000ff", "#D2691E", "#BB0000"),
          colorkey=list(
            space='right'
          ))
#



levelplot(Probs,par.settings = RdBuTheme, layers=c(2,6),
#          layout=c(3, 2),
          xlab=NULL, ylab=NULL, scales=list(draw=FALSE)
#          names.attr= nam
          )


levelplot(Probs, layers = 1,margin = list(FUN = 'mean'), contour=TRUE,par.settings = magmaTheme)

levelplot(Probs, layers = 1, margin = list(FUN = 'mean'), contour=TRUE,par.settings = magmaTheme)

levelplot(Probs, layers = 5,par.settings = plasmaTheme)
levelplot(Probs, layers = 5,par.settings = magmaTheme)
levelplot(Probs, layers = 5,par.settings = rasterTheme)
levelplot(Probs, layers = 5,par.settings = PuOrTheme)
levelplot(Probs, layers = 5,par.settings = RdBuTheme)
levelplot(Probs, layers = 5,par.settings = streamTheme)
levelplot(Probs, layers = 5,par.settings = viridisTheme)
#levelplot(Probs, layers = 5,par.settings = xscale.raster)
#levelplot(Probs, layers = 5,par.settings = xscale.raster.subticks)
levelplot(Probs, layers = 5,par.settings = YlOrRdTheme)
#levelplot(Probs, layers = 5,par.settings = yscale.raster)
#levelplot(Probs, layers = 5,par.settings = yscale.raster.subticks)



###Resultados
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), 
       widths=c(1,1), heights=c(1,1))

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE), 
       widths=c(1,1), heights=c(1,1))

################################################

####Balanceo alternativo
#set.seed(9560)
#down_train <- downSample(x = cali[,5:25],
#                         y = cali$ST)
set.seed(9560)
up_train <- upSample(x = cali[,5:25],
                     y = cali$ST)

#library(DMwR)
#set.seed(111)
#smote_train  <- SMOTE(ST~.,data=cali[,4:25])                         
#table(smote_train$Class) 

#library(ROSE)
#set.seed(9560)
#rose_train <- ROSE(ST~.,data=cali[,4:25])




