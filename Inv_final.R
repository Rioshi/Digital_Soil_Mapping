library(sp)
library(raster)
library(e1071)
library(caret)
library(dplyr)
library(ggplot2)
library(rpart)
library(rattle)
########Lectura de datos
cali <- shapefile('Shape/Suelo_cal.shp')
cali
summary(cali)
cali2 <- as.data.frame(cali)
cali2 <- cali2[,-c(2:3,10,12)]
cali2$X <- as.numeric(cali2$X)
cali2$Y <- as.numeric(cali2$Y)
cali2 <- cali2[,1:12]


#Resumen por taxon
cali2 %>% 
  group_by(cali2$ST) %>% 
  summarise(Promedio=mean(slope),sd=sd(slope),max=max(slope),min=min(slope))

cali2 %>% 
  group_by(cali2$ST) %>% 
  summarise(Promedio=mean(MDE),sd=sd(MDE),max=max(MDE),min=min(MDE))

#Bosplot por taxon
box <- ggplot(cali2, aes(x=ST, y=slope)) + 
  geom_boxplot() + theme_classic() + theme(legend.position="bottom") +
  xlab("") + ylab("Pendiente (%)")
box

box <- ggplot(cali2, aes(x=ST, y=MDE)) + 
  geom_boxplot() + theme_classic() + theme(legend.position="bottom") +
  xlab("") + ylab("Elevación (m)")
box

box <- ggplot(cali2, aes(x=ST, y=Aspect)) + 
  geom_boxplot() + theme_classic() + theme(legend.position="bottom") +
  xlab("") + ylab("Orientación (°)")
box

box <- ggplot(cali2, aes(x=ST, y=Aspect)) + 
  geom_boxplot() + theme_classic() + theme(legend.position="bottom") +
  xlab("") + ylab("Orientación (°)")
box

box <- ggplot(cali2, aes(x=ST, y=SAVI_mean)) + 
  geom_boxplot() + theme_classic() + theme(legend.position="bottom") +
  xlab("") + ylab("SAVI")
box

box <- ggplot(cali2, aes(x=ST, y=CI_mean)) + 
  geom_boxplot() + theme_classic() + theme(legend.position="bottom") +
  xlab("") + ylab("Costra Biológica")
box



box <- ggplot(cali2, aes(x=ST, y=FER_mean)) + 
  geom_boxplot() + theme_classic() + theme(legend.position="bottom") +
  xlab("") + ylab("FER")
box

box <- ggplot(cali2, aes(x=ST, y=CARB_mean)) + 
  geom_boxplot() + theme_classic() + theme(legend.position="bottom") +
  xlab("") + ylab("CARB")
box

box <- ggplot(cali2, aes(x=ST, y=ARCI_mean)) + 
  geom_boxplot() + theme_classic() + theme(legend.position="bottom") +
  xlab("") + ylab("ARC")
box

###Modelling
covar <- read.table("Shape/covars.txt",header = TRUE)
model <- naiveBayes(ST~.,data=cali2[,c(4,5:12)])
model$tables
covar$TAXA<-predict(model,newdata=covar[,3:10],type="class")
covar$TAXA2 <-as.numeric(covar$TAXA)
write.table(covar,"covars/classify.txt", sep="\t", row.names=F)

##Modelling Rtree

model <- rpart(ST~.,data=cali2[,c(4,5:12)])
covar$TAXA<-predict(model,newdata=covar[,3:10],type ="class")
covar$TAXA2 <-as.numeric(covar$TAXA)

drawTreeNodes(model,cex=.6,pch=11,size=4*.8, col=NULL,nodeinfo=TRUE,
              units="",cases="obs",digits=getOption("digits"),decimals=2,print.levels=TRUE,
              new=TRUE)

write.table(covar,"covars/classify.txt", sep="\t", row.names=F)
levels(covar$TAXA)
