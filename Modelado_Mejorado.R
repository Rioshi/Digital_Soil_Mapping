library(caret)
library(doParallel)
library(googlesheets)
library(raster)
library(rasterVis)
library(spdep)
library(dplyr)
library(maptools)
library(agricolae)
library(VGAM)
library(nnet)

rm(list = ls())
##################################
### TRABAJO EN PARALELO ##########
##################################
detectCores()
doParallel::registerDoParallel(7)



##################################
### LECTURA DE DATOS ##########
##################################
setwd("F:/TESIS/2018")
covariables <- readRDS(file = "F:/TESIS/2018/RDATA/covariables.rds")
cali <- read.csv("calic.csv",header = TRUE)
cali <- cali[,1:4]

##################################
### Extraccion Covariables #######
##################################
cali.sp <- cali
coordinates(cali.sp) = ~X+Y
proj4string(cali.sp) <- CRS("+proj=utm +zone=18 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
cali <- cbind(cali, extract(covariables, cali.sp))
cali <- cali[,4:25]




##################################
### Modelamiento LOGISTIC ########
##################################
hiperparametros <- data.frame(parameter = "none")
control_train <- trainControl(method = "repeatedcv", number = particiones,
                              repeats = repeticiones,
                              classProbs = TRUE, 
                              summaryFunction = twoClassSummary, 
                              returnResamp = "final", verboseIter = FALSE, 
                              allowParallel = TRUE) #para que trabaje en paralelo

modelo_logistic <- train(ST ~ ., data = cali,
                         method = "glm",
                         tuneGrid = hiperparametros,
                         metric = "ROC",
                         trControl = control_train,
                         family = "binomial")
##################################
### Modelamiento nnet ############
##################################
md1 <- multinom(ST ~ ., data = cali)
summary(md1)
exp(coef(md1))
