setwd("H:/Worldclim2")

## Leer lista de paquetes ##
library(raster)
library(car)
library(rgdal)
library(RStoolbox)

#######################
####Lectura de archivos
star.dir <- getwd()


#Temperatura maxima
setwd("H:/Worldclim2/Tmax")
files <- list.files(pattern=".tif$")
t_max <- list()
for(i in 1:length(files)) {
  t_max[[i]] <- stack(files[i])
  }
t_max <- do.call(stack,t_max)


#Temperatura minima
setwd("H:/Worldclim2/Tmin")
files <- list.files(pattern=".tif$")
t_min <- list()
for(i in 1:length(files)) {
  t_min[[i]] <- stack(files[i])}
t_min <- do.call(stack,t_min)

#Temperatura promedio
T_mean <- list()
for(i in 1:12) {
  T_mean[[i]] <- stack((t_max[[i]]+t_min[[i]])/2)}
T_mean <- do.call(stack,T_mean)


#Precipitacion
setwd("H:/Worldclim2/Prec")
files <- list.files(pattern=".tif$")
files <- files[-13]
prec <- list()
for(i in 1:length(files)) {
  prec[[i]] <- stack(files[i])}
prec <- do.call(stack,prec)

#Evapotranspiracion 
setwd("H:/Worldclim2/ET_hagreveas")
files <- list.files(pattern=".tif$")
HS <- list()
for(i in 1:length(files)) {
  HS[[i]] <- stack(files[i])}
HS <- do.call(stack,HS)

#Radiacion
setwd("H:/Worldclim2/Rad")
files <- list.files(pattern=".tif$")
Rd <- list()
for(i in 1:length(files)) {
  Rd[[i]] <- stack(files[i])}
Rd <- do.call(stack,Rd)

#Vapor de agua
setwd("H:/Worldclim2/Humedad")
files <- list.files(pattern=".tif$")
Hd <- list()
for(i in 1:length(files)) {
  Hd[[i]] <- stack(files[i])}
Hd <- do.call(stack,Hd)

#Velocidad de tiento
setwd("H:/Worldclim2/Viento")
files <- list.files(pattern=".tif$")
Wind <- list()
for(i in 1:length(files)) {
  Wind[[i]] <- stack(files[i])}
Wind <- do.call(stack,Wind)




####################CALCULOS NO VIABLES#################################
#Latitud
latit <- raster("H:/Worldclim2/latitud.tif")

#Evapotranspiracion con Rsaga
library(RSAGA)
env <- rsaga.env("C:/Program Files (x86)/SAGA-GIS")
rsaga.get.libraries(path = env$modules)
rsaga.get.modules(libs = "climate_tools", env = env)
rsaga.get.usage(lib = "climate_tools", module = "ETpot (after Hargreaves, Grid)", env = env)

rsaga.geoprocessor(lib = "climate_tools", module = "ETpot (after Hargreaves, Grid)",
                   param = list(T = T_mean[[1]],
                                T_MIN = t_min[[1]],
                                T_MAX = t_max[[1]],
                                PET = "pet",
                                LAT = latit,
                                TIME = 1,
                                MONTH = 0
                                ),
                   env = env)
####################FIN DE CALCULOS NO VIABLES#################################

#Corte de capas
aoi <- shapefile("H:/TESIS/2018/Limites/extension_corte/Grid System Extent.shp")
aoi <- spTransform(aoi, crs(t_max[[1]]))

t_max <- crop(t_max,aoi)
t_min <- crop(t_min,aoi)
T_mean <- crop(T_mean,aoi)
prec <- crop(prec,aoi)
HS <- crop(HS,aoi)
Rd <- crop(Rd,aoi)
Hd <- crop(Hd,aoi)
Wind <- crop(Wind,aoi)

#Balance hidrico
BH <- list()
for(i in 1:12) {
  BH[[i]] <- stack(prec[[i]]-HS[[i]])}
BH <- do.call(stack,BH)
BH[BH < 0] <- 0

#Asegurar nombre de layers
names(t_max)
names(t_min)
names(T_mean)
nam_tmean <- paste("tmean",c(1,10:12,2:9),sep = "")
names(T_mean) <- nam_tmean
names(prec)
names(HS)
nam_BH <- paste("BH",c(1,10:12,2:9),sep="")
names(BH) <- nam_BH


###Guardar los objetos
save(HS,prec,T_mean,t_min,t_max,Rd,Hd,Wind, file = "H:/TESIS/2018/RDATA/clima.RData")
