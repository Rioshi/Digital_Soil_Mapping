#Lectura de covariables litologia
library(raster)
startdir <- getwd()
setwd("H:/TESIS/2018/ASTER L1T/Procesado-07-DIC-2015")

limite <- shapefile("H:/TESIS/2018/Limites/extension_corte/Grid System Extent.shp")
limite <- spTransform(limite, crs("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84"))


aster <- list.files(pattern="tif$")
stack1 <- list()
for(i in 1:length(aster)) {
  stack1[[i]] <- crop(raster(aster[i]),limite)
  }
aster <- do.call(stack, stack1)
rm(i,limite,stack1,startdir)
save(aster, file = "H:/TESIS/2018/RDATA/aster.RData")
