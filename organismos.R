library(water)
aoi <- shapefile("H:/TESIS/2018/Limites/extension_corte/Grid System Extent.shp")
l8p <- raster("H:/TESIS/2018/SRL8/LC08_L1TP_007068_20140603_20170422_01_T1_sr_band1.tif")
aoi <- spTransform(aoi, crs(l8p))
rm(l8p)

l8 <- loadImageSR(path="H:/TESIS/2018/SRL8",aoi)
values(l8)[values(l8) <= 0] = NA


#Indice de vegetacion
ndvi <- (l8$NIR-l8$R)/(l8$NIR+l8$R)
values(ndvi)[values(ndvi) <= 0] = NA


#Costra biologica de Kraemeri
ci <- 1 - ((l8$R-l8$B)/(l8$R+l8$B))

#Biological soil crust index de Chen
bsci <- (1-(2*abs(l8$R-l8$G)))/((l8$R+l8$G+l8$NIR)/3)
values(bsci)[values(bsci) > 6] = NA

rm(aoi,l8)
save(ndvi,ci, file = "H:/TESIS/2018/RDATA/organismos.RData")
