###########3 MASK ###############
maskILCYM <- function(strInput, strMask, strOutput){
  r1 <- raster(strInput)
  sh1 <- readShapeSpatial(strMask)

  e <- extent(sh1@bbox[[1]], sh1@bbox[[3]], sh1@bbox[[2]], sh1@bbox[[4]])
  r1 <- crop(r1,e)
  
  r2<- rasterize(sh1,r1)
  
  rm1 <- mask(r1,r2, reverse = TRUE)
  
  writeRaster(rm1,filename=strOutput, overwrite=TRUE)
  xc <- rm1@extent@xmin
  yc <- rm1@extent@ymin
  nr <- rm1@nrows
  nc <- rm1@ncols
  
  cs <- (r1@extent@xmax - r1@extent@xmin) / r1@ncols
  
  rm(r1)
  rm(r2)
  return(list(xc=xc,yc=yc,nr=nr,nc=nc, cs=cs))
}
#####################################

############### REACLASS ###########
reclassILCYM <- function(strInput, strOutput, strReclass){
  r1 <- raster(strInput)
  m <- strReclass
  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  rc <- reclass(r1,rclmat, include.lowest=FALSE, right= NA)
  writeRaster(rc,filename=strOutput, overwrite=TRUE)
  
  xc <- rc@extent@xmin
  yc <- rc@extent@ymin
  nr <- rc@nrows
  nc <- rc@ncols

  cs <- (r1@extent@xmax - r1@extent@xmin) / r1@ncols

  rm(r1)
  return(list(xc=xc,yc=yc,nr=nr,nc=nc, cs=cs))

}
####################################

############ OVERLAY ##############
overlayILCYM <- function(strInput, strOverlay, strOutput, strFunc){
  r1 <- raster(strInput)
  r2 <- raster(strOverlay)
  
  if(strFunc == "suma"){
    ro <- overlay(r1, r2, fun= sum)
  }
  if(strFunc == "resta"){
    ro <- r1 - (r2)
  }
  if(strFunc == "por"){
    ro <- r1 * r2
  }
  if(strFunc == "entre"){
    ro <- r1 / r2
  }
  if(strFunc == "min"){
    ro <- overlay(r1, r2, fun= min)
  }
  if(strFunc == "max"){
    ro <- overlay(r1, r2, fun= max)
  }
  

  
  writeRaster(ro,filename=strOutput, overwrite=TRUE)

  xc <- ro@extent@xmin
  yc <- ro@extent@ymin
  nr <- ro@nrows
  nc <- ro@ncols

  cs <- (r1@extent@xmax - r1@extent@xmin) / r1@ncols

  rm(r1)
  return(list(xc=xc,yc=yc,nr=nr,nc=nc, cs=cs))
}
####################################

############ AGGREGATE ##############
aggregateILCYM <- function(strInput, strOutput, strKind, strFact, strExpand){
  r1 <- raster(strInput)
  ra <- aggregate(r1,fact= strFact , fun= strKind, expand= strExpand,na.rm=TRUE);

  writeRaster(ra,filename=strOutput, overwrite=TRUE)

  xc <- ra@extent@xmin
  yc <- ra@extent@ymin
  nr <- ra@nrows
  nc <- ra@ncols

  cs <- (r1@extent@xmax - r1@extent@xmin) / r1@ncols

  rm(r1)
  return(list(xc=xc,yc=yc,nr=nr,nc=nc, cs=cs))
}
####################################

############ DISAGGREGATE ##############
disaggregateILCYM <- function(strInput, strOutput, strFactor){
  r1 <- raster(strInput)
  rd <- disaggregate(r1, fact= strFactor)

  writeRaster(rd,filename=strOutput, overwrite=TRUE)

  xc <- rd@extent@xmin
  yc <- rd@extent@ymin
  nr <- rd@nrows
  nc <- rd@ncols

  cs <- (r1@extent@xmax - r1@extent@xmin) / r1@ncols

  rm(r1)
  return(list(xc=xc,yc=yc,nr=nr,nc=nc, cs=cs))
}

####################################

############ CUT ##############
cutILCYM <- function(strInput, strOutput, strMinX, strMaxX, strMinY, strMaxY){
  r1 <- raster(strInput)
  e <- extent(strMinX, strMaxX, strMinY, strMaxY)
  rd <- crop(r1, e)
  
  writeRaster(rd,filename=strOutput, overwrite=TRUE)

  xc <- rd@extent@xmin
  yc <- rd@extent@ymin
  nr <- rd@nrows
  nc <- rd@ncols

  cs <- (r1@extent@xmax - r1@extent@xmin) / r1@ncols

  rm(r1)
  return(list(xc=xc,yc=yc,nr=nr,nc=nc, cs=cs))
}
####################################

############ CUT ##############
mergeILCYM <- function(strInput, strInput2, strOutput){
  r1 <- raster(strInput)
  r2 <- raster(strInput2)
  
  rd <- merge(r1,r2)
  
  writeRaster(rd,filename=strOutput, overwrite=TRUE)

  xc <- rd@extent@xmin
  yc <- rd@extent@ymin
  nr <- rd@nrows
  nc <- rd@ncols

  cs <- (r1@extent@xmax - r1@extent@xmin) / r1@ncols

  rm(r1)
  return(list(xc=xc,yc=yc,nr=nr,nc=nc, cs=cs))
}
####################################

############### CALC #################
calcILCYM <- function(strOutput, dataFinal, xmin=xmin,xmax=xmax, ymin=ymin,ymax=ymax){
  r1 <- raster(dataFinal, xmn=xmin,xmx=xmax, ymn=ymin,ymx=ymax, crs=NA)

  writeRaster(r1,filename=strOutput, overwrite=TRUE)

  xc <- r1@extent@xmin
  yc <- r1@extent@ymin
  nr <- r1@nrows
  nc <- r1@ncols

  cs <- (r1@extent@xmax - r1@extent@xmin) / r1@ncols

  rm(r1)
  return(list(xc=xc,yc=yc,nr=nr,nc=nc, cs=cs))
}
####################################
terrainFunction <- function(strInput, strOutput, stropt){
library(raster)
r1 <- raster(strInput)
t1 <- terrain(r1, opt=stropt,unit="degrees")

writeRaster(t1,filename=strOutput, overwrite=TRUE)
}
###################################
overlayStack <- function(strInput, strOutput, fun1){
library(raster)
s1 <- stackOpen(strInput)
ro <- overlay(s1, fun=fun1)
writeRaster(ro,filename=strOutput, overwrite=TRUE)
}

###################################
rasterToPoly <- function(r1, strOutput){
pol1 <- rasterToPolygons(r1, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=TRUE)
writePolyShape(pol1, fn=strOutput, factor2char = TRUE, max_nchar=254)
}
###################################

