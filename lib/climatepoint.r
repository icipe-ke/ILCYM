getTemps <- function(strPathMin, strPathMax, rowIndex, colIndex, numCols, numRows){
  library(sp)
  library(raster)
  index <-  rowIndex * numCols + colIndex +1
  temps <- matrix(NA, 12, 2)


  for(i in 1:length(rMin)){
    valMin <- rMin[[i]][index]
    valMax <- rMax[[i]][index]

    temps[i,] <-c(valMin, valMax)

  }
  
  return(temps=temps)
}
##########################################################################
extractFunction <- function(inputShpPath, strOutput, bolClima, strRaster, strRasterName){

	library(raster)
	library(maptools)


	inputShapeFile <- readShapePoints(inputShpPath, proj4string = CRS(as.character(NA)), verbose = FALSE, repair=FALSE)

	shpSpatialPoint <- SpatialPoints(inputShapeFile)

	#outputFile <- data.frame(shpSpatialPoint@coords[,1], shpSpatialPoint@coords[,2])
		outputFile <- inputShapeFile@data
	
	nombresTmin <- ""
	nombresTmax <- ""
  if(bolClima){
	  for(i in 1:length(rMin)){
		  tempsValues <- extract(rMin[[i]], shpSpatialPoint)/10
		  outputFile <- cbind(outputFile, tempsValues)	
		  nombresTmin <- c(nombresTmin, paste("tmin_",i))
	  }

	  for(i in 1:length(rMax)){
		  tempsValues <- extract(rMax[[i]], shpSpatialPoint)/10
		  outputFile <- cbind(outputFile, tempsValues)	
		  nombresTmax <- c(nombresTmax, paste("tmax_",i))
	  }
	  	nombreColumnas <- c(colnames(inputShapeFile@data),nombresTmin[2:length(nombresTmin)], nombresTmax[2:length(nombresTmax)])
	}else{
	  r1 <- raster(strRaster)
	  tempsValues <- extract(r1, shpSpatialPoint)
	  outputFile <- cbind(outputFile, tempsValues)	
	  nombreColumnas <- c(colnames(inputShapeFile@data),strRasterName)
	  
	}

	#nombreColumnas <- c(nameLong, nameLat,nombresTmin[2:length(nombresTmin)], nombresTmax[2:length(nombresTmax)])

	colnames(outputFile) <- nombreColumnas

	write.table(outputFile, strOutput, row.names=FALSE, sep="\t")

}
############################################
extractFunctionE <- function(nameLong, nameLat, inputShpPath, strOutput){

	library(raster)
	library(maptools)


	inputShapeFile <- readShapePoints(inputShpPath, proj4string = CRS(as.character(NA)), verbose = FALSE, repair=FALSE)

	shpSpatialPoint <- SpatialPoints(inputShapeFile)

	outputFile <- data.frame(shpSpatialPoint@coords[,1], shpSpatialPoint@coords[,2])
	
	nombresTmin <- ""
	nombresTmax <- ""

	for(i in 1:length(rMin)){
		tempsValues <- extract(rMin[[i]], shpSpatialPoint)/10
		outputFile <- cbind(outputFile, tempsValues)	
		nombresTmin <- c(nombresTmin, paste("tmin_",i))
	}

	for(i in 1:length(rMax)){
		tempsValues <- extract(rMax[[i]], shpSpatialPoint)/10
		outputFile <- cbind(outputFile, tempsValues)	
		nombresTmax <- c(nombresTmax, paste("tmax_",i))
	}

	nombreColumnas <- c(nameLong, nameLat,nombresTmin[2:length(nombresTmin)], nombresTmax[2:length(nombresTmax)])
	colnames(outputFile) <- nombreColumnas

	write.table(outputFile, strOutput, row.names=FALSE, sep="\t")

}
###############################################
climatePointFunction <- function(dbLong, dbLat){
    xy <- data.frame(dbLong, dbLat)
    shpSpatialPoint <- SpatialPoints(xy)
    outputFile <- data.frame(0,0)
    
  	for(i in 1:length(rMin)){
      tempsValuesMin <- extract(rMin[[i]], shpSpatialPoint)/10
      tempsValuesMax <- extract(rMax[[i]], shpSpatialPoint)/10
      outputFile <- rbind(outputFile, c(tempsValuesMin, tempsValuesMax))
    }  
    
    outputFile <- outputFile[2:nrow(outputFile),]
    
    return(data.frame(outputFile))
}

