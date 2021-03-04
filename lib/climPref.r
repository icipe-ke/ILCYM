loadRasters <- function(strPathMin, strPathMax){
  archivosMin=list.files(strPathMin,pattern="flt")
  archivosMax=list.files(strPathMax,pattern="flt")

  archivosMin=paste(strPathMin,"/",archivosMin,sep="")
  archivosMax=paste(strPathMax,"/",archivosMax,sep="")

  rMin<- list(1)
  rMax <- list(1)
  for(i in 1:length(archivosMax)){
    #paste("rMin",i,sep="") <- raster(archivosMin[i])
    #paste("rMax",i,sep="") <- raster(archivosMax[i])
    rMin[[i]] <-raster(archivosMin[i])
    rMax[[i]] <-raster(archivosMax[i])
  }

  #return(list(rMin1=rMin1, rMin2=rMin2,rMin3=rMin3,rMin4=rMin4,rMin5=rMin5,rMin6=rMin6,rMin7=rMin7,rMin8=rMin8,rMin9=rMin9,rMin10=rMin10,rMin11=rMin11,rMin12=rMin12,
  #            rMax1=rMax1, rMax2=rMax2,rMax3=rMax3,rMax4=rMax4,rMax5=rMax5,rMax6=rMax6,rMax7=rMax7,rMax8=rMax8,rMax9=rMax9,rMax10=rMax10,rMax11=rMax11,rMax12=rMax12,))
  return(list(rMin=rMin,rMax=rMax))
}
