######### DEVELOPMENT RATE ###############
comp.species<-function(dirs,especies,ni,corrx,corry,labx=labx,laby=laby,titulo,scaleX,scaleY,tam,leyen=FALSE)
{
 col1=c("#6BAED6","#5FAE27","darkred","orange","green1","#5844ff") ## se puede aumentar mas colores pero consideramos solo hasta 6 especies
 N1=length(dirs)

 par(cex=tam)
 corrx2=seq(corrx[1],corrx[2],scaleX)
 corry2=seq(corry[1],corry[2],scaleY)

 plot(c(corrx[1],corrx[2]),c(corry[1],corry[2]),frame=F,type="n",pch=19,xlim=corrx,ylim=corry,xlab=list(labx,font=2,cex=1.3),ylab=list(laby,font=2,cex=1.3),axes=F,xaxt = "n", main=titulo) ## 4
 for(i in 1:length(corry2)){lines(c(0,corrx[2]),c(corry2[i],corry2[i]),lty=1,lwd=1.5,col="gray60",type = "l")}
 axis(1, corrx2)  ## cambio
 axis(2, corry2,las=2) ## cambio


 for(j in 1:N1)
 {
	load(dirs[j])
	ini<-params$hfeno$pdv_dv[[ni]]
	expre<-params$hfeno$fdv_dv[[ni]]
	for (i in names(ini)){ temp <- ini[[i]];storage.mode(temp) <- "double";assign(i, temp)	}
	if(params$modelim[ni]<=14){x<-seq(0,corrx[2],length=1000)+273.15}else{x<-seq(0,corrx[2],length=1000)}
	x0<-seq(0,corrx[2],length=1000)
	y1<-eval(expre[[3]]);y1[y1<0]=0

	points(x0,y1,type="l",col=col1[j],lwd=2)
 }
 if(leyen){legend("topright",especies,col =col1,lty = 1,lwd=2,cex=0.8,bg = 'white')}
}
##########################################

######### SENESCENCIA HEMBRAS ###############
comp.speciesH<-function(dirs,especies,corrx,corry,labx=labx,laby=laby,titulo,scaleX,scaleY,tam,leyen=FALSE)
{
 col1=c("#6BAED6","#5FAE27","darkred","orange","green1","#5844ff") ## se puede aumentar mas colores pero consideramos solo hasta 6 especies
 N1=length(dirs)

 par(cex=tam)
 corrx2=seq(corrx[1],corrx[2],scaleX)
 corry2=seq(corry[1],corry[2],scaleY)

 plot(c(corrx[1],corrx[2]),c(corry[1],corry[2]),frame=F,type="n",pch=19,xlim=corrx,ylim=corry,xlab=list(labx,font=2,cex=1.3),ylab=list(laby,font=2,cex=1.3),axes=F,xaxt = "n", main=titulo) ## 4
 for(i in 1:length(corry2)){lines(c(0,corrx[2]),c(corry2[i],corry2[i]),lty=1,lwd=1.5,col="gray60",type = "l")}
 axis(1, corrx2)  ## cambio
 axis(2, corry2,las=2) ## cambio


 for(j in 1:N1)
 {
	load(dirs[j])
	ini<-params$hfeno$pfh_h
	expre<-params$hfeno$fh_h
	for (i in names(ini)){ temp <- ini[[i]];storage.mode(temp) <- "double";assign(i, temp)	}
	if(params$modelm[1]<=14){x<-seq(0,corrx[2],length=1000)+273.15}else{x<-seq(0,corrx[2],length=1000)}
	x0<-seq(0,corrx[2],length=1000)
	y1<-eval(expre[[3]]);y1[y1<0]=0

	points(x0,y1,type="l",col=col1[j],lwd=2)
 }
 if(leyen){legend("topright",especies,col =col1,lty = 1,lwd=2,cex=0.8,bg = 'white')}
}
##########################################

######### SENESCENCIA MACHOS ###############
comp.speciesM<-function(dirs,especies,corrx,corry,labx=labx,laby=laby,titulo,scaleX,scaleY,tam,leyen=FALSE)
{
 col1=c("#6BAED6","#5FAE27","darkred","orange","green1","#5844ff") ## se puede aumentar mas colores pero consideramos solo hasta 6 especies
 N1=length(dirs)

 par(cex=tam)
 corrx2=seq(corrx[1],corrx[2],scaleX)
 corry2=seq(corry[1],corry[2],scaleY)

 plot(c(corrx[1],corrx[2]),c(corry[1],corry[2]),frame=F,type="n",pch=19,xlim=corrx,ylim=corry,xlab=list(labx,font=2,cex=1.3),ylab=list(laby,font=2,cex=1.3),axes=F,xaxt = "n", main=titulo) ## 4
 for(i in 1:length(corry2)){lines(c(0,corrx[2]),c(corry2[i],corry2[i]),lty=1,lwd=1.5,col="gray60",type = "l")}
 axis(1, corrx2)  ## cambio
 axis(2, corry2,las=2) ## cambio


 for(j in 1:N1)
 {
	load(dirs[j])
	ini<-params$hfeno$pfm_m
	expre<-params$hfeno$fm_m
	for (i in names(ini)){ temp <- ini[[i]];storage.mode(temp) <- "double";assign(i, temp)	}
	if(params$modelm[2]<=14){x<-seq(0,corrx[2],length=1000)+273.15}else{x<-seq(0,corrx[2],length=1000)}
	x0<-seq(0,corrx[2],length=1000)
	y1<-eval(expre[[3]]);y1[y1<0]=0

	points(x0,y1,type="l",col=col1[j],lwd=2)
 }
 if(leyen){legend("topright",especies,col =col1,lty = 1,lwd=2,cex=0.8,bg = 'white')}
}
##################################

######### MORTALIDAD ###############
mort.speciesM<-function(dirs,especies,ni,corrx,corry,labx=labx,laby=laby,titulo,scaleX,scaleY,tam,leyen=FALSE)
{
 col1=c("#6BAED6","#5FAE27","darkred","orange","green1","#5844ff") ## se puede aumentar mas colores pero consideramos solo hasta 6 especies
 N1=length(dirs)

 par(cex=tam)
 corrx2=seq(corrx[1],corrx[2],scaleX)
 corry2=seq(corry[1],corry[2],scaleY)

 plot(c(corrx[1],corrx[2]),c(corry[1],corry[2]),frame=F,type="n",pch=19,xlim=corrx,ylim=corry,xlab=list(labx,font=2,cex=1.3),ylab=list(laby,font=2,cex=1.3),axes=F,xaxt = "n", main=titulo) ## 4
 for(i in 1:length(corry2)){lines(c(0,corrx[2]),c(corry2[i],corry2[i]),lty=1,lwd=1.5,col="gray60",type = "l")}
 axis(1, corrx2)  ## cambio
 axis(2, corry2,corry2*100,las=2) ## cambio


 for(j in 1:N1)
 {
	load(dirs[j])
	ini<-params$hfeno$pmortal[[ni]]
	expre<-params$hfeno$mortal[[ni]]
	for (i in names(ini)){ temp <- ini[[i]];storage.mode(temp) <- "double";assign(i, temp)	}
	x<-seq(0,corrx[2],length=1000)
	x0<-seq(0,corrx[2],length=1000)
	y1<-eval(expre[[3]]);y1[y1<0]=0

	points(x0,y1,type="l",col=col1[j],lwd=2)
 }
 if(leyen){legend("topright",especies,col =col1,lty = 1,lwd=2,cex=0.8,bg = 'white')}
}
##################################

######### FECUNDIDAD ##############
fecun.speciesM<-function(dirs,especies,corrx,corry,labx=labx,laby=laby,titulo,scaleX,scaleY,tam,leyen=FALSE)
{
 col1=c("#6BAED6","#5FAE27","darkred","orange","green1","#5844ff") ## se puede aumentar mas colores pero consideramos solo hasta 6 especies
 N1=length(dirs)

 par(cex=tam)
 corrx2=seq(corrx[1],corrx[2],scaleX)
 corry2=seq(corry[1],corry[2],scaleY)

 plot(c(corrx[1],corrx[2]),c(corry[1],corry[2]),frame=F,type="n",pch=19,xlim=corrx,ylim=corry,xlab=list(labx,font=2,cex=1.3),ylab=list(laby,font=2,cex=1.3),axes=F,xaxt = "n", main=titulo) ## 4
 for(i in 1:length(corry2)){lines(c(0,corrx[2]),c(corry2[i],corry2[i]),lty=1,lwd=1.5,col="gray60",type = "l")}
 axis(1, corrx2)  ## cambio
 axis(2, corry2,las=2) ## cambio


 for(j in 1:N1)
 {
	load(dirs[j])
	ini<-params$hfeno$ptazaeh_h
	expre<-params$hfeno$ftazaeh_h
	for (i in names(ini)){ temp <- ini[[i]];storage.mode(temp) <- "double";assign(i, temp)	}
	x<-seq(0,corrx[2],length=1000)
	x0<-seq(0,corrx[2],length=1000)
	y1<-eval(expre[[3]]);y1[y1<0]=0

	points(x0,y1,type="l",col=col1[j],lwd=2)
 }
 if(leyen){legend("topright",especies,col =col1,lty = 1,lwd=2,cex=0.8,bg = 'white')}
}

##################################


crearReporte <- function(dir1, dirOut, nameFile, nameProject, spProject, autorName, strDate){
load(dir1) ## solo necesito una fenologia

 estadios=params$estadios
 n1=length(params$modelim);c1<-c2<-c3<-c(1);valor1<-valor2<-valor3<-list(1);etiq1<-etiq2<-etiq3<-list(1);#eti2=c(1);val2=c(1)
 StEr1<-StEr2<-StEr3<-list(1); StE2=c(1)										#Agregado por DM

 for(ni in 1:n1){
 c1[ni]=length(params$hfeno$pdv_dv[[ni]]);c2[ni]=length(params$hfeno$pmortal[[ni]]); c3[ni]=length(params$hfeno$slope_dv[[ni]])}
 nn1=max(c1);nn2=max(c2); nn3=max(c3)											

 datF1=data.frame(v1=rep(NA,nn1))
 datF3=data.frame(v1=rep(NA,nn2))
 datF2=data.frame(v1=rep(NA,nn3))

 DRate=params$hfeno$pdv_dv; Mort=params$hfeno$pmortal; SE_DR=params$hfeno$StEr_DRate; 	SE_M=params$hfeno$StEr_Mort		#Agregado por DM 
 et2=params$hfeno$distri_dv; DTime=params$hfeno$slope_dv; SE_DT=params$hfeno$StEr_DTime
 
 for(ni in 1:n1)
 {

  for(i in 1:length(DRate[[ni]])) { if(abs(DRate[[ni]][i]) < 0.00001) {DRate[[ni]][i]<-0.00001}; if(SE_DR[[ni]][i]!="NaN") {	if(SE_DR[[ni]][i]<0.00001) {SE_DR[[ni]][i]<-0.00001}}	}	#Agregado por DM

  if(nn1 > c1[ni])
  {
	valor1[[ni]]=c(round(DRate[[ni]],5),rep("",nn1-c1[ni]))
	etiq1[[ni]]=labels(valor1[[ni]])
	StEr1[[ni]]=as.character(c(paste(" ± ",round(SE_DR[[ni]],5)),rep("",nn1-c1[[ni]]))) 			#Agregado por DM
  }
  else
  {
	valor1[[ni]]=as.character(round(DRate[[ni]],5))
	etiq1[[ni]]=labels(DRate[[ni]])
	StEr1[[ni]]=as.character(paste(" ± ",round(SE_DR[[ni]],5))) 			#Agregado por DM
  }
  
  for(i in 1:length(DRate[[ni]])) if( SE_DR[[ni]][i]!="NaN") {	if(DRate[[ni]][i] == SE_DR[[ni]][i]) {StEr1[[ni]][i]<-as.character(paste(" ± ",0.000001))}	}	else	  StEr1[[ni]][i]<-" "	 	#Agregado por DM

  datF1=data.frame(datF1,etiqueta=etiq1[[ni]],valor=paste(valor1[[ni]],StEr1[[ni]]),"|")			#Cambie

#########################################

  for(i in 1:length(DTime[[ni]])) { if(abs(DTime[[ni]][i]) < 0.00001) {DTime[[ni]][i]<-0.00001}; if(SE_DT[[ni]][i]<0.00001) {SE_DT[[ni]][i]<-0.00001} }			#Agregado por DM

  if(nn3 > c3[ni])
  {
	valor3[[ni]]=c(round(DTime[[ni]],5),rep("",nn3-c3[ni]))
	etiq3[[ni]]=labels(valor3[[ni]])
	StEr3[[ni]]=as.character(c(paste(" ± ",round(SE_DT[[ni]],5)),rep("",nn3-c3[[ni]]))) 			#Agregado por DM
  }
  else
  {
	valor3[[ni]]=as.character(round(DTime[[ni]],5))
	etiq3[[ni]]=labels(DTime[[ni]])
	StEr3[[ni]]=as.character(paste(" ± ",round(SE_DT[[ni]],5))) 			#Agregado por DM
  }

  for(i in 1:length(DTime[[ni]])) { if(DTime[[ni]][i] == SE_DT[[ni]][i]) {StEr3[[ni]][i]<-as.character(paste(" ± ",0.000001))};   	 	#Agregado por DM
	if(etiq3[[ni]][i] == 1 || etiq3[[ni]][i] == 2 || etiq3[[ni]][i] == 3 || etiq3[[ni]][i] == 4) etiq3[[ni]][i]="slope" } #ALTERNATIVA A SLOPE ################################################	
  datF2=data.frame(datF2,etiqueta=etiq3[[ni]],valor=paste(valor3[[ni]],StEr3[[ni]]),"|")			#Cambie 


###
#  if(abs(DTime[[ni]])<0.00001) DTime[[ni]]<-0.00001; if(SE_DT[[ni]]<0.00001) SE_DT[[ni]]<-0.00001; if(DTime[[ni]]==SE_DT[[ni]]) SE_DT[[ni]]<-0.000001
#  eti2[ni]=et2[[ni]]; val2[ni]=round(DTime[[ni]],5)
#  StE2[ni]=round(SE_DT[[ni]],5)
########################################Error in et2[[ni]] : subscript out of bounds


  for(i in 1:length(Mort[[ni]])) { if(abs(Mort[[ni]][i]) < 0.00001) {Mort[[ni]][i]<-0.00001}; if(SE_M[[ni]][i]<0.00001) {SE_M[[ni]][i]<-0.00001} }			#Agregado por DM

  if(nn2 > c2[ni])
  {
	valor2[[ni]]=c(round(Mort[[ni]],5),rep("",nn2-c2[ni]))
	etiq2[[ni]]=labels(valor2[[ni]])
	StEr2[[ni]]=as.character(c(paste(" ± ",round(SE_M[[ni]],5)),rep("",nn2-c2[[ni]]))) 			#Agregado por DM
  }
  else
  {
	valor2[[ni]]=as.character(round(Mort[[ni]],5))
	etiq2[[ni]]=labels(Mort[[ni]])
	StEr2[[ni]]=as.character(paste(" ± ",round(SE_M[[ni]],5))) 			#Agregado por DM
  }

  for(i in 1:length(Mort[[ni]])) if(Mort[[ni]][i] == SE_M[[ni]][i]) {StEr2[[ni]][i]<-as.character(paste(" ± ",0.000001))}  	 	#Agregado por DM

  datF3=data.frame(datF3,etiqueta=etiq2[[ni]],valor=paste(valor2[[ni]],StEr2[[ni]]),"|")			#Cambie 

 }

datF1=datF1[,-1];datF2=datF2[,-1];datF3=datF3[,-1];	datF1=datF1[,-ncol(datF1)]; datF2=datF2[,-ncol(datF2)]; datF3=datF3[,-ncol(datF3)]

################################### INICIO ETIQUETAS ###################################

 labcol=c(1); for(i in 1:n1) labcol=c(labcol,paste(estadios[i],"label"),paste(estadios[i],"value"),"|"); labcol=labcol[-1]; labcol=labcol[-length(labcol)]                      #Cambie
# labcol=paste(rep(estadios[1:n1],rep(2,n1)),rep(c("label","value"),n1))
 colnames(datF1)=labcol;colnames(datF2)=labcol;colnames(datF3)=labcol

################################### FIN ETIQUETAS ######################################

################################### INICIO DEVELOPMENT_TIME MADUROS ##################################

 DTimeH=params$hfeno$slope_snh; SE_DTH=params$hfeno$StEr_DTimeH 
 DTimeM=params$hfeno$slope_snm; SE_DTM=params$hfeno$StEr_DTimeM

 for(i in 1:length(DTimeH)) { if(abs(DTimeH[i])<0.00001) DTimeH[i]<-0.00001; if(SE_DTH[i]<0.00001) SE_DTH[i]<-0.00001; if(DTimeH[i]==SE_DTH[i]) SE_DTH[i]<-0.000001} 
 eti6_H=labels(DTimeH); val6_H=round(DTimeH,5); StE6_H=round(SE_DTH,5)
 datF6_H=data.frame(Label=eti6_H,Value=paste(val6_H," ± ",StE6_H)); rownames(datF6_H)=1:length(eti6_H)			#Agregado por DM

 for(i in 1:length(DTimeM)) { if(abs(DTimeM[i])<0.00001) DTimeM[i]<-0.00001; if(SE_DTM[i]<0.00001) SE_DTM[i]<-0.00001; if(DTimeM[i]==SE_DTM[i]) SE_DTM[i]<-0.000001} 
 eti6_M=labels(DTimeM); val6_M=round(DTimeM,5); StE6_M=round(SE_DTM,5)
 datF6_M=data.frame(Label=eti6_M,Value=paste(val6_M," ± ",StE6_M)); rownames(datF6_M)=1:length(eti6_M)			#Agregado por DM

 maxDT=max(nrow(datF6_H),nrow(datF6_M))
 if(maxDT > nrow(datF6_H)) { blank=data.frame(matrix((rep(c(" "," "),maxS-nrow(datF6_H))),maxS-nrow(datF6_H),2)); colnames(blank)=c("Label","Value"); datF6_H=rbind(datF6_H,blank) }
 sepDT=cbind(rep("|",maxDT))
 datF6=cbind(datF6_H,sepDT,datF6_M)
 
 name=c(paste(estadios[length(estadios)-1],"label"),paste(estadios[length(estadios)-1],"value"),"|",paste(estadios[length(estadios)],"label"),paste(estadios[length(estadios)],"value"))
 colnames(datF6)=name
 if(length(params$modelm)==1) datF6=datF6[,c(1,2)]

################################### FIN DEVELOPMENT_TIME MADUROS #####################################

################################### INICIO SENESCENCIA ##################################

 SenH=params$hfeno$pfh_h; SE_SenH=params$hfeno$StEr_DRateH
 SenM=params$hfeno$pfm_m; SE_SenM=params$hfeno$StEr_DRateM

 for(i in 1:length(SenH)) { if(abs(SenH[i])<0.00001) SenH[i]<-0.00001; if(SE_SenH[i]<0.00001) SE_SenH[i]<-0.00001; if(SenH[i]==SE_SenH[i]) SE_SenH[i]<-0.000001} 
 eti7_H=labels(SenH); val7_H=round(SenH,5); StE7_H=round(SE_SenH,5)
 datF7_H=data.frame(Label=eti7_H,Value=paste(val7_H," ± ",StE7_H)); rownames(datF7_H)=1:length(eti7_H)			#Agregado por DM

 for(i in 1:length(SenM)) { if(abs(SenM[i])<0.00001) SenM[i]<-0.00001; if(SE_SenM[i]<0.00001) SE_SenM[i]<-0.00001; if(SenM[i]==SE_SenM[i]) SE_SenM[i]<-0.000001} 
 eti7_M=labels(SenM); val7_M=round(SenM,5); StE7_M=round(SE_SenM,5)
 datF7_M=data.frame(Label=eti7_M,Value=paste(val7_M," ± ",StE7_M)); rownames(datF7_M)=1:length(eti7_M)			#Agregado por DM

 maxS=max(nrow(datF7_H),nrow(datF7_M))
 if(maxS > nrow(datF7_H)) { blank=data.frame(matrix((rep(c(" "," "),maxS-nrow(datF7_H))),maxS-nrow(datF7_H),2)); colnames(blank)=c("Label","Value"); datF7_H=rbind(datF7_H,blank) }
 if(maxS > nrow(datF7_M)) { blank=data.frame(matrix((rep(c(" "," "),maxS-nrow(datF7_M))),maxS-nrow(datF7_M),2)); colnames(blank)=c("Label","Value"); datF7_M=rbind(datF7_M,blank) }
 sepS=cbind(rep("|",maxS))
 datF7=cbind(datF7_H,sepS,datF7_M)
 
 name=c(paste(estadios[length(estadios)-1],"label"),paste(estadios[length(estadios)-1],"value"),"|",paste(estadios[length(estadios)],"label"),paste(estadios[length(estadios)],"value"))
 colnames(datF7)=name
 if(length(params$modelm)==1) datF7=datF7[,c(1,2)]

################################### FIN SENESCENCIA ##################################

# rownames(datF1)=1:nn1;rownames(datF3)=1:nn2
# datF2=data.frame(Stage=estadios[1:n1],Link=eti2,Slope=paste(val2," ± ",StE2))				#Cambie

 TOvi=params$hfeno$ptazaeh_h; SE_TO=params$hfeno$StEr_ToviH
 for(i in 1:length(TOvi)) { if(abs(TOvi[i])<0.00001) TOvi[i]<-0.00001; if(SE_TO[i]<0.00001) SE_TO[i]<-0.00001; if(TOvi[i]==SE_TO[i]) SE_TO[i]<-0.000001} 
 eti4=labels(params$hfeno$ptazaeh_h); val4=round(TOvi,5); StE4=round(SE_TO,5)					#Agregado por DM
 datF4=data.frame(Label=eti4,Value=paste(val4," ± ",StE4)); rownames(datF4)=1:length(eti4)			#Agregado por DM
#datF4=data.frame(t(round(params$hfeno$ptazaeh_h,3)))

 ROvi=params$hfeno$povih_h; SE_RO=params$hfeno$StEr_RoviH
 for(i in 1:length(ROvi)) { if(abs(ROvi[i])<0.00001) ROvi[i]<-0.00001; if(SE_RO[i]<0.00001) SE_RO[i]<-0.00001; if(ROvi[i]==SE_RO[i]) SE_RO[i]<-0.000001} 
 eti5=labels(params$hfeno$povih_h); val5=round(ROvi,5); StE5=round(SE_RO,5)						#Agregado por DM
 datF5=data.frame(Label=eti5,Value=paste(val5," ± ",StE5)); rownames(datF5)=1:length(eti5)			#Agregado por DM
# datF5=data.frame(t(round(params$hfeno$povih_h,3))

spname <- paste("<b>Species name: </b>","<EM>",spProject,"</EM>")
 pname <- paste("<b>Project name: </b>",nameProject)
 autor <- paste("<b>Author name: </b>",autorName)
 sDate <- paste("<b>Compilation date: </b>",strDate)
 
 
 target <- HTMLInitFile(file.path(dirOut),filename=nameFile)

 HTML.title("Parameter Values of the Functions Used", HR=1)
HTML(spname,file=target)
HTML(pname,file=target)
HTML(autor,file=target)
HTML(sDate,file=target)
HTML("<br></br>",file=target)


 HTML.title("Development Rate", HR=3)
tmp00=datF1
HTML(tmp00,file=target)
HTML.title("Senescence", HR=3)
tmp06=datF7
HTML(tmp06,file=target)
HTML.title("Distribution function of development", HR=3)
tmp01=datF2; tmp05=datF6
HTML(tmp01,file=target)
HTML(tmp05,file=target)
HTML.title("Mortality", HR=3)
tmp02=datF3
HTML(tmp02,file=target)
HTML.title("Total eggs per female", HR=3)
tmp03=datF4
HTML(tmp03,file=target)												#Cambie
#HTML(tmp03,file=target,row.names = FALSE)
HTML.title("Proportion of progeny production with age", HR=3)
tmp04=datF5
HTML(tmp04,file=target)												#Cambie
#HTML(tmp04,file=target,row.names = FALSE)
HTMLEndFile()
}
#################################################################
#################################################################
#################################################################
crearReporteAnt <- function(dir1, dirOut, nameFile, nameProject, spProject, autorName, strDate){
load(dir1) ## solo necesito una fenologia

 estadios=params$estadios
 n1=length(params$modelim);c1<-c2<-c3<-c(1);valor1<-valor2<-valor3<-list(1);etiq1<-etiq2<-etiq3<-list(1);#eti2=c(1);val2=c(1)
 StEr1<-StEr2<-StEr3<-list(1); StE2=c(1)										#Agregado por DM

 for(ni in 1:n1){
 c1[ni]=length(params$hfeno$pdv_dv[[ni]]);c2[ni]=length(params$hfeno$pmortal[[ni]]); c3[ni]=length(params$hfeno$slope_dv[[ni]])}
 nn1=max(c1);nn2=max(c2); nn3=max(c3)											

 datF1=data.frame(v1=rep(NA,nn1))
 datF3=data.frame(v1=rep(NA,nn2))
 datF2=data.frame(v1=rep(NA,nn3))

 DRate=params$hfeno$pdv_dv; Mort=params$hfeno$pmortal; SE_DR=params$hfeno$StEr_DRate; 	SE_M=params$hfeno$StEr_Mort		#Agregado por DM 
 et2=params$hfeno$distri_dv; DTime=params$hfeno$slope_dv; SE_DT=params$hfeno$StEr_DTime
 for(ni in 1:n1)
 {

  for(i in 1:length(DRate[[ni]])) { if(abs(DRate[[ni]][i]) < 0.00001) {DRate[[ni]][i]<-0.00001}; if(SE_DR[[ni]][i]!="NaN") {	if(SE_DR[[ni]][i]<0.00001) {SE_DR[[ni]][i]<-0.00001}}	}	#Agregado por DM

  if(nn1 > c1[ni])
  {
	valor1[[ni]]=c(round(DRate[[ni]],5),rep("",nn1-c1[ni]))
	etiq1[[ni]]=labels(valor1[[ni]])
	StEr1[[ni]]=as.character(c(paste(" ± ",round(SE_DR[[ni]],5)),rep("",nn1-c1[[ni]]))) 			#Agregado por DM
  }
  else
  {
	valor1[[ni]]=as.character(round(DRate[[ni]],5))
	etiq1[[ni]]=labels(DRate[[ni]])
	StEr1[[ni]]=as.character(paste(" ± ",round(SE_DR[[ni]],5))) 			#Agregado por DM
  }
  
  for(i in 1:length(DRate[[ni]])) if( SE_DR[[ni]][i]!="NaN") {	if(DRate[[ni]][i] == SE_DR[[ni]][i]) {StEr1[[ni]][i]<-as.character(paste(" ± ",0.000001))}	}	else	  StEr1[[ni]][i]<-" "	 	#Agregado por DM

  datF1=data.frame(datF1,etiqueta=etiq1[[ni]],valor=paste(valor1[[ni]],StEr1[[ni]]),"|")			#Cambie

#########################################

  for(i in 1:length(DTime[[ni]])) { if(abs(DTime[[ni]][i]) < 0.00001) {DTime[[ni]][i]<-0.00001}; if(SE_DT[[ni]][i]<0.00001) {SE_DT[[ni]][i]<-0.00001} }			#Agregado por DM

  if(nn3 > c3[ni])
  {
	valor3[[ni]]=c(round(DTime[[ni]],5),rep("",nn3-c3[ni]))
	etiq3[[ni]]=labels(valor3[[ni]])
	StEr3[[ni]]=as.character(c(paste(" ± ",round(SE_DT[[ni]],5)),rep("",nn3-c3[[ni]]))) 			#Agregado por DM
  }
  else
  {
	valor3[[ni]]=as.character(round(DTime[[ni]],5))
	etiq3[[ni]]=labels(DTime[[ni]])
	StEr3[[ni]]=as.character(paste(" ± ",round(SE_DT[[ni]],5))) 			#Agregado por DM
  }

  for(i in 1:length(DTime[[ni]])) if(DTime[[ni]][i] == SE_DT[[ni]][i]) {StEr3[[ni]][i]<-as.character(paste(" ± ",0.000001))}  	 	#Agregado por DM 

  datF2=data.frame(datF2,etiqueta=etiq3[[ni]],valor=paste(valor3[[ni]],StEr3[[ni]]),"|")			#Cambie 


###
#  if(abs(DTime[[ni]])<0.00001) DTime[[ni]]<-0.00001; if(SE_DT[[ni]]<0.00001) SE_DT[[ni]]<-0.00001; if(DTime[[ni]]==SE_DT[[ni]]) SE_DT[[ni]]<-0.000001
#  eti2[ni]=et2[[ni]]; val2[ni]=round(DTime[[ni]],5)
#  StE2[ni]=round(SE_DT[[ni]],5)
#################################Error in et2[[ni]] : subscript out of bounds


  for(i in 1:length(Mort[[ni]])) { if(abs(Mort[[ni]][i]) < 0.00001) {Mort[[ni]][i]<-0.00001}; if(SE_M[[ni]][i]<0.00001) {SE_M[[ni]][i]<-0.00001} }			#Agregado por DM

  if(nn2 > c2[ni])
  {
	valor2[[ni]]=c(round(Mort[[ni]],5),rep("",nn2-c2[ni]))
	etiq2[[ni]]=labels(valor2[[ni]])
	StEr2[[ni]]=as.character(c(paste(" ± ",round(SE_M[[ni]],5)),rep("",nn2-c2[[ni]]))) 			#Agregado por DM
  }
  else
  {
	valor2[[ni]]=as.character(round(Mort[[ni]],5))
	etiq2[[ni]]=labels(Mort[[ni]])
	StEr2[[ni]]=as.character(paste(" ± ",round(SE_M[[ni]],5))) 			#Agregado por DM
  }

  for(i in 1:length(Mort[[ni]])) if(Mort[[ni]][i] == SE_M[[ni]][i]) {StEr2[[ni]][i]<-as.character(paste(" ± ",0.000001))}  	 	#Agregado por DM

  datF3=data.frame(datF3,etiqueta=etiq2[[ni]],valor=paste(valor2[[ni]],StEr2[[ni]]),"|")			#Cambie 

 }
 
 datF1=datF1[,-1];datF2=datF2[,-1];datF3=datF3[,-1];	datF1=datF1[,-ncol(datF1)]; datF2=datF2[,-ncol(datF2)]; datF3=datF3[,-ncol(datF3)]

 labcol=c(1); for(i in 1:n1) labcol=c(labcol,paste(estadios[i],"label"),paste(estadios[i],"value"),"|"); labcol=labcol[-1]; labcol=labcol[-length(labcol)]                      #Cambie
# labcol=paste(rep(estadios[1:n1],rep(2,n1)),rep(c("label","value"),n1))
 colnames(datF1)=labcol;colnames(datF2)=labcol;colnames(datF3)=labcol

# rownames(datF1)=1:nn1;rownames(datF3)=1:nn2
# datF2=data.frame(Stage=estadios[1:n1],Link=eti2,Slope=paste(val2," ± ",StE2))				#Cambie

 TOvi=params$hfeno$ptazaeh_h; SE_TO=params$hfeno$StEr_ToviH
 for(i in 1:length(TOvi)) { if(abs(TOvi[i])<0.00001) TOvi[i]<-0.00001; if(SE_TO[i]<0.00001) SE_TO[i]<-0.00001; if(TOvi[i]==SE_TO[i]) SE_TO[i]<-0.000001} 
 eti4=labels(params$hfeno$ptazaeh_h); val4=round(TOvi,5); StE4=round(SE_TO,5)					#Agregado por DM
 datF4=data.frame(Label=eti4,Value=paste(val4," ± ",StE4)); rownames(datF4)=1:length(eti4)			#Agregado por DM
#datF4=data.frame(t(round(params$hfeno$ptazaeh_h,3)))

 ROvi=params$hfeno$povih_h; SE_RO=params$hfeno$StEr_RoviH
 for(i in 1:length(ROvi)) { if(abs(ROvi[i])<0.00001) ROvi[i]<-0.00001; if(SE_RO[i]<0.00001) SE_RO[i]<-0.00001; if(ROvi[i]==SE_RO[i]) SE_RO[i]<-0.000001} 
 eti5=labels(params$hfeno$povih_h); val5=round(ROvi,5); StE5=round(SE_RO,5)						#Agregado por DM
 datF5=data.frame(Label=eti5,Value=paste(val5," ± ",StE5)); rownames(datF5)=1:length(eti5)			#Agregado por DM


 spname <- paste("<b>Species name: </b>","<EM>",spProject,"</EM>")
 pname <- paste("<b>Project name: </b>",nameProject)
 autor <- paste("<b>Author name: </b>",autorName)
 sDate <- paste("<b>Compilation date: </b>",strDate)

 target <- HTMLInitFile(file.path(dirOut),filename=nameFile)
HTML.title("Parameter Values of the Functions Used", HR=1)
HTML(spname,file=target)
HTML(pname,file=target)
HTML(autor,file=target)
HTML(sDate,file=target)
HTML("<br></br>",file=target)

HTML.title("Development Rate", HR=3)
tmp00=datF1
HTML(tmp00,file=target)
HTML.title("Distribution function of development", HR=3)
tmp01=datF2
HTML(tmp01,file=target)
HTML.title("Mortality", HR=3)
tmp02=datF3
HTML(tmp02,file=target)
HTML.title("Total eggs per female", HR=3)
tmp03=datF4
HTML(tmp03,file=target)												#Cambie
HTML.title("Proportion of progeny production with age", HR=3)
tmp04=datF5
HTML(tmp04,file=target)												#Cambie
HTMLEndFile()

}

