##########################################################################
# TITULO: METODO DE LA VALIDACION CRUZADA CON EL PRESS EN THIN PLATE #####
# datt:  Se refiere a la data original de regresión
# y  :  Columna que donde esta la variable dependiente
##########################################################################

# THIN PLATE
PRESS.TP=function(datt,y)
{

n=dim(datt)[1];p=dim(datt)[2]
datt=as.matrix(datt)

estim=rep(0,n)
nk=round(sqrt(nrow(datt)))-2

for (i in 1:n)
{
   model=Tps(datt[-i,-y],datt[-i,y])
   X=t(data.frame(datt[i,-y]))
   estim[i]= predict(model,X)
}
press=sum((y-estim)^2)
pressm=press/n
#cat("\n","VALIDACION CRUZADA PRESS","\n")
return(pressm)
}


##############################################################
## Interpolacion Individual y creacion del raster predecido ##
##############################################################

 TP.unique<-function(mat,X,Y,ag,name.ag,fileI="tpred.flt",CV=FALSE, m=2, Cons=1)
 {
   Temp=data.frame(mat[,X]) ## observados
   Temp2=data.frame(coordinates(ag),alt=ag@data[,c(name.ag)]) ## raster
   colnames(Temp)[1:2] <- colnames(coordinates(ag))

   # Thin plate
   tps <- Tps(Temp, mat[,Y],m)
   spl_pred <- Cons*predict(tps, Temp2[,colnames(Temp)]) ## verifica
   tempor=data.frame(coordinates(ag),Tpred=spl_pred)
   rm(spl_pred)

   if(CV){dat= data.frame(Temp,Y=mat[,Y]); CVi=PRESS.TP(dat,ncol(dat))}else{CVi=NA}

   gridded(tempor) = ~s1+s2 ## Creando el objeto Grid
   #writeGDAL(tempor["Tpred"],  fileI, drivername = "EHdr")
   writeAsciiGrid(tempor["Tpred"],na.value=-9999,  fileI)
   return(CVi)
 }

###########################
## Ploteo en Google Maps ##
###########################


 RasInt.PoinObs.GMaps<-function(Raster,MObs,color="Oranges",nivc=8, nameF='combination.htm', nameRaster)
 {
  proj4string(Raster) <- CRS("+init=epsg:4326")
  m1=plotGoogleMaps(Raster, zcol=nameRaster, mapTypeId='ROADMAP',colPalette=brewer.pal(nivc,color), add=TRUE)

  obs_point=MObs
  colnames(obs_point)[1:2] <- colnames(coordinates(Raster))
  coordinates(obs_point) <- colnames(coordinates(Raster))
  proj4string(obs_point) <- CRS("+init=epsg:4326")
  m<-plotGoogleMaps(obs_point, previousMap= m1, filename=nameF)
 }



##################################################################################         .
##  GENERADOR DE ARCHIVOS .flt AL APLICAR EL METODO DE INTERPOLACION THIN PLATE ##
##################################################################################


# daT   : data completa que contiene todas las variables y caracteristicas mencionadas abajo observadas
# clon  : un numero que determina la columna de longitudes
# clat  : un numero que determina la columna de latitudes
# calt  : un numero que determina la columna de la altura
# cnday : un numero que determina la columna de la posicion del dia
# ctmin : un numero que determina la columna de la temperatura minima
# ctmax : un numero que determina la columna de la temperatura maxima
# ag    : raster de la altura (tambien utilizamos su longitud y latitud como covariables de prediccion)
# ndias : un vector que indica los dias que se han de evaluar


TP.generator<-function(dat,clon,clat,calt,cnday,ctmin,ctmax,ag,ndias,dir1,dir2,CV=FALSE,m=2,Cons=10)
{
 name.ag=colnames(ag@data)[1]
 colnames(dat)[c(clon,clat)]=colnames(coordinates(ag))
 orden=c(paste("00",1:9,sep=""),paste("0",10:99,sep=""),100:366) ## tiene 366 por que es un maximo probable

 CVmin<-CVmax<-c(1)
 for(i in ndias)
 {
  ptm=dat[dat[,cnday]==i,] ## filtrando solo para el 1er dia
  X=c(clon,clat,calt)
  fileI1=paste(dir1,"/tmin",orden[i],".flt",sep="")
  fileI2=paste(dir2,"/tmax",orden[i],".flt",sep="")

  CVmin[i]=TP.unique(ptm,X,Y=ctmin,ag,name.ag,fileI=fileI1,CV=CV,m=m,Cons=Cons)
  CVmax[i]=TP.unique(ptm,X,Y=ctmax,ag,name.ag,fileI=fileI2,CV=CV,m=m,Cons=Cons)
 }
 return(list(CVmin=CVmin,CVmax=CVmax))
}



###########################################################
# Cálculo de la tasa de desarrollo senescencia y mortalidad

RateI<-function(vec,Table2,Ki,modelim,parametrosc,parametrosm=NULL,funciont,funcionm=NULL,nmax,steps,J=NA)
{
	for (i in names(parametrosc)){temp <- parametrosc[i];storage.mode(temp) <- "double";assign(i, temp)}
	i=0:(steps-1)
	T1=((vec[3]-vec[2])/2)*cos(pi*(i+0.5)/steps) + (vec[3]+vec[2])/2 ##  [ (Tmax_i -Tmin_i)/2 ]x Coseno(pi x (0.5) ) + [ (Tmax_i + Tmin_i)/2 ]
	if(vec[1]!=nmax){T2<-((vec[3]-Table2[vec[1]+1,2])/2)*cos(pi*(i+0.5)/steps) + (vec[3]+Table2[vec[1]+1,2])/2}else{  ##  [ (Tmax_i -Tmin_(i+1))/2 ]x Coseno(pi x (0.5) ) + [ (Tmax_i + Tmin_(i+1))/2 ]
		T2<-((vec[3]-Table2[1,2])/2)*cos(pi*(i+0.5)/steps) + (vec[3]+Table2[1,2])/2}

	if((modelim[Ki]==1 || modelim[Ki]==2 || modelim[Ki]==3 ||
				modelim[Ki]==4 || modelim[Ki]==5 || modelim[Ki]==6 ||
				modelim[Ki]==7 || modelim[Ki]==8 || modelim[Ki]==9 ||
				modelim[Ki]==10 || modelim[Ki]==11 || modelim[Ki]==12 || modelim[Ki]==13 ||
				modelim[Ki]==14) & is.na(J)){x = T1 + 273.15;x2 = T2 + 273.15}else{x = T1;x2 = T2}


	Rat=eval(funciont);
	#if(!is.na(J)){Rat[Rat<=0]=0}
	Rat[Rat<=0]=0
	Ratetot1<-(sum(Rat))/steps    ### aqui evalua la tasa de desarrollo de todas las divisiones del dia (i)
	x=x2;Rat=eval(funciont)
	#if(!is.na(J)){Rat[Rat<=0]=0}
	Rat[Rat<=0]=0
	Ratetot2<-(sum(Rat))/steps    ### aqui evalua la tasa de desarrollo de todas las divisiones del dia (i+1)
	Rate=(Ratetot1+Ratetot2)/2    ## aqui se calcula la tasa de desarrollo representativa del dia (i)

	## Si evaluamos temperatura mensual, entonces la temperatura T(i) = T(i+1), hasta que sea el ultimo dia del mes por lo que Ratetot1 es su T.D. representativa

	if(!is.null(funcionm))
	{
		for (i in names(parametrosm)){temp <- parametrosm[i];storage.mode(temp) <- "double";assign(i, temp)}
		x=T1;M1=eval(funcionm);M1[M1>1]=1;M1[M1<0]=0; Mortality1 <- (sum(M1))/steps    ### aqui evalua la Mortalidad de todas las divisiones
		x=T2;M2=eval(funcionm);M2[M2>1]=1;M2[M2<0]=0; Mortality2 <- (sum(M2))/steps    ### aqui evalua la Mortalidad de todas las divisiones del segundo dia
		Mortality=(Mortality1+Mortality2)/2
		return(c(Rate=Rate,Mortality=Mortality))
	}else{return(c(Rate=Rate))}
}



#####################################################################################
# Cálculos de los Indices observados (las temperaturas no estan multiplicadas por 10)

GenActIndex.uni<-function(Table,modelim=modelim,modelm=modelm,estadios=estadios,xi=xi,steps=steps,filtro=NULL,DL=NULL)
{
	Table2=cbind(id=1:nrow(Table),Table) ## filtrar en esta temperatura

	#print(posic)

	if(length(Table[is.na(Table)])==0)
	{
	   if(!is.null(filtro)){tmm=apply(Table,2,mean);if(tmm[1] > filtro[1] && tmm[2] < filtro[2]){filtroin=TRUE}else{filtroin=FALSE}}else{filtroin=TRUE}
	   if(filtroin)
	   {
		inmaduros <-  estadios[-(length(estadios)-1):-(length(estadios))]
		maduros   <-  estadios[(length(estadios)-1):(length(estadios))]
		nmax=nrow(Table)
		matriz<-matrix(0,ncol=length(inmaduros)*2+length(maduros)+1,nrow=nrow(Table))

		for(K in 1:length(inmaduros))
		{
			#  Desarrollo  # extraccion de funciones y parametros

			parametrosc <- hfeno$pdv_dv[[K]]
			funciont <- as.expression(hfeno$fdv_dv[[K]][[3]])

			#   Mortalidad # extraccion de funciones y parametros

			parametrosm <- hfeno$pmortal[[K]]
			funcionm <- as.expression(hfeno$moratl[[K]][[3]])

			RM=apply(Table2,1,RateI,Table2,K,modelim,parametrosc,parametrosm,funciont,funcionm,nmax,steps) ## procesamiento de tasa de desarrollo y mortalidad por cada temperatura

			matriz[,2*K-1]=RM[1,];matriz[,2*K]=RM[2,]
		}
		#  Hembras

		parametrosc <- hfeno$pfh_h
		for (i in names(parametrosc)){temp <- parametrosc[i];storage.mode(temp) <- "double";assign(i, temp)}
		formulac <- hfeno$fh_h
		funciont <- as.expression(formulac[[3]])
		RM=apply(Table2,1,RateI,Table2,(K+1),modelim,parametrosc,funciont=funciont,nmax=nmax,steps=steps,J=NA)
		matriz[,2*K+1]=RM

		#  Machos

		parametrosc <- hfeno$pfm_m;
		for (i in names(parametrosc)){temp <- parametrosc[i];storage.mode(temp) <- "double";assign(i, temp)}
		formulac <- hfeno$fm_m
		funciont <- as.expression(formulac[[3]])
		RM=apply(Table2,1,RateI,Table2,(K+2),modelim,parametrosc,funciont=funciont,nmax=nmax,steps=steps,J=NA)
		matriz[,2*(K+1)]=RM

		#  Fecundidad

		parametrosc <- hfeno$ptazaeh_h
		for (i in names(parametrosc)){temp <- parametrosc[i];storage.mode(temp) <- "double";assign(i, temp)}
		formulac <- hfeno$ftazaeh_h
		funciont <- as.expression(formulac[[3]])
		RM=apply(Table2,1,RateI,Table2,(K+3),modelim,parametrosc,funciont=funciont,nmax=nmax,steps=steps,J=3)
		matriz[,2*(K+1)+1]=RM

		TF=data.frame(Ind=1:nmax)
		for(i in 1:K){TF1=data.frame(DaiSurv1=(1-matriz[,2*i])^(matriz[,2*i-1]+1e-14));TF=cbind(TF,TF1)}  ## al aumentar el decimal 1e14 me dará (1-1)^(1e14)=0 , de lo contrario me daria (1-1)^0=1
		for(i in 1:K){Gl0=data.frame(Gl01=1/matriz[,2*i-1]);TF=cbind(TF,Gl0)}
		TF$Gl11=(1/matriz[,2*K+1])/(1/xi);TF=TF[,-1]
		temp1=TF[,(K+1):ncol(TF)]
		TF$Gl=apply(temp1,1,sum)
		temp1<-temp2<-1;for(j in 1:K){temp1=temp1*(1-matriz[,2*j])}
		TF$Is=temp1
		TF$Ro=(TF$Is*matriz[,2*K+3])/2;TF$Ro[TF$Ro<0]=0
		TF$rm=log(TF$Ro)/TF$Gl
		TF$Fr=exp(TF$rm)
		TF$Dt=log(2)/TF$rm

		GI=sum(1/TF[,"Gl"])
		temp.fr=TF[,"Fr"];ind1=temp.fr<1;temp.fr[ind1]=1
		AI=log(prod(temp.fr),10)
		if(is.null(DL)){DL=0} ## La diferencia limite segun su origen siempre ha de ser positivo
		#ERI=1;for(s in 1:K){p0=(length(TF[TF[,s]<=DL,s]))/nmax;ERI=ERI*(1-p0)}

		ERI1=apply(TF[,1:K]*TF[,"Ro"],1,mean,na.rm = TRUE)
		#ERI=mean(ERI1)/quantile(ERI1,0.95);ERI[is.na(ERI)]=0;ERI[ERI>1]=1  ## Habria que hacer pruebas bajo distintas condiciones
		ERI=mean(ERI1)/max(ERI1);ERI[is.na(ERI)]=0;ERI[ERI>1]=1

		indices=c(AI=AI,GI=GI,ERI=ERI)
		return(list(indices=indices))
	   }else{indices=c(AI=0,GI=0,ERI=0);return(list(indices=indices))}
	}else{

		indices=c(AI=NA,GI=NA,ERI=NA)
		return(list(indices=indices))
	}

}


#########################################################################
# Extraccion y validacion de interpolaciones segun los valores observados

overlay.valid<-function(ptm1,clon,clat,calt,cnday,ctmin,ctmax,ccod,ndias,dir1,dir2,Map1=FALSE,Div=1)
{

 archivos1=list.files(dir1,pattern="flt");archivos1=paste(dir1,"/",archivos1,sep="")
 archivos2=list.files(dir2,pattern="flt");archivos2=paste(dir2,"/",archivos2,sep="")
 sqe1<-sqe2<-r1<-r2<-rep(NA,length(ndias))
 j=1
 ptm2=data.frame(x=NA,x=NA,x=NA,x=NA,x=NA,x=NA,x=NA,x=NA,x=NA)
 for(i in ndias)
 {
	Tmin1=readGDAL(archivos1[j])
	Tmax1=readGDAL(archivos2[j])
	ptm=ptm1[ptm1[,cnday]==i,]
 	colnames(ptm)[c(clon,clat)]=colnames(coordinates(Tmin1))
	Temp0=data.frame(ptm[,c(clon,clat,calt,ctmin,ctmax,ccod)]) ## observados (6 y 7: tienen la columna de la temp min y max)

	Temp1=Temp0
	coordinates(Temp1) <- colnames(coordinates(Tmin1))
	Tmin1.interp <- overlay(Tmin1, Temp1) ## es un "SpatialPointsDataFrame" en ambos imputs
	Tmax1.interp <- overlay(Tmax1, Temp1) ## es un "SpatialPointsDataFrame" en ambos imputs

	yo1=Temp0[,4] # tmin
	yo2=Temp0[,5] # tmax
	yl1=Tmin1.interp$band1/Div
	yl2=Tmax1.interp$band1/Div
	sqe1[j]=sum((yl1-yo1)^2,na.rm = TRUE)
	sqe2[j]=sum((yl2-yo2)^2,na.rm = TRUE)
	r1[j]<-1-sqe1[j]/sum((yo1-mean(yo1))^2,na.rm = TRUE) # este es el R^2
	r2[j]<-1-sqe2[j]/sum((yo2-mean(yo2))^2,na.rm = TRUE) # este es el R^2
	Temp0$tmin.pred=yl1
	Temp0$tmax.pred=yl2
	Temp0$dia=rep(i,length(yo1))

	colnames(ptm2)=colnames(Temp0)
	ptm2=rbind(ptm2,Temp0)

	if(Map1){
		RasInt.PoinObs.GMaps(Raster=Tmin1,MObs=Temp0,nivc=8,color="BrBG")
		RasInt.PoinObs.GMaps(Raster=Tmax1,MObs=Temp0,nivc=8,color="BrBG")
	}
 j=j+1
 }
 RR=data.frame(SSR_Tmin=sqe1,R2_Tmin=r1,SSR_Tmax=sqe2,R2_Tmax=r2);row.names(RR)=ndias
 return(list(RR=RR,ptm2=ptm2[-1,]))
}




#########################################
# Graficos de las interpolaciones diarias

graphs.interps<-function(tempor,direc,TG=FALSE)
{
 colnames(tempor)[c(4,5)]=c("tmin","tmax")

 ### Tmin
 tempor1=data.frame(tempor[,-5],Temp=rep("Tmin",nrow(tempor)));colnames(tempor1)[4]="Tmin"
 tempor2=data.frame(tempor[,-4],Temp=rep("Tmmin.Pred",nrow(tempor)));colnames(tempor2)[4]="Tmin";tempor2[,4]=tempor[,7]
 temporF=rbind(tempor1,tempor2)

 if(TG)
 {
	 p = ggplot(data=temporF, aes(x=dia)) + geom_line(aes(y=Tmin, group=Temp, color=Temp)) + facet_wrap(~Codigo, scales='free')
	 dev.new(width=8, height=4)
	 p + geom_ribbon(aes(ymin=tmin, ymax=tmin.pred),  alpha=0.3, data=tempor)
 }

 nombres=unique(tempor[,6])

 for(i in 1:length(nombres))
 {
  temporF0=tempor[tempor[,6]==nombres[i],]
  temporF1=temporF[temporF[,5]==nombres[i],]
  xvalues <- temporF1[temporF1[,9]=="Tmin",4]
  yvalues <- temporF1[temporF1[,9]=="Tmmin.Pred",4]
  R2=(summary(lm(yvalues~xvalues)))[8]
  d <- data.frame(xvalues,yvalues)

  dev.new(width=5, height=5)
  p=qplot(xvalues,yvalues)+geom_smooth(method="lm") + ggtitle(paste("R Squared =",round(R2[[1]],4)))
  ggsave(filename=paste(direc,"/Tmin - ",nombres[i],".jpeg",sep=""), plot=p)
  rm(p)
  dev.off()

  dev.new(width=12, height=4)
  p = ggplot(data=temporF1, aes(x=dia)) + geom_line(aes(y=Tmin, group=Temp, color=Temp)) + facet_wrap(~Codigo, scales='free')  + geom_ribbon(aes(ymin=tmin, ymax=tmin.pred),  alpha=0.3, data=temporF0)
  ggsave(filename=paste(direc,"/Tmin - Time series- ",nombres[i],".jpeg",sep=""), plot=p)
  rm(p)
  dev.off()
 }



 ### Tmax
 tempor1=data.frame(tempor[,-4],Temp=rep("Tmax",nrow(tempor)));colnames(tempor1)[4]="Tmax"
 tempor2=data.frame(tempor[,-5],Temp=rep("Tmax.Pred",nrow(tempor)));colnames(tempor2)[4]="Tmax";tempor2[,4]=tempor[,8]
 temporF=rbind(tempor1,tempor2)

 if(TG)
 {
	 p = ggplot(data=temporF, aes(x=dia)) + geom_line(aes(y=Tmax, group=Temp, color=Temp)) + facet_wrap(~Codigo, scales='free')
	 dev.new(width=5, height=4)
	 p + geom_ribbon(aes(ymin=tmax, ymax=tmax.pred),  alpha=0.3, data=tempor)
 }

 for(i in 1:length(nombres))
 {
  temporF0=tempor[tempor[,6]==nombres[i],]
  temporF1=temporF[temporF[,5]==nombres[i],]
  xvalues <- temporF1[temporF1[,9]=="Tmax",4]
  yvalues <- temporF1[temporF1[,9]=="Tmax.Pred",4]
  R2=(summary(lm(yvalues~xvalues)))[8]
  d <- data.frame(xvalues,yvalues)

  dev.new(width=5, height=5)
  p=qplot(xvalues,yvalues)+geom_smooth(method="lm") + ggtitle(paste("R Squared =",round(R2[[1]],4)))
  ggsave(filename=paste(direc,"/Tmax - ",nombres[i],".jpeg",sep=""), plot=p)
  rm(p)
  dev.off()

  dev.new(width=12, height=4)
  p = ggplot(data=temporF1, aes(x=dia)) + geom_line(aes(y=Tmax, group=Temp, color=Temp)) + facet_wrap(~Codigo, scales='free')  + geom_ribbon(aes(ymin=tmax, ymax=tmax.pred),  alpha=0.3, data=temporF0)
  ggsave(filename=paste(direc,"/Tmax - Time series - ",nombres[i],".jpeg",sep=""), plot=p)
  rm(p)
  dev.off()
 }

}



#####################################################################
# Interpolacion de Indices y cálculo con otra temperatura interpolada

# ptm     : tabla que ha sido generada previamente con la funcion "overlay.valid()"
# modelim : caracteristicas de la fenologia
# modelm  :
# estadios:
# xi      :
# steps   :
# m       : grado de la interpolacion en el spline (por defecto es 3)
maps.indic<-function(ptm,modelim=modelim,modelm=modelm,estadios=estadios,xi=xi,steps=steps,m=2,filedem=filedem,INDpath=INDpath)
{
 ## Usando las temperaturas observadas
 Tables<-list(1)
 loc=unique(ptm[,6])
 for(i in 1:length(unique(ptm[,6]))){Tables[[i]]=ptm[ptm[,6]==loc[i],4:5]}

 ptm1= aggregate(ptm[,3],list(x=ptm[,1],y=ptm[,2],loc=ptm[,6]),mean)

 indic2=sapply(Tables,GenActIndex.uni,modelim=modelim,modelm=modelm,estadios=estadios,xi=xi,steps=steps)
 IND1=data.frame(ptm1[,-3],loc=loc,t(data.frame(indic2)));rownames(IND1)=1:nrow(IND1);colnames(IND1)[3]="alt"

 write.table(IND1,paste(INDpath,"Indices_Temp_Observ.txt",sep=""))

 ## Usando las temperaturas interpoladas
 #Tables=list(1)
 #loc=unique(ptm[,6])
 #for(i in 1:length(unique(ptm[,6]))){Tables[[i]]=ptm[ptm[,6]==loc[i],7:8]} ## 6:7 indica la Tmin y Tmax observados

 #indic2=sapply(Tables,GenActIndex.uni,modelim=modelim,modelm=modelm,estadios=estadios,xi=xi,steps=steps)
 #IND2=data.frame(ptm1[,-3],loc=loc,t(data.frame(indic2)));rownames(IND2)=1:nrow(IND2);colnames(IND2)[3]="alt"

 #write.table(IND2,paste(INDpath,"Indices_de_Temp_Interp.txt",sep=""))



 ##############################################################################################

 ### Para el AI
 X=c(clon=1,clat=2,calt=3)
 Y=5 # la temperatura minima
 ag = readAsciiGrid(filedem) ## raster altura 
 name.ag=colnames(ag@data)[1]
 fileI=paste(INDpath,"AI",1,".asc",sep="")

 TP.unique(IND1,X,Y,ag,name.ag,fileI,CV=TRUE,m=m)

 AI=readGDAL(fileI)
 AI$band1[AI$band1<0]=0
 writeAsciiGrid(AI["band1"],na.value=-9999,  fileI)
 #spplot(AI[c("band1")])
 #Temp01=data.frame(IND1[,c(X,Y)],AI.Temp_Intp=IND2[,Y],Location=IND2[,4]) ## observados
 Temp01=data.frame(IND1[,c(X,Y)],Location=IND1[,4]) ## observados

 Temp1=Temp01
 coordinates(Temp1) <- colnames(coordinates(AI)) 
 AI.interp <- overlay(AI, Temp1)
 Temp01$AI.Index.Interp=AI.interp$band1


 ### Para el GI

 X=c(clon=1,clat=2,calt=3)
 Y=6 # la temperatura minima
 name.ag=colnames(ag@data)[1]
 fileI=paste(INDpath,"GI",1,".asc",sep="")

 TP.unique(IND1,X,Y,ag,name.ag,fileI,CV=TRUE,m=m)

 GI=readGDAL(fileI)
 GI$band1[GI$band1<0]=0
 writeAsciiGrid(GI["band1"],na.value=-9999,  fileI)
 #spplot(GI[c("band1")])
 #Temp02=data.frame(IND1[,c(X,Y)],GI.Temp_Intp=IND2[,Y],Location=IND2[,4]) ## observados
 Temp02=data.frame(IND1[,c(X,Y)],Location=IND1[,4]) ## observados


 Temp1=Temp02
 coordinates(Temp1) <- colnames(coordinates(GI)) 
 GI.interp <- overlay(GI, Temp1)
 Temp02$GI.Index.Interp=GI.interp$band1


 ### Para el ERI

 X=c(clon=1,clat=2,calt=3)
 Y=7 # la temperatura minima
 name.ag=colnames(ag@data)[1]
 fileI=paste(INDpath,"ERI",1,".asc",sep="")

 TP.unique(IND1,X,Y,ag,name.ag,fileI,CV=TRUE,m=m)

 ERI=readGDAL(fileI)
 ERI$band1[ERI$band1<0]=0
 writeAsciiGrid(ERI["band1"],na.value=-9999,  fileI)
 #spplot(ERI[c("band1")])
 #Temp03=data.frame(IND1[,c(X,Y)],ERI.Temp_Intp=IND2[,Y],Location=IND2[,4]) ## observados
 Temp03=data.frame(IND1[,c(X,Y)],Location=IND1[,4]) ## observados

 Temp1=Temp03
 coordinates(Temp1) <- colnames(coordinates(ERI)) 
 ERI.interp <- overlay(ERI, Temp1)
 Temp03$ERI.Index.Interp=ERI.interp$band1

 return(list(Temp.AI=Temp01, Temp.GI=Temp02, Temp.ERI=Temp03))

}
