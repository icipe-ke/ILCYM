#########################################################

# Uniendo de Indices por zonas separadas segun la latitud

#########################################################
#########################################################
zone.div<-function(dir1,dir2,ilon,ilat,R,dir.out,name.out=name.out,modelim=modelim,modelm=modelm,estadios=estadios,xi=xi,steps=steps,filtro=NULL,DL=NULL){


	lats=p.area(ilat,R)

	#############################################
	# Corriendo por cada area y almacenar su data
	#############################################

	#cl <- makeCluster(4)
	#registerDoParallel(cl)

	for(j in R:1){

		TempsAll=GenMatriz(dir1,dir2,ilon,lats[j,])

		coords=TempsAll$coords;nfil=nrow(coords)

		x1=TempsAll$x1

		y1=TempsAll$y1


		RFinal=matrix(NA,nfil,4)

		################################

		# Corriendo por punto del Area j

		################################

		#w <- startWorkers(4) # La PC tiene 4 nucleos
		#registerDoSMP(w)

		#system.time(

				foreach(i=1:nfil,.options.RNG = 123) %dorng% {

					RFinal[i,]=GenActIndex(i,TempsAll=TempsAll,coords=coords,x1=x1,y1=y1,modelim=modelim,modelm=modelm,estadios=estadios,xi=xi,steps=steps,filtro,DL)$indices
				}

		#)


		rm(TempsAll);rm(x1);rm(y1)

		Inds=data.frame(coords[1:nfil,],RFinal)

		rm(RFinal)

		##################################################################

		# Generando los archivos con encabezado: Lon - Lat - GI - AI - ERI

		##################################################################

		write.table(Inds,paste(dir.out,"file",j,".txt",sep=""),row.names = F)

		rm(Inds)

		#stopWorkers(w)
		#rmSessions(all.names=TRUE)

	}


	###########################################

	# Corriendo por cada area y uniendo su data

	###########################################

	Inds=read.table(paste(dir.out,"file",R,".txt",sep=""),header=T)

	if(R!=1)
	{
		for(j in (R-1):1){

			TempInds=read.table(paste(dir.out,"file",j,".txt",sep=""),header=T)
			Inds=rbind(Inds,TempInds)
		}
	}

	#mins<-apply(Inds[,-c(1,2)], 2,min, na.rm=T)
	#maxs<-apply(Inds[,-c(1,2)], 2,max, na.rm=T)

	rm(TempInds)
	gridded(Inds) = ~x+y ## Creando el objeto Grid

	writeAsciiGrid(Inds["X1"], na.value = -9999,paste(dir.out,"AI.asc",sep=""))

	writeAsciiGrid(Inds["X2"], na.value = -9999,paste(dir.out,"GI.asc",sep=""))

	writeAsciiGrid(Inds["X3"], na.value = -9999,paste(dir.out,"ERI.asc",sep=""))

	writeAsciiGrid(Inds["X4"], na.value = -9999,paste(dir.out,"Ro.asc",sep=""))
	
	### creacion d FLT(binario) ###
  	#writeGDAL(Inds["X3"], paste(dir.out,name.out,"_ERI.flt",sep=""), drivername = "EHdr")
	#writeGDAL(Inds["X1"], paste(dir.out,name.out,"_AI.flt",sep=""), drivername = "EHdr")
	#writeGDAL(Inds["X2"], paste(dir.out,name.out,"_GI.flt",sep=""), drivername = "EHdr")
  
	#return(list(mins=mins,maxs=maxs))
}

#############################################
########################################
# Funcion generadora de Indices Mejorado
########################################
GenActIndex<-function(posic,TempsAll=TempsAll,coords=coords,x1=x1,y1=y1,modelim=modelim,modelm=modelm,estadios=estadios,xi=xi,steps=steps,filtro=NULL,DL=NULL)
{
	d1<-1:length(x1);d2<-1:length(y1);filtroin=TRUE
	plon=d1[x1==coords[posic,1]]
	plat=d2[y1==coords[posic,2]]

	## ahora se extrae la temperatura desde listas en paralelo
	#Table=data.frame(mini=c(sapply(TempsAll$Tmin,extr1,plon,plat)),maxi=c(sapply(TempsAll$Tmax,extr1,plon,plat)))
	Table=data.frame(mini=c(sapply(TempsAll$Tmin,extr1,plon,plat)),maxi=c(sapply(TempsAll$Tmax,extr1,plon,plat))+3)

	Table=Table/10;
	Table2=cbind(id=1:nrow(Table),Table) ## filtrar en esta temperatura

	#print(posic)

	if(length(Table[is.na(Table)])==0)
	{
	   if(!is.null(filtro)){tmm=apply(Table,2,mean);if(tmm[1] > filtro[1] && tmm[2] < filtro[2]){filtroin=TRUE}else{filtroin=FALSE}}
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
			funcionm <- as.expression(hfeno$mortal[[K]][[3]])

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
		Ro=mean(TF[,"Ro"])

		GI=sum(1/TF[,"Gl"])
		temp.fr=TF[,"Fr"];ind1=temp.fr<1;temp.fr[ind1]=1
		AI=log(prod(temp.fr),10) 
		if(is.null(DL)){DL=0} ## La diferencia limite segun su origen siempre ha de ser positivo
		#ERI=1;for(s in 1:K){p0=(length(TF[TF[,s]<=DL,s]))/nmax;ERI=ERI*(1-p0)}


		ERI1=apply(TF[,1:K]*TF[,"Ro"],1,mean,na.rm = TRUE)
		#ERI=mean(ERI1)/quantile(ERI1,0.95);ERI[is.na(ERI)]=0;ERI[ERI>1]=1  ## Habria que hacer pruebas bajo distintas condiciones
		ERI=mean(ERI1)/max(ERI1);ERI[is.na(ERI)]=0;ERI[ERI>1]=1


		indices=c(AI=AI,GI=GI,ERI=ERI,Ro=Ro)
		return(list(indices=indices))
	   }else{indices=c(AI=0,GI=0,ERI=0,Ro=0);return(list(indices=indices))}
	}else{

		indices=c(AI=NA,GI=NA,ERI=NA,Ro=NA)
		return(list(indices=indices))
	}

}

##########################################################
##########################################################
extr1=function(mat,plon,plat){mat[plon,plat]}
#########################################################
GenMatriz<-function(dir1,dir2,ilon,ilat){

	archivos1=list.files(dir1,pattern="flt");archivos1=paste(dir1,"/",archivos1,sep="")
	archivos2=list.files(dir2,pattern="flt");archivos2=paste(dir2,"/",archivos2,sep="")
	
	Tmin1=readGDAL(archivos1[1])
	
	# para la extracion de los datos usamos la funcion:
	
	geodat=data.frame(coordinates(Tmin1))
	
	###########################################
	## Extraccion de las longitudes y latitudes
	
	x <- c(geodat[,1]);x <- unique(x)
	n1=length(x) ## tamaño
	
	y <- c(geodat[,2]); y <- unique(y)
	n2=length(y)
	
	#rm(geodat)
	
	#######################################################################
	## Extraccion de las resoluciones tanto para longitud como para latitud
	
	v1=Tmin1@grid@cellsize[1]/2
	v2=Tmin1@grid@cellsize[2]/2
	
	
	##############################################################
	## Definiendo los rangos en la longitud y latitud con posicion
	
	r11=ilon[1];r12=ilon[2];r21=ilat[1];r22=ilat[2]
	
	k1=rownames(geodat[geodat[,1] >= (r11-v1) & geodat[,1] <= (r12+v1),])
	k2=rownames(geodat[geodat[,2] >= (r21-v2) & geodat[,2] <= (r22+v2),])
	
	sector=intersect(k1,k2)
	sector=as.numeric(sector)
	
	coords=geodat[sector,]
	
	rm(geodat)
	
	
	#################################################
	## Definiendo los rangos en la longitud y latitud 
	
	ind1<-1:n1
	d1=x>=(r11-v1) & x<=(r12+v1)
	x1=x[d1]
	ind1=ind1[d1]
	
	
	ind2<-1:n2
	d2=y>=(r21-v2) & y<=(r22+v2)
	y1=y[d2]
	ind2=ind2[d2]
	
	
	####################################################################################
	## Creacion de la matriz que contiene los valores de la variable por cada coordenada
	
	
	Nps=length(archivos1)
	Tmin<-Tmax<-list(1)
	z <- c(Tmin1@data[,1]);zTmin1=matrix(z,n1,n2);rm(Tmin1)
	Tmin[[1]]=zTmin1[ind1,ind2]; rownames(Tmin[[1]])=x1; colnames(Tmin[[1]])=y1

	Tmax1=readGDAL(archivos2[1])
	z <- c(Tmax1@data[,1]);zTmax1=matrix(z,n1,n2);rm(Tmax1)
	Tmax[[1]]=zTmax1[ind1,ind2]; rownames(Tmax[[1]])=x1; colnames(Tmax[[1]])=y1


	for(i in 2:Nps) ## En caso de que no aguante la memoria se pueden generar archivos transfomados en matrices por cada flt
	{
		Tmin2=readGDAL(archivos1[i])
		z <- c(Tmin2@data[,1]);zTmin2=matrix(z,n1,n2);rm(Tmin2)
		Tmin[[i]]=zTmin2[ind1,ind2]; rownames(Tmin[[i]])=x1; colnames(Tmin[[i]])=y1
		rm(z)

		Tmax2=readGDAL(archivos2[i])
		z <- c(Tmax2@data[,1]);zTmax2=matrix(z,n1,n2);rm(Tmax2)
		Tmax[[i]]=zTmax2[ind1,ind2]; rownames(Tmax[[i]])=x1; colnames(Tmax[[i]])=y1
	}

	
	return(list(Tmin=Tmin,Tmax=Tmax,coords=coords,x1=x1,y1=y1))
}

##########################################
# Funcion generadora de Tasas y mortalidad
##########################################
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

######################################
######################################
matrizA <- function(estadios, hfeno, Table, steps){
	inmaduros <-  estadios[-(length(estadios)-1):-(length(estadios))]
	maduros   <-  estadios[(length(estadios)-1):(length(estadios))]
	nmax=nrow(Table)
	matriz<-matrix(0,ncol=2*(length(inmaduros)*2+length(maduros)+1),nrow=nrow(Table))
	Table2=cbind(id=1:nrow(Table),Table)
	
	for(K in 1:length(inmaduros))
	{
		#  Desarrollo  # extraccion de funciones y parametros
		
		parametrosc <- hfeno$pdv_dv[[K]]
		funciont <- as.expression(hfeno$fdv_dv[[K]][[3]])
		
		#   Mortalidad # extraccion de funciones y parametros
		
		parametrosm <- hfeno$pmortal[[K]]
		funcionm <- as.expression(hfeno$moratl[[K]][[3]])
		
		RM=apply(Table2,1,RateI,Table2,K,parametrosc,parametrosm,funciont,funcionm,nmax,steps) ## procesamiento de tasa de desarrollo y mortalidad por cada temperatura
		matriz[,4*K-3]=RM[1,];matriz[,4*K-2]=RM[2,]
		matriz[,4*K-1]=RM[3,];matriz[,4*K]=RM[4,]
	}
	#  Hembras
	
	parametrosc <- hfeno$pfh_h
	for (i in names(parametrosc)){temp <- parametrosc[i];storage.mode(temp) <- "double";assign(i, temp)}
	formulac <- hfeno$fh_h
	funciont <- as.expression(formulac[[3]])
	RM=apply(Table2,1,RateI,Table2,(K+1),parametrosc,funciont=funciont,nmax=nmax,steps=steps,J=NA)
	matriz[,4*(K+1)-3]=RM[1,];matriz[,4*(K+1)-2]=RM[2,]
	
	#  Machos
	
	parametrosc <- hfeno$pfm_m;
	for (i in names(parametrosc)){temp <- parametrosc[i];storage.mode(temp) <- "double";assign(i, temp)}
	formulac <- hfeno$fm_m
	funciont <- as.expression(formulac[[3]])
	RM=apply(Table2,1,RateI,Table2,(K+2),parametrosc,funciont=funciont,nmax=nmax,steps=steps,J=NA)
	matriz[,4*(K+1)-1]=RM[1,];matriz[,4*(K+1)]=RM[2,]
	matriz[matriz>1]<-1;matriz[matriz<0]<-0
	
	#  Fecundidad
	
	parametrosc <- hfeno$ptazaeh_h
	for (i in names(parametrosc)){temp <- parametrosc[i];storage.mode(temp) <- "double";assign(i, temp)}
	formulac <- hfeno$ftazaeh_h
	funciont <- as.expression(formulac[[3]])
	RM=apply(Table2,1,RateI,Table2,(K+3),parametrosc,funciont=funciont,nmax=nmax,steps=steps,J=3)
	matriz[,4*(K+1)+1]=RM[1,];matriz[,4*(K+1)+2]=RM[2,]
	
	vectorPares = (1:(ncol(matriz)/2))*2
	matrizA = matriz[,vectorPares]
	return(matrizA)
}
####################################################################
####################################################################
p.area<-function(ilat,R){
	R=R+1
	lats=seq(ilat[1],ilat[2],length.out=R)+0.0000000001
	mat.lat=matrix(NA,R-1,2)
	for(i in 2:R) 
	{
		mat.lat[i-1,]=c(lats[i-1],lats[i])
	}
	return(mat.lat)
}
######################################
######################################
lifeTableParameters <- function(estadios, hfeno, Table, steps){

	Table2=cbind(id=1:nrow(Table),Table) ## filtrar en esta temperatura

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
			funcionm <- as.expression(hfeno$mortal[[K]][[3]])

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
		return(list(indices=indices,TF=TF))
}
#############################################
#############################################

showParameterGraph <- function(TdV, parm, corry, grises){
	if(grises==TRUE){ccol="gray20"}else{ccol=2}	
	
	TdV[is.infinite(TdV)] = NA
	TdV[is.nan(TdV)] = NA
	laby = c("Generation length", "Net reproduction rate", "Intrinsic rate of increase", "Finite rate of increase", "Doubling time")
	plot(1:length(TdV),TdV,ylim=corry,xlim=c(1,length(TdV)),type="l",col=ccol,frame=F,axes=F,ylab=list(laby[parm], cex=0.9),xlab="", lwd=2)
	axis(2, at=seq(round(corry[1],2), round(corry[2],2), length.out=5), las=2, cex.axis=0.7)
	
	par(new=TRUE)
	mes<-1:length(TdV)
	plot(mes,TdV,xlab=list("Julian day",cex=0.8),ylab="",type="n",axes=F)
	#dias <- c(0,"","",100,"","",200,"","",300,"","")
	axis(1,at=cumsum(c(1,31,28,31,30,31,30,31,31,30,31,30)), cumsum(c(1,31,28,31,30,31,30,31,31,30,31,30)), cex.axis=0.7)
	par(new=TRUE)
} 
#############################################################################################################
graff_index<-function(TdV, Table1,grises)#,leg)
{
	TdV[is.infinite(TdV)] = NA#vefiricar a temperaturas bajas
	TdV[is.nan(TdV)] = NA
	if(grises==TRUE){ccol=c("gray40", "gray60", "gray20")}else{ccol=1:3}
	plot(1:length(TdV),TdV,ylim=c(0,max(TdV)),xlim=c(1,length(TdV)),type="l",col=ccol[3],frame=F,axes=F,ylab="",xlab="")
	axis(4)
	#legend(10.5,min(TdV),leg,col="blue",text.col="blue",lty=1)
	par(new=TRUE)
	mes<-1:length(TdV)
	plot(mes,TdV,xlab="month",ylab="",type="n",axes=F)
	meses<-c("J", "F", "M", "A", "M", "J", "J", "A" ,"S", "O" ,"N", "D")
	axis(1,at=cumsum(c(1,31,28,31,30,31,30,31,31,30,31,30)),lab=meses)
	par(new=TRUE)
	matplot(1:length(TdV),Table1[,1:2],type="l",col=ccol[2],frame=F,axes=F,ylab="T ºC",xlab="",lwd=2)
	axis(2,seq(round(min(Table1[,1])),max(round(Table1[,2])),length=6))  
	text(2.5,max(Table1[,1])+0.5,"min",col=ccol[1])
	text(2.5,max(Table1[,2])+0.5,"max",col=ccol[1])


}

###########################################
###########################################
matriz<-function(data)
{
	mes<-c(31,28,31,30,31,30,31,31,30,31,30,31)
	mat<-cbind(round(seq(data[,1][1],data[,1][2],length=mes[1]),1),round(seq(data[,2][1],data[,2][2],length=mes[1]),1))
	for(i in 2:nrow(data))
	{
		if(i != nrow(data)) mat<-rbind(mat,cbind(round(seq(data[,1][i],data[,1][i+1],length=mes[i]),1),round(seq(data[,2][i],data[,2][i+1],length=mes[i]),1)))
		else mat<-rbind(mat,cbind(round(seq(data[,1][i],data[,1][i-11],length=mes[i]),1),round(seq(data[,2][i],data[,2][i-11],length=mes[i]),1)))
	}
	return(list(temperaturas=mat))
}

###########################################
###########################################
extractTempsFromDB<-function(long, lat, dir){
	archivos1=list.files(dir1,pattern="flt");archivos1=paste(dir1,"/",archivos1,sep="")
	archivos2=list.files(dir2,pattern="flt");archivos2=paste(dir2,"/",archivos2,sep="")
	
	Tmin1=readGDAL(archivos1[1])
	
	# exttrae resolucion
	v1=Tmin1@grid@cellsize[1]/2
	v2=Tmin1@grid@cellsize[2]/2
	
	## Extraccion de las longitudes y latitudes
	geodat=data.frame(coordinates(Tmin1))
	x <- c(geodat[,1]);x <- unique(x)
	n1=length(x) ## tamaño
	
	y <- c(geodat[,2]); y <- unique(y)
	n2=length(y)
	rm(geodat)
	
	
	#se indica coordenadas
	r11=long;r12=long;r21=lat;r22=lat
	
	#extrae posicion en la matriz
	ind1<-1:n1
	d1=x>=(r11-v1) & x<=(r12+v1)
	x1=x[d1]
	ind1=ind1[d1]
	
	ind2<-1:n2
	d2=y>=(r21-v2) & y<=(r22+v2)
	y1=y[d2]
	ind2=ind2[d2]
	
	z <- c(Tmin1@data[,1]); zTmin1=matrix(z,n1,n2); rm(Tmin1)
	valueT1=zTmin1[ind1,ind2] / 10
	Tmin2=readGDAL(archivos1[2])
	z <- c(Tmin2@data[,1]); zTmin2=matrix(z,n1,n2); rm(Tmin2)
	valueT2=zTmin2[ind1,ind2] / 10
	Tmin3=readGDAL(archivos1[3])
	z <- c(Tmin3@data[,1]); zTmin3=matrix(z,n1,n2); rm(Tmin3)
	valueT3=zTmin3[ind1,ind2] / 10
	Tmin4=readGDAL(archivos1[4])
	z <- c(Tmin4@data[,1]); zTmin4=matrix(z,n1,n2); rm(Tmin4)
	valueT4=zTmin4[ind1,ind2] / 10
	Tmin5=readGDAL(archivos1[5])
	z <- c(Tmin5@data[,1]); zTmin5=matrix(z,n1,n2); rm(Tmin5)
	valueT5=zTmin5[ind1,ind2] / 10
	Tmin6=readGDAL(archivos1[6])
	z <- c(Tmin6@data[,1]); zTmin6=matrix(z,n1,n2); rm(Tmin6)
	valueT6=zTmin6[ind1,ind2] / 10
	Tmin7=readGDAL(archivos1[7])
	z <- c(Tmin7@data[,1]); zTmin7=matrix(z,n1,n2); rm(Tmin7)
	valueT7=zTmin7[ind1,ind2] / 10
	Tmin8=readGDAL(archivos1[8])
	z <- c(Tmin8@data[,1]); zTmin8=matrix(z,n1,n2); rm(Tmin8)
	valueT8=zTmin8[ind1,ind2] / 10
	Tmin9=readGDAL(archivos1[9])
	z <- c(Tmin9@data[,1]); zTmin9=matrix(z,n1,n2); rm(Tmin9)
	valueT9=zTmin9[ind1,ind2] / 10
	Tmin10=readGDAL(archivos1[10])
	z <- c(Tmin10@data[,1]); zTmin10=matrix(z,n1,n2); rm(Tmin10)
	valueT10=zTmin10[ind1,ind2] / 10
	Tmin11=readGDAL(archivos1[11])
	z <- c(Tmin11@data[,1]); zTmin11=matrix(z,n1,n2); rm(Tmin11)
	valueT11=zTmin11[ind1,ind2] / 10
	Tmin12=readGDAL(archivos1[12])
	z <- c(Tmin12@data[,1]); zTmin12=matrix(z,n1,n2); rm(Tmin12)
	valueT12=zTmin12[ind1,ind2] / 10
	
	Tmax1=readGDAL(archivos2[1])
	z <- c(Tmax1@data[,1]); zTmax1=matrix(z,n1,n2); rm(Tmax1)
	valueTm1=zTmax1[ind1,ind2] / 10
	Tmax2=readGDAL(archivos2[2])
	z <- c(Tmax2@data[,1]); zTmax2=matrix(z,n1,n2); rm(Tmax2)
	valueTm2=zTmax2[ind1,ind2] / 10
	Tmax3=readGDAL(archivos2[3])
	z <- c(Tmax3@data[,1]); zTmax3=matrix(z,n1,n2); rm(Tmax3)
	valueTm3=zTmax3[ind1,ind2] / 10
	Tmax4=readGDAL(archivos2[4])
	z <- c(Tmax4@data[,1]); zTmax4=matrix(z,n1,n2); rm(Tmax4)
	valueTm4=zTmax4[ind1,ind2] / 10
	Tmax5=readGDAL(archivos2[5])
	z <- c(Tmax5@data[,1]); zTmax5=matrix(z,n1,n2); rm(Tmax5)
	valueTm5=zTmax5[ind1,ind2] / 10
	Tmax6=readGDAL(archivos2[6])
	z <- c(Tmax6@data[,1]); zTmax6=matrix(z,n1,n2); rm(Tmax6)
	valueTm6=zTmax6[ind1,ind2] / 10
	Tmax7=readGDAL(archivos2[7])
	z <- c(Tmax7@data[,1]); zTmax7=matrix(z,n1,n2); rm(Tmax7)
	valueTm7=zTmax7[ind1,ind2] / 10
	Tmax8=readGDAL(archivos2[8])
	z <- c(Tmax8@data[,1]); zTmax8=matrix(z,n1,n2); rm(Tmax8)
	valueTm8=zTmax8[ind1,ind2] / 10
	Tmax9=readGDAL(archivos2[9])
	z <- c(Tmax9@data[,1]); zTmax9=matrix(z,n1,n2); rm(Tmax9)
	valueTm9=zTmax9[ind1,ind2] / 10
	Tmax10=readGDAL(archivos2[10])
	z <- c(Tmax10@data[,1]); zTmax10=matrix(z,n1,n2); rm(Tmax10)
	valueTm10=zTmax10[ind1,ind2] / 10
	Tmax11=readGDAL(archivos2[11])
	z <- c(Tmax11@data[,1]); zTmax11=matrix(z,n1,n2); rm(Tmax11)
	valueTm11=zTmax11[ind1,ind2] / 10
	Tmax12=readGDAL(archivos2[12])
	z <- c(Tmax12@data[,1]); zTmax12=matrix(z,n1,n2); rm(Tmax12)
	valueTm12=zTmax12[ind1,ind2] / 10
	
	TMinVal <- c(valueT1,valueT2,valueT3,valueT4,valueT5,valueT6,valueT7,valueT8,valueT9,valueT10,valueT11,valueT12)
	TMaxVal <- c(valueTm1,valueTm2,valueTm3,valueTm4,valueTm5,valueTm6,valueTm7,valueTm8,valueTm9,valueTm10,valueTm11,valueTm12)
	
	return(list(Tmin=TMinVal, Tmax=TMaxVal))
}
##########################################
###########################################
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
			funcionm <- as.expression(hfeno$mortal[[K]][[3]])

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
############################################
overlay.valid<-function(coords,ndias,dir1,dir2)
{
 archivos1=list.files(dir1,pattern="flt");archivos1=paste(dir1,"/",archivos1,sep="")
 archivos2=list.files(dir2,pattern="flt");archivos2=paste(dir2,"/",archivos2,sep="")
 ptm2=data.frame(lon=NA,lat=NA,loc=NA,tmin=NA,tmax=NA,dia=NA)

 for(i in ndias)
 {
	Tmin1=readGDAL(archivos1[i])
	Tmax1=readGDAL(archivos2[i])

	Temp0=data.frame(coords)

	Temp1=Temp0
	coordinates(Temp1) <- colnames(coordinates(Tmin1))
	Tmin1.interp <- overlay(Tmin1, Temp1) ## es un "SpatialPointsDataFrame" en ambos imputs
	Tmax1.interp <- overlay(Tmax1, Temp1) ## es un "SpatialPointsDataFrame" en ambos imputs

	yl1=Tmin1.interp$band1
	yl2=Tmax1.interp$band1
	Temp0$tmin=yl1
	Temp0$tmax=yl2
	Temp0$dia=rep(i,length(yl1))

	colnames(Temp0)=colnames(ptm2)
	ptm2=rbind(ptm2,Temp0)
 }
 return(ptm2=ptm2[-1,])
}
