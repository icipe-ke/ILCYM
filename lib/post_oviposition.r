##############################

# Reemplazo de un valor por NA

##############################



rm.dead<-function(datos,muerte)

{
	
	la <- as.matrix(datos);n1=nrow(datos);n2=ncol(datos)
	
	la[la==muerte] <- NA
	
	la=as.numeric(la)
	
	datos=matrix(la,n1,n2)
	
	return(datos)
	
}

##############################

# Tiempo de vida de un insecto

##############################



tmvida<-function(v1){length(v1[!is.na(v1)]-1)}





#####################

# Prueba chi cuadrado

#####################



chi.test.rate<-function(oviF,oviE)

{
	
	x2=sum(((oviF-oviE)^2)/oviE)
	
	pval=1-pchisq(x2,length(oviF)-1)
	
	Res=data.frame(X2=x2,P_value=pval)
	
	return(Res)
	
}





###########################

# Modelamiento de las Tasas

###########################





model.rate<-function(x,y,model)

{
	
	
	plot(x,y,xlim=c(0,max(x)+0.25*max(x)),xlab="Temperature in °C",ylab="Rate of oviposition",type="n")
	points(x,y)
	
	cols=c("royalblue","#FB6A4A","lightblue", "red", "green")
	
	
	
	modelos=c("Cuadrático","Exponencial 1","Exponencial 2","Taylor 1","Taylor 2")
	
	title("Curves in the rate of oviposition") 
	
	legend(max(x), .5, modelos, lty = 1,col=cols,cex=0.6) 
	
	
	
	modsel<-c(R2=NA,R2_Adj=NA,AIC=NA)
	
	
	
	for(i in 1:5)
	
	{
		
		modls=i
		
		## Aqui se definen las funciones de los modelos
		
		if(modls==1){f1=as.formula(y ~ b0 + b1*x + b2*x^2);inis0=coef(lm(y~x+I(x^2)));names(inis0)=NULL;inis=list(b0=inis0[1],b1=inis0[2],b2=inis0[3]);ff <- function(x){b0 + b1*x + b2*x^2}}
		
		if(modls==2){f1=as.formula(y ~ (exp(a+b*x))/(1+exp(a+b*x)));inis=list(a=-0.127141939,b=-0.004245373);ff <- function(x){(exp(a+b*x))/(1+exp(a+b*x))}}
		
		if(modls==3){f1=as.formula(y ~ 1-(exp(a+b*x))/(1+exp(a+b*x)));inis=list(a=-0.127141939,b=-0.004245373);ff <- function(x){1-(exp(a+b*x))/(1+exp(a+b*x))}}
		
		if(modls==4){f1=as.formula(y ~  c*exp((-0.5)*((x-u)^2)/(t^2)));inis=list(c=0.505,u=24.5,t=8);ff <- function(x){c*exp((-0.5)*((x-u)^2)/(t^2))}}
		
		#if(modls==5){f1=as.formula(y ~ 1-c*exp((-0.5)*((x-u)^2)/(t^2)));inis=list(c=0.505,u=24.5,t=8);ff <- function(x){1-c*exp((-0.5)*((x-u)^2)/(t^2))}} este modelo esta fallando..!!!
		
		
		
		## Se estiman los parametros usando el algoritmo de newton
		
		out <- nls(f1, start = inis,trace=F)
		
		if(model==i){outF=out}
		
		
		
		yl=fitted(out)
		art=0.00000000000000000000000000001
		sqe=sum(residuals(out)^2)+art
		
		sal <- coef(out)
		
		
		
		r<-1-sqe/sum((y-mean(y))^2) # este es el R^2
		
		r_ajus<- 1 - ((length(x) - 1) / (length(x) - length(sal))) * (1-r)  # este es el R^2 adjus
		
		AC<-AIC(out) # este es el AIC 
		
		
		
		modsel0<-data.frame(R2=round(r,3),R2_Adj=round(r_ajus,3),AIC=round(AC,3))
		
		modsel=rbind(modsel,modsel0)
		
		
		
		for (k in names(sal))
		
		{ 
			
			temp <- sal[k] 
			
			storage.mode(temp) <- "double"
			
			assign(k, temp)
			
		}
		
		curve(ff,add=TRUE,col=cols[i])  ## 2
		
	}
	
	modsel=modsel[-1,];rownames(modsel)=modelos
	
	return(list(outF=outF,modsel=modsel))
	
}





#######################################################################

# Control de calidad entre las tasas de oviposicion de hembras y machos

#######################################################################



q.c.rate<-function(dato1,dato2,muerte="Dead",model=1)

{
	
	dF=dim(dato1)
	
	dM=dim(dato1)
	
	
	
	#################################################################
	
	# viendo si hay diferente numero de insectos por cada temperatura
	
	
	
	if(dM[1]!=dF[1]){stop("the number of columns are different")}
	
	if(dM[2]!=dF[2]){stop("the number of rows are different")}
	
	
	
	#################################################################
	
	# viendo si hay diferente numero de insectos por cada temperatura
	
	
	
	datF=rm.dead(dato1,muerte) 
	
	datM=rm.dead(dato2,muerte) 
	
	
	
	tF=table(datF[,1])
	
	tM=table(datM[,1])
	
	
	
	ndif=length(tF[tF!=tM])
	
	if(ndif>=1){stop(paste("The number of repetitions in ",ndif," different temperatures, both females and males",sep=""))}
	
	
	
	#####################################################
	
	# viendo si un mismo insecto muere en dias diferentes
	
	
	
	tmvF=apply(datF,1,tmvida)
	
	tmvM=apply(datM,1,tmvida)
	
	
	
	ind=1:length(tmvF)
	
	#cat("\n","These insects are killed on different days, and being the same insect","\n","\n")
cat("\n","These insects are killed on different days, and being", "\n","the same insect","\n","\n")
	
	print(data.frame(Position=ind[tmvF!=tmvM],Temperature=datF[tmvF!=tmvM,1]))
	
	
	
	#######################################################################################
	
	# viendo la relacion de la oviposicion total por cada temperatura respecto a cada grupo
	
	
	
	temps=unique(datF[,1])
	
	oviM<-oviF<-c(1)
	
	for(i in 1:length(tF))
	
	{
		
		oviF[i]=sum(datF[datF[,1]==temps[i],],na.rm=T)
		
		oviM[i]=sum(datM[datM[,1]==temps[i],],na.rm=T)
		
	}
	
	
	
	Tasa=oviF/(oviF+oviM)
	
	
	
	######################################
	
	# caso1: prueba chi cuadrado tasa fija
	
	
	
	cat("\n","\n","Chi square test for a fixed rate of oviposition","\n")
	
	TasaT=sum(oviF)/sum(oviF+oviM)
	
	oviE=(oviF+oviM)*TasaT
	
	print(chi.test.rate(oviF,oviE))
	
	
	
	###############################################
	
	# caso1: prueba chi cuadrado en ajuste de tasas
	
	
	
	cat("\n","\n","Chi square test for a adjusted rate of oviposition","\n")
	
	x=temps
	
	y=Tasa
	
	
	
	for(j in 1:3)
	
	{
		
		out=model.rate(x,y,model=j)$outF
		
		TasaT=predict(out)
		
		
		
		cat("\n","** Using model",j,"\n")
		
		
		
		oviE=(oviF+oviM)*TasaT
		
		print(chi.test.rate(oviF,oviE))
		
	}
	
}





###############

# Desacumulador

###############



descum<-function(vec){dc<-c(1);dc[1]=vec[1];for(i in 2:length(vec)){dc[i]=vec[i]-vec[i-1]};return(dc)}


######################
tasas.ovip<-function(dato1,dato2,muerte) #,f1,f2,f3,f4
{
 datF=rm.dead(dato1,muerte) 
 datM=rm.dead(dato2,muerte)

 tmvF=apply(datF,1,tmvida)
 tmvM=apply(datM,1,tmvida)
 tmvF=data.frame(Temp=dato1[,1],Tvida=tmvF)
 tmvM=data.frame(Temp=dato1[,1],Tvida=tmvM)
 tm=by(tmvF[,2],tmvF[,1],median)


 temps=unique(datF[,1])
 TF=data.frame(Temp=NA,Day=NA,Tmort=NA,NormAge=NA,NeggF=NA,TeggF=NA,NeggM=NA,TeggM=NA,Negg=NA,Tegg=NA,TDeggF=NA,NormTDeggF=NA,NormTDeggF_Sim=NA,TDeggF_Sim=NA,TDegg=NA,TDegg_Des=NA,DPeggF=NA)

 for(i in 1:length(temps))
 {
    d1=datF[datF[,1]==temps[i],-1]
    d2=datM[datF[,1]==temps[i],-1]
    t0=(1:ncol(d1))/tm[i]
    t1=(nrow(d1)-apply(d1,2,tmvida))/nrow(d1)
    t21=apply(d1,2,sum,na.rm=T); t22=apply(d2,2,sum,na.rm=T); TT=t21+t22
    t31=sum(t21); t32=sum(t22) ; tT=t31+t32
    dp=t21/TT;dp[is.nan(dp)]=0
    

    f1=function(x){a*x^2+b*x+c};a=-0.00782385862094685;b=0.373335767277453;c=-3.94635389373034;te=f1(temps)
    f2=function(x){b1-b1/(1+b2*exp((-1)*(b3*x+b4*(x^2)+b5*(x^3))))};b1=1.12233427441224;b2=52.9518584742297;b3=1.78640099268727;b4=3.66422146033451;b5=-0.724881563625798;tsim=f2(t0)
    f3=function(x){1-exp((-1)*(b1*x+b2*(x^2)+b3*(x^3)))};b1=2.19238585109792;b2=-2.1998761139546;b3=5.21343140284426;tde=f3(t0);t4=descum(tde)
    f4=function(x){b1*1/x+b2*x+b3};b1=-4613.26455566904;b2=-16.0885801838611;b3=666.283125638406;tep=f4(temps)

    TF1=data.frame(Temp=temps[i],Day=1:ncol(d1),Tmort=t1,NormAge=t0,NeggF=t21,TeggF=cumsum(t21/t31),NeggM=t22,TeggM=cumsum(t22/t32),Negg=TT,Tegg=cumsum(TT/tT),TDeggF=dp,NormTDeggF=dp/te[i],NormTDeggF_Sim=tsim,TDeggF_Sim=tsim*te[i],TDegg=tde,TDegg_Des=t4,DPeggF=t4*tep[i])
    TF=rbind(TF,TF1)

 }
 TF=TF[-1,]

 NeggObs=TF[,9]/50;

 #t(TF[TF[,1]==temps[1],])
 return(list(temps=temps,TF=TF, datF=datF, tep=tep))
}



#####################################################################################################################
eval.dato<-function(dato,muerte){datoF=dato[dato[,2]!=muerte,];return(datoF)}
	
	