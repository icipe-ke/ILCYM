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
#############################################################################################################
#############################################################################################################
life.table <- function (data, estad){
    maduros   <-  estad[(length(estad)-1):(length(estad))]
    inmaduros <-  estad[-(length(estad)-1):-(length(estad))]
    estadios <- estad
    datat <- as.matrix(data)
    num <- matrix(as.numeric(datat), nrow(datat), ncol(datat))
    numc <- num
    for (i in 1:nrow(numc)) for (j in 1:ncol(numc)) if (is.na(numc[i, j])) numc[i, j] <- 0
    newegg <- matrix(apply(numc, 1, sum), nrow(numc))
    day <- 0:(nrow(num) - 1)
    oviposition <- apply(numc, 2, sum)
    estados1 <- matrix(rep(0, length(day)), nrow = length(day),
                                            ncol = length(estadios))
    #for (i in 1:day) for (j in 1:length(estadios)) estados1[, j] <- Estado(data, estad, estadios[j])
	for (j in 1:length(estadios)) estados1[, j] <- Estado(data, estad, estadios[j]) #verificar
    life.table <- data.frame(estados1, newegg)
    nombres <- rep(0, length(estadios))
    for (i in 1:length(estadios)) nombres[i] <- estadios[i]
    colnames(life.table) <- c(nombres, "new.Egg")
    trues <- rep(TRUE, nrow(life.table))
    for (i in 1:nrow(life.table)) trues[i] <- sum(life.table[i, ]) != 0
    life.table <- life.table[trues, ]
    print(life.table)
    return(list(life.table=life.table))
}
#############################################################################################################
#############################################################################################################
                                                                     
                                                                     
                                                                     
                                             
statist <- function (data, estad,poli=1){
    estadios <- estad
    esadult <- estadios[(length(estadios) - 1):length(estadios)]
    datat <- as.matrix(data)
    estadoa <- as.list(as.data.frame((datat)))
    estadob <- as.list(rep(0, ncol(data)))
    num <- matrix(as.numeric(datat), nrow(datat), ncol(datat))
    numc <- num
    for (i in 1:nrow(numc)) for (j in 1:ncol(numc)) if (is.na(numc[i, j]))
                                                    numc[i, j] <- 0
    newegg <- matrix(apply(numc, 1, sum), nrow(numc))
    est1 <- subset(estadios, estadios != esadult[2])
    nopupa1 <- nopupa2 <- rep(0, (length(est1) - 2))
    estados3 <- huevo <- matrix(rep(0, length(nopupa1) * ncol(data)),
                         nrow = ncol(data), ncol = length(nopupa1))
    for (is in 1:(length(est1) - 1)) {
        est2 <- est1[is + 1]
        ifelse(est2 == esadult[1], caso <- "caso2", caso <- "caso1")
        if (caso == "caso1") {
            for (i in 1:ncol(data)) {
                estadoa[[i]] <- as.matrix(estadoa[[i]])
                estadob[[i]] <- length(subset(estadoa[[i]], estadoa[[i]] ==  est2))
            }
            estados3[, is] <- as.numeric(as.character(as.matrix(estadob)))
            for (i in 1:length(estados3[, 1])) for (j in 1:length(nopupa1)) if (estados3[i, j] != 0)
                                                                            estados3[i, j] <- 1
            nopupa1[is] <- sum(estados3[, is])
            huevo[, is] <- Estado(t(data), estad, est1[is])
            for (i in 1:length(estados3[, 1])) for (j in 1:length(nopupa1)) ifelse(estados3[i, j] == 1,
                                                                            huevo[i, j] <- huevo[i, j], huevo[i, j] <- 0)
            for (i in 1:length(estados3[, 1])) for (j in 1:length(nopupa1)) ifelse(huevo[i, j] == 0,
                                                                            huevo[i, j] <- NA, huevo[i, j] <- huevo[i, j])
            nopupa2[is] <- mean(huevo[, is], na.rm = T)
        }
    }
    if (caso == "caso2") {
        for (i in 1:ncol(data)) {
            estadoa[[i]] <- as.matrix(estadoa[[i]])
            estadob[[i]] <- length(subset(estadoa[[i]], estadoa[[i]] == esadult[2]))
        }
        estados2 <- as.numeric(as.character(as.matrix(estadob)))
        for (i in 1:length(estados2)) if (estados2[i] != 0) estados2[i] <- 1
        male <- sum(estados2)
        num <- matrix(as.numeric(datat), nrow(datat), ncol(datat))
        a <- as.list(as.data.frame((num)))
        b <- as.list(rep(0, ncol(data)))
        for (i in 1:ncol(data)) b[[i]] <- sum(as.matrix(table(a[[i]])))
        female1 <- as.numeric(as.character(as.matrix(b)))
        for (i in 1:length(female1)) if (female1[i] != 0) female1[i] <- 1
        female <- sum(female1)
        for (i in 1:ncol(data)) {
            estadoa[[i]] <- as.matrix(estadoa[[i]])
            estadob[[i]] <- length(subset(estadoa[[i]], estadoa[[i]] == (estad[is + 3])))
        }
        estados5 <- as.numeric(as.character(as.matrix(estadob)))
        for (i in 1:length(estados5)) if (estados5[i] != 0) estados5[i] <- 1
        dead <- sum(estados5)
        pupa1 <- male + female
        capullo <- Estado(t(data), estad, est1[is])
        capullo2 <- rep(0, length(female1))
        for (i in 1:length(female1)) ifelse(estados2[i] == 0 &
            female1[i] == 0, capullo2[i] <- 0, capullo2[i] <- 1)
        for (i in 1:length(female1)) ifelse(capullo2[i] == 1,
            capullo[i] <- capullo[i], capullo[i] <- 0)
        for (i in 1:length(female1)) ifelse(capullo[i] == 0,
            capullo[i] <- NA, capullo[i] <- capullo[i])
        pupa <- mean(capullo, na.rm = T)
    }

    f <- c(ncol(data), female/(female + male), poli*male, poli*female,
        ncol(data) - male - female, male + female + dead, sum(newegg)/female)
    observations <- data.frame(f)

    if(poli>1){rownames(observations) <- c("Number of hosts parasitized :", "Sex ratio of parasitoid:","Males of parasitoid :", "Females of parasitoid :", "Immature death of Host:", "", "Total_Eggs/Total_Females  of parasitoid :")}else{
    rownames(observations) <- c("Number of insects  :", "Sex ratio :","Males :", "Females :", "Immature death :", "", "Total_Eggs/Total_Females:")}

    colnames(observations) <- c("Observations")
    todo <- rep(0, length(female1))
    for (i in 1:length(female1)) todo[i] <- sum(capullo[i], huevo[i, ])
    todos <- mean(todo, na.rm = T)
    time <- c(round(nopupa2, 3), round(pupa, 3), round(todos, 3))
    insect <- c(nopupa1, pupa1, "")
    d1 <- rep(0, (length(nopupa2) - 1))
    for (i in 1:(length(nopupa2) - 1)) d1[i] <- nopupa1[i + 1] * 100/nopupa1[i]
    d2 <- pupa1 * 100/nopupa1[-1:-(length(nopupa2) - 1)]
    percente <- c(round((nopupa1[1]/(ncol(data))) * 100, 2), round(d1, 2), round(d2, 2))/100
    percent <- c(paste(percente * 100, "%"), "")
    acc <- rep(1, length(nopupa1) + 1)
    acc[1] <- percente[1]
    for (i in 2:(length(nopupa1) + 1)) acc[i] <- acc[i - 1] * percente[i]
    accume <- acc
    accum <- c(paste(round(100 * accume, 3), "%"), "")
    survival <- data.frame(cbind(time, insect, percent, accum))
    rownames(survival) <- c(est1[-(length(est1))], "Total")
    colnames(survival) <- c("Time", "Insect", "Percent", "Accumm")
    inmor <- rep(0, (length(nopupa1) - 1))
    for (i in 1:(length(nopupa1) - 1)) inmor[i] <- nopupa1[i] - nopupa1[i + 1]
    insectmor <- c(ncol(data) - nopupa1[1], inmor, nopupa1[-1:-(length(nopupa2) - 1)] - pupa1)
    pmor <- rep(0, (length(nopupa1) - 1))
    for (i in 1:(length(nopupa1) - 1)) pmor[i] <- (nopupa1[i] -
        nopupa1[i + 1]) * 100/nopupa1[i]
    permore <- c(round(((ncol(data) - nopupa1[1])/(ncol(data))) *
        100, 2), round(pmor, 2), round((nopupa1[-1:-(length(nopupa2) - 1)] - pupa1)
        * 100/nopupa1[-1:-(length(nopupa2) - 1)], 2))
    permor <- paste(round(permore, 3), "%")
    acc <- c(paste(round((1 - accume) * 100, 3), "%"))
    mortality <- data.frame(cbind(insectmor, permor, acc))
    rownames(mortality) <- est1[-(length(est1))]
    colnames(mortality) <- c("Insect", "Percent", "Accumm")
    print(observations)
    if(poli==1){cat("\n","Survival","\n")}else{cat("\n","Survival of Host","\n")}
    print(survival)
    if(poli==1){cat("\n","Mortality","\n")}else{cat("\n","Mortality of Host","\n")}
    print(mortality)
    
    if(poli>1){salida <- list(Observ = observations, Survival = survival, Mortality = mortality)
    names(salida)=c("Observ","Survival of Host","Mortality of Host")
    }else{salida <- list(Observ = observations, Survival = survival, Mortality = mortality)}

    return(salida)
}
#############################################################################################################
#############################################################################################################
parameters <- function (N, estad, ltb,poli=1){
 if(nrow(ltb)!=0)
 {
    estadios <- estad
    lifetable <- ltb
    s.x <- m.x <- l.x <- rep(0, nrow(lifetable))
    for (i in 1:nrow(lifetable)) l.x[i] <- sum(lifetable[i, 1:(length(estadios))])/N  ## porcentaje de vivos por dias
    # for (i in 1:nrow(lifetable)) if(l.x[i]==0) m.x[i] <- 0 else m.x[i] <- lifetable[,(length(estad)+1)][i]/N/l.x[i]
    for (i in 1:nrow(lifetable)) if(l.x[i]==0) m.x[i] <- 0 else m.x[i] <- poli*lifetable[,(length(estad)+1)][i]/N/l.x[i]
    ## el objeto m.x no es afectado por el efecto poliembrionico pero asumimos la mortalidad del host es la unica que afecta ala mortalidad del parasitoide
    s.x1 <- c(l.x[1], l.x[2])
    s.x2 <- rep(0, nrow(lifetable) - 2)
    for (i in 1:(nrow(lifetable) - 2)) s.x2[i] <- l.x[i + 2]/l.x[i + 1]
    s.x <- c(s.x1, s.x2) ## es el porcentaje de vivos segun la cantidad del dia anterior
    Ro <- sum(l.x * m.x)  ## este parametro depende del procentaje de vivos de los inmaduros
    Day <- 0:(nrow(lifetable)-1)
    Tg <- sum(Day * l.x * m.x)/Ro
    GRR <- sum(m.x)
    r <- log(Ro)/Tg
    euler <- sum(l.x * m.x * exp(-r * (Day + 1)))
    r1 <- r - 1e-06
    for (k in 1:100) {
        euler1 <- sum(l.x * m.x * exp(-r1 * (Day + 1)))
        r2 <- ((1 - euler) * (r1 - r)/(euler1 - euler)) + r
        r1 <- r2
    }
    rf <- r1
    Dt <- log(2)/rf
    lamda <- exp(rf)
    T <- log(Ro)/rf
    parameters <- t(data.frame(r = rf, Ro = Ro, GRR = GRR, T = T, lambda = lamda, Dt = Dt))
    colnames(parameters) <- "Parameters"
    print(parameters)
    if(poli>1)
    {
       matrizOut=data.frame(ltb[,1:(length(estadios)-2)],poli*ltb[,(ncol(ltb)-2):ncol(ltb)])
       nombres=colnames(ltb);nombres2=c(paste(nombres[1:(length(estadios)-2)],"Host",sep="_"),paste(nombres[(ncol(ltb)-2):ncol(ltb)],"Paras",sep="_"))
       colnames(matrizOut)=nombres2
    }else{matrizOut=ltb}
    return(list(ltb=matrizOut,m.x=m.x,s.x=s.x,l.x=as.matrix(l.x),parametro=parameters))
 }else
 {
    parameters <- t(data.frame(r = NA, Ro = NA, GRR = NA, T = NA, lambda = NA, Dt = NA))
    return(list(parametro=parameters))
 }
}

#############################################################################################################
#########1####################################################################################################
grafflife1 <- function (data, estad,parametros,ltb,corrx,labx,laby,lgx,lgy,tam,titulo = NULL, directory = NULL){
    estadios <- estad
    lifetable <- ltb
    Day <- 0:(nrow(lifetable) - 1)
    rf <- parametros[1, 1]
    s.x <- m.x <- l.x <- rep(0, nrow(lifetable))
    for (i in 1:nrow(lifetable)) l.x[i] <- sum(lifetable[i, 1:(length(estadios))])/ncol(data)
    Cx <- rep(0, nrow(lifetable))
    for (i in 1:nrow(lifetable)) Cx[i] <- l.x[i] * exp(-rf *
        (Day[i] + 1))/sum(l.x * exp(-rf * (Day + 1)))
    par(cex=tam)
    plot(Day, lifetable[,1], type = "l", col = 1, frame = F, xlim=corrx, axes=F,xaxt = "n",xlab =labx , ylab =laby,main=titulo)
    axis(1, xaxp=c(corrx,5))
    axis(2,las=2)
    for (j in 1:length(estadios)) lines(Day, lifetable[1:(length(estadios))][, j],
        type = "l", col = j)
    lines(Day, (lifetable[,length(estadios)-1] + lifetable[,length(estadios)-2]), type = "l", col = j + 1)
    estadios2 <- c(estadios, "Adults")
   legend(lgx,lgy, estadios2, cex = 0.8, col = 1:(j + 1), lty = 1)
 }
 #############################################################################################################
grafflife2 <- function (data, estad,parametros,ltb,corrx,corry,labx,laby,lgx,lgy,tam,titulo = NULL, directory = NULL){
    estadios <- estad
    lifetable <- ltb
    Day <- 0:(nrow(lifetable) - 1)
    rf <- parametros[1, 1]
    s.x <- m.x <- l.x <- rep(0, nrow(lifetable))
    for (i in 1:nrow(lifetable)) l.x[i] <- sum(lifetable[i, 1:(length(estadios))])/ncol(data)
    Cx <- rep(0, nrow(lifetable))
    for (i in 1:nrow(lifetable)) Cx[i] <- l.x[i] * exp(-rf *
        (Day[i] + 1))/sum(l.x * exp(-rf * (Day + 1)))
    par(cex=tam)
    plot(Day, Cx, type = "l", xlim =corrx, ylim =corry , frame = F, xlab =labx , ylab =laby,main=titulo)
    lineas <- matrix(rep(0, length(estadios) * nrow(lifetable)),
        nrow = nrow(lifetable), ncol = length(estadios))
    for (j in 1:length(estadios)) for (i in 1:nrow(lifetable)) {
        lineas[i, j] <- Cx[i] * (lifetable[1:(length(estadios))][i, j]/
        sum(lifetable[1:(length(estadios))][i, ]))
    }
    for (j in 1:length(estadios)) lines(Day, lineas[, j], type = "l", col = j)
    legend(lgx,lgy, estadios, cex = 0.8, col = 1:(j + 1), lty = 1)
}
 #############################################################################################################
grafflife3 <- function (data, estad,parametros,ltb,corrx,labx,laby,tam,titulo = NULL, directory = NULL){
    estadios <- estad
    lifetable <- ltb
    Day <- 0:(nrow(lifetable) - 1)
    rf <- parametros[1, 1]
    s.x <- m.x <- l.x <- rep(0, nrow(lifetable))
    for (i in 1:nrow(lifetable)) l.x[i] <- sum(lifetable[i, 1:(length(estadios))])/ncol(data)
    par(cex=tam)
    plot(Day, l.x, type = "s", xlim = corrx,ylim = c(0, 1), frame = F, col = 4, xlab =labx,
        ylab = laby,main=titulo)
}
#################################################
##################################################
simultemp<-function(N,sexratio,isFixed,Rs,temp,resps,xi,steps,poli=1, ..){
	n1=length(temp)
	pars=data.frame(Temper=NA,r=NA,Ro=NA,GRR=NA,T=NA,lambda=NA,Dt=NA)

	for(i in 1:n1)
	{
		Table=data.frame(V1=rep(temp[i],365),V2=rep(temp[i],365))

		if(isFixed){
			Rs=rep(sexratio,365)
		}else{
			Table2=cbind(id=1:nrow(Table),Table)
			paramR=params$paramR
			ff=params$ff
			Rs=apply(Table2,1,RateR,Table2,nmax=nrow(Table),steps=steps,ff,paramR)
		}


		for(j in 1:reps[i])
		{
			#cuadro<-Ratios(Table,modelim,modelme,estadios, xi,steps)$MATRIZ

			cuadro <- matrizA(estadios, hfeno, Table, steps)
			#insectos<-simulacion(cuadro,estadios,N,M,hfeno,sexratio)
			insectos<-simulacion(cuadro,estadios,N,M,hfeno,Rs)$mat
			insects<-insectos
			lftable<-life.table (insects, estadios)
			ltb<-lftable$life.table
			#summ_stat<-statist(insects, estadios,poli)			
			#obs_stat<-summ_stat$Observations
			#survi_stat<-summ_stat$Survival
			#mort_stat<-summ_stat$Mortality
			pravida<-parameters(N, estadios, ltb,poli)
			parametros<-pravida$parametro
			pars=rbind(pars,c(Temper=temp[i],t(parametros)))
		}
	}
	pars=pars[-1,]
	return(pars)
}

#################################################
##################################################
graphsimulp<-function(pars,modls=c(1,1,1,1,1,1),direc,ecua,nomec,ejex,ejey,tecu,grises,tit,ax,ay){
	
	if(grises==TRUE){ccol=c("gray20","gray50")}else{ccol=c("royalblue","red")}        
	
	
	
	etip=c("Rm","Ro","GRR","GL","Lambda","Dt")
	
	nfile=etip[!is.na(modls)]
	
	NMS=paste(nfile,collapse="_")
	
	pathIm = paste(direc,NMS,".jpg",sep="")
	
	jpeg(pathIm, width = 720, height = 800,quality = 100)##
	
	
	
	
	
	mx<-mmx<-max(pars[,1],na.rm = TRUE);mnx<-min(pars[,1],na.rm = TRUE)
	
	
	
	if(mx<=40){mx=40}else{mx=mx+10}
	
	N=length(modls[!is.na(modls)]) 
	
	if(N==1){par(mfrow=c(1,1));  if(is.na(ejex)){ejex=-1.73};  if(is.na(ejey)){ejey=-1.61}  }
	
	if(N==2){par(mfrow=c(2,1));  if(is.na(ejex)){ejex=-0.69};  if(is.na(ejey)){ejey=-1.60}  }
	
	if(N==3 || N==4){par(mfrow=c(2,2));  if(is.na(ejex)){ejex=-0.86};  if(is.na(ejey)){ejey=-0.84}  }
	
	if(N==5 || N==6){par(mfrow=c(3,2));  if(is.na(ejex)){ejex=-0.66};  if(is.na(ejey)){ejey=-1.25}  }
	
	
	
	
	
	vars=c(2,3,4,5,6,7)
	
	nomv=c(Temper="Temperature",r = expression(Intrinsic~rate~group("(",r[m],")")),Ro = "Net reproduction rate (Ro)",GRR = "Gross reproduction rate (GRR)",T = "Generation length in days (GL)",lambda = expression(Finite~rate~of~increase~group("(",symbol("l"),")") ),Dt = "Doubling time (Dt)") ## nombre de los parametros
	
	modelss=c("cubic","quadratic","logarithmic","exponential")  ## tipo de modelos propuestos
	
	
	
	modsel<-matrix(NA,6,4)
	
	colnames(modsel)<-c("R2","R2_Adj","AIC","Deviance")
	
	sal2<-et<-list(1);ubix=20
	
	
	
	for(k in vars) ## gráfico por cada parámetro
	
	{ 
		
		if(ecua==TRUE){if(is.na(modls[k-1])){next}} ## este codigo hace que solo ejecute los parametros que tienen un modelo escogido
		
		my=max(pars[,k],na.rm = TRUE);mny=min(pars[,k],na.rm = TRUE)
		
		if(mny<(my-mny) & my>0){mny=0}
		
		plot(pars[,1],pars[,k],xlim=c(0,mx),ylim=c(mny-0.05*(abs(mny)),my+0.13*(abs(my))),xlab="Temperature °C",ylab=nomv[k],pch=19,col=ccol[1],axes=FALSE)  ## 1
		
		axis(1,at=seq(0,mx,by=10),line=ejex,cex.axis=ax) ## falta incrementar los tamaños de los subtitulos a 12
		
		axis(2,las=2,line=ejey,cex.axis=ay)
		
		
		y=pars[,k]
		
		x=pars[,1];x=x[is.finite(y)];y=y[is.finite(y)]
		
		
		
		if(ecua==TRUE)
		
		{
			
			## Aqui se definen las funciones de los modelos
			
			if(modls[k-1]==1){f1=as.formula(y ~ b0 + b1*x + b2*x^2 + b3*x^3);inis=list(b0=1,b1=1,b2=2,b3=2);ff <- function(x){b0 + b1*x + b2*x^2 + b3*x^3}}
			
			if(modls[k-1]==2){f1=as.formula(y ~ b0 + b1*x + b2*x^2);inis=list(b0=1,b1=1,b2=2);ff <- function(x){b0 + b1*x + b2*x^2}}
			
			if(modls[k-1]==3){f1=as.formula(y ~ b0 + b1*x + b2*log(x));inis0=coef(lm(y~x+I(log(x))));names(inis0)=NULL;inis=list(b0=inis0[1],b1=inis0[2],b2=inis0[3]);ff <- function(x){b0 + b1*x + b2*log(x)}}
			
			if(modls[k-1]==4){f1=as.formula(y ~ b0 + b1*x + b2*exp(x));inis0=coef(lm(y~x+I(exp(x))));names(inis0)=NULL;inis=list(b0=inis0[1],b1=inis0[2],b2=inis0[3]);ff <- function(x){b0 + b1*x + b2*exp(x)}}
			
			
			
			if(tit==TRUE){title(paste("Using the",modelss[modls[k-1]],"model"))}
			
			
			
			## Se estiman los parametros usando el algoritmo marc quart
			
			out <- nls(f1, start = inis,trace = TRUE)
			
			yl=fitted(out)
			
			sqe=sum(residuals(out)^2)
			
			sal <- coef(out);sal2[[k-1]]=sal
			
			
			
			r<-1-sqe/sum((y-mean(y))^2) # este es el R^2
			
			r_ajus<- 1 - ((length(x) - 1) / (length(x) - length(sal))) * (1-r)  # este es el R^2 adjus
			
			AC<-AIC(out) # este es el AIC 
			
			Dev=deviance(out) # este es el indicador  Deviance
			
			
			
			modsel0<-c(R2=round(r,3),R2_Adj=round(r_ajus,3),AIC=round(AC,3),Deviance=round(Dev,3))
			
			modsel[k-1,]<-modsel0
			
			
			
			## En esta parte solo se agrega la acuacion al gráfico
			
			if((k-1)==1)
			
			{
				
				if(modls[k-1]==1){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(bolditalic(r[m]) == list(sal1) + list(sal2)*T + list(sal3)*T^2 + list(sal4)*T^3,list(sal1=sal[1],sal2=sal[2],sal3=sal[3],sal4=sal[4])),cex=tecu))}
				
				if(modls[k-1]==2){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(bolditalic(r[m]) == list(sal1) + list(sal2)*T + list(sal3)*T^2,list(sal1=sal[1],sal2=sal[2],sal3=sal[3])),cex=tecu))}
				
				if(modls[k-1]==3){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(bolditalic(r[m]) == list(sal1) + list(sal2)*T + list(sal3)*log(T),list(sal1=sal[1],sal2=sal[2],sal3=sal[3])),cex=tecu))}
				
				if(modls[k-1]==4){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(bolditalic(r[m]) == list(sal1) + list(sal2)*T + list(sal3)*e^(T),list(sal1=sal[1],sal2=sal[2],sal3=sal[3])),cex=tecu))}
				
			}
			
			
			
			if((k-1)==2 || (k-1)==3 || (k-1)==4 || (k-1)==6)
			
			{
				
				if(modls[k-1]==1){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(bolditalic(list(etip1)) == list(sal1) + list(sal2)*T + list(sal3)*T^2 + list(sal4)*T^3,list(etip1=etip[k-1],sal1=sal[1],sal2=sal[2],sal3=sal[3],sal4=sal[4])),cex=tecu))}
				
				if(modls[k-1]==2){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(bolditalic(list(etip1)) == list(sal1) + list(sal2)*T + list(sal3)*T^2,list(etip1=etip[k-1],sal1=sal[1],sal2=sal[2],sal3=sal[3])),cex=tecu))}
				
				if(modls[k-1]==3){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(bolditalic(list(etip1)) == list(sal1) + list(sal2)*T + list(sal3)*log(T),list(etip1=etip[k-1],sal1=sal[1],sal2=sal[2],sal3=sal[3])),cex=tecu))}
				
				if(modls[k-1]==4){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(bolditalic(list(etip1)) == list(sal1) + list(sal2)*T + list(sal3)*e^(T),list(etip1=etip[k-1],sal1=sal[1],sal2=sal[2],sal3=sal[3])),cex=tecu))}
				
			} 
			
			
			
			if((k-1)==5)
			
			{
				
				if(modls[k-1]==1){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(symbol(list(etip1)) == list(sal1) + list(sal2)*T + list(sal3)*T^2 + list(sal4)*T^3,list(etip1=etip[k-1],sal1=sal[1],sal2=sal[2],sal3=sal[3],sal4=sal[4])),cex=tecu))}
				
				if(modls[k-1]==2){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(symbol(list(etip1)) == list(sal1) + list(sal2)*T + list(sal3)*T^2,list(etip1=etip[k-1],sal1=sal[1],sal2=sal[2],sal3=sal[3])),cex=tecu))}
				
				if(modls[k-1]==3){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(symbol(list(etip1)) == list(sal1) + list(sal2)*T + list(sal3)*log(T),list(etip1=etip[k-1],sal1=sal[1],sal2=sal[2],sal3=sal[3])),cex=tecu))}
				
				if(modls[k-1]==4){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(symbol(list(etip1)) == list(sal1) + list(sal2)*T + list(sal3)*e^(T),list(etip1=etip[k-1],sal1=sal[1],sal2=sal[2],sal3=sal[3])),cex=tecu))}
				
			}
			
			if(nomec==TRUE){eval(et[[k-1]])}
			
			for (i in names(sal))
			
			{ 
				
				temp <- sal[i] 
				
				storage.mode(temp) <- "double"
				
				assign(i, temp)
				
			}
			
			curve(ff,add=TRUE,from=mnx,to=mmx,col=ccol[2],lwd=2)  ## 2
			
		}
		
	}
	
	dev.off()
	
	
	
	if(ecua==TRUE) 
	
	{
		
		jpeg(paste(direc,"Indicadores",NMS,".jpg",sep=""), width = 770, height = 500,quality = 100)##
		
		
		
		plot(1:12,1:12,type="n",axes=FALSE,xlab="",ylab="")
		
		grid(0,6,lty=1)
		
		ubix=8.5;ubix2=8.5
		
		lines(c(4.5,4.5),c(0,13),lty=1,lwd=2,col="gray90",type = "l")
		
		box(lty = "solid", col = "gray30")
		
		
		
		if(!is.na(modls[1])){
			
			sal=sal2[[1]]
			
			my=(12/1.09);k=2;if(modls[1]==1){ubix=8.5};if(modls[1]==2){ubix=7.6};if(modls[1]==3 || modls[1]==4){ubix=6.7}
			
			text(1.8,12,nomv[2], col="gray40") ## para rm
			
			eval(et[[1]])
			
			text(ubix2-3, my*(1.09)-0.94,  substitute(bolditalic(R^{2})==list(r2),list(r2=modsel[1,1])))
			
			text(ubix2-1.2, my*(1.09)-0.97,  substitute(bolditalic(R^{2}~Adj)==list(r2a),list(r2a=modsel[1,2])))
			
			text(ubix2+0.6, my*(1.09)-1,  substitute(bolditalic(AIC)==list(aic),list(aic=modsel[1,3])))
			
			text(ubix2+2.6, my*(1.09)-1,  substitute(bolditalic(Deviance)==list(dev),list(dev=modsel[1,4])))
			
		}
		
		
		
		if(!is.na(modls[2])){
			
			sal=sal2[[2]]
			
			my=(12/1.09)-1.8;k=3;if(modls[2]==1){ubix=8.5};if(modls[2]==2){ubix=7.6};if(modls[2]==3 || modls[2]==4){ubix=6.7}
			
			text(2.2,10,nomv[3], col="gray40") ## para Ro
			
			eval(et[[2]])
			
			text(ubix2-3, my*(1.09)-0.94,  substitute(bolditalic(R^{2})==list(r2),list(r2=modsel[2,1])))
			
			text(ubix2-1.2, my*(1.09)-0.97,  substitute(bolditalic(R^{2}~Adj)==list(r2a),list(r2a=modsel[2,2])))
			
			text(ubix2+0.6, my*(1.09)-1,  substitute(bolditalic(AIC)==list(aic),list(aic=modsel[2,3])))
			
			text(ubix2+2.6, my*(1.09)-1,  substitute(bolditalic(Deviance)==list(dev),list(dev=modsel[2,4])))
			
		}
		
		
		
		if(!is.na(modls[3])){
			
			sal=sal2[[3]]
			
			my=(12/1.09)-3.6;k=4;if(modls[3]==1){ubix=8.5};if(modls[3]==2){ubix=7.6};if(modls[3]==3 || modls[3]==4){ubix=6.7}
			
			text(2.5,8,nomv[4], col="gray40") ## para GRR
			
			eval(et[[3]])
			
			text(ubix2-3, my*(1.09)-0.94,  substitute(bolditalic(R^{2})==list(r2),list(r2=modsel[3,1])))
			
			text(ubix2-1.2, my*(1.09)-0.97,  substitute(bolditalic(R^{2}~Adj)==list(r2a),list(r2a=modsel[3,2])))
			
			text(ubix2+0.6, my*(1.09)-1,  substitute(bolditalic(AIC)==list(aic),list(aic=modsel[3,3])))
			
			text(ubix2+2.6, my*(1.09)-1,  substitute(bolditalic(Deviance)==list(dev),list(dev=modsel[3,4])))
			
		}
		
		
		
		if(!is.na(modls[4])){
			
			sal=sal2[[4]]
			
			my=(12/1.09)-5.4;k=5;if(modls[4]==1){ubix=8.5};if(modls[4]==2){ubix=7.6};if(modls[4]==3 || modls[4]==4){ubix=6.7}
			
			text(2.5,6,nomv[5], col="gray40") ## para GL
			
			eval(et[[4]])
			
			text(ubix2-3, my*(1.09)-0.94,  substitute(bolditalic(R^{2})==list(r2),list(r2=modsel[4,1])))
			
			text(ubix2-1.2, my*(1.09)-0.97,  substitute(bolditalic(R^{2}~Adj)==list(r2a),list(r2a=modsel[4,2])))
			
			text(ubix2+0.6, my*(1.09)-1,  substitute(bolditalic(AIC)==list(aic),list(aic=modsel[4,3])))
			
			text(ubix2+2.6, my*(1.09)-1,  substitute(bolditalic(Deviance)==list(dev),list(dev=modsel[4,4])))
			
		}
		
		
		
		if(!is.na(modls[5])){
			
			sal=sal2[[5]]
			
			my=(12/1.09)-7.2;k=6;if(modls[5]==1){ubix=8.5};if(modls[5]==2){ubix=7.6};if(modls[5]==3 || modls[5]==4){ubix=6.7}
			
			text(2.4,4,nomv[6], col="gray40") ## para lambda
			
			eval(et[[5]])
			
			text(ubix2-3, my*(1.09)-0.94,  substitute(bolditalic(R^{2})==list(r2),list(r2=modsel[5,1])))
			
			text(ubix2-1.2, my*(1.09)-0.97,  substitute(bolditalic(R^{2}~Adj)==list(r2a),list(r2a=modsel[5,2])))
			
			text(ubix2+0.6, my*(1.09)-1,  substitute(bolditalic(AIC)==list(aic),list(aic=modsel[5,3])))
			
			text(ubix2+2.6, my*(1.09)-1,  substitute(bolditalic(Deviance)==list(dev),list(dev=modsel[5,4])))
			
		}
		
		
		
		if(!is.na(modls[6])){
			
			sal=sal2[[6]]
			
			my=(12/1.09)-9;k=7;if(modls[6]==1){ubix=8.5};if(modls[6]==2){ubix=7.6};if(modls[6]==3 || modls[6]==4){ubix=6.7}
			
			text(2,2,nomv[7], col="gray40") ## para Dt
			
			eval(et[[6]])
			
			text(ubix2-3, my*(1.09)-0.94,  substitute(bolditalic(R^{2})==list(r2),list(r2=modsel[6,1])))
			
			text(ubix2-1.2, my*(1.09)-0.97,  substitute(bolditalic(R^{2}~Adj)==list(r2a),list(r2a=modsel[6,2])))
			
			text(ubix2+0.6, my*(1.09)-1,  substitute(bolditalic(AIC)==list(aic),list(aic=modsel[6,3])))
			
			text(ubix2+2.6, my*(1.09)-1,  substitute(bolditalic(Deviance)==list(dev),list(dev=modsel[6,4])))
			
		}
		
		
		
		dev.off()
		
	}
	
	return(list(coefs=sal2,modsel=modsel,et=et, pathIm=pathIm))
	
	
 }

 #############################################################################################################
##########################################################################################################
distrimodel<-function(vec,sll){  # sl es 1, 2 ,3,4 ...)

  # probit
  if(hfeno$distri_dv[[sll]]=="probit")
  {
    WW<-pnorm(log(vec)*hfeno$slope_dv[[sll]])
  }
  # logit
  if(hfeno$distri_dv[[sll]]=="logit")
  {
    WW<-1/(1+exp(-(log(vec)*hfeno$slope_dv[[sll]])))
  }
  # cloglog
  if(hfeno$distri_dv[[sll]]=="cloglog")
  {
    WW<-1-exp(-exp(log(vec)*hfeno$slope_dv[[sll]]))
  }
return(WW)
}
#############################################################################################################
#############################################################################################################
Estado <- function (data, estad,est){
   maduros   <-  estad[(length(estad)-1):(length(estad))]
   if(est==maduros[1]){
        datat <- as.matrix(data)
       num <- matrix(as.numeric(datat), nrow(datat), ncol(datat))
        a <- as.list(as.data.frame(t(num)))
        b <- as.list(rep(0, nrow(data)))
        for (i in 1:nrow(data)) b[[i]] <- sum(as.matrix(table(a[[i]])))
        female <- as.numeric(as.character(as.matrix(b)))
        return(female)

   }
   if(est!=maduros[1]){
        datat <- as.matrix(data)
       estadoa <- as.list(as.data.frame(t(datat)))
       estadob <- as.list(rep(0, nrow(data)))
       for (i in 1:nrow(data)) {
           estadoa[[i]] <- as.matrix(estadoa[[i]])
           estadob[[i]] <- length(subset(estadoa[[i]], estadoa[[i]] == est))
        }
        estados <- as.numeric(as.character(as.matrix(estadob)))
        return(estados)
   }
}
##########################################
distriModelAdults<-function(vec,k, sizeInmaduros){  # sl es 4 o 5
    if(k == (sizeInmaduros +1)){
        # probit
        if(hfeno$distri_snh=="probit"){
            WW<-pnorm(log(vec)*hfeno$slope_snh)
        }
        # logit
        if(hfeno$distri_snh=="logit"){
            WW<-1/(1+exp(-(log(vec)*hfeno$slope_snh)))
        }
        # cloglog
        if(hfeno$distri_snh=="cloglog"){
            WW<-1-exp(-exp(log(vec)*hfeno$slope_snh))
        }
    }
    
    if(k == (sizeInmaduros+2)){
        # probit
        if(hfeno$distri_snm=="probit"){
            WW<-pnorm(log(vec)*hfeno$slope_snm)
        }
        # logit
        if(hfeno$distri_snm=="logit"){
            WW<-1/(1+exp(-(log(vec)*hfeno$slope_snm)))
        }
        # cloglog
        if(hfeno$distri_snm=="cloglog"){
            WW<-1-exp(-exp(log(vec)*hfeno$slope_snm))
        }
    }
    return(WW)
}
#########################################################
                                                                     
                                                                     
                                                                     
                                             
simulacion<-function(cuadro,estadios,n,m,hfeno,Rs){
	mat<-matrix("dead",ncol=n,nrow=(m))
	for(z in 1:n){
		k <- 1
		Insect_age <- 0
		survival_p <- 1
		dead <- 0
		Day <- 1
		slope<- 0
		djx <- 0
		mjx <- 0
		kj <- NA
		rand_d <- round(runif(1),3) #' random number between 0 and 1 with 3 digits
		rand_m <- round(runif(1),3)
		rand_ovi <- round(runif(1),3)
		rand_sex = 0
		
		if(rand_ovi == 0){
			rand_ovi <- 0.0001
		}
		
		sizeInmaduros = length(estadios[-(length(estadios)-1):-(length(estadios))])
		numFemale = sizeInmaduros + 1
		numMale = sizeInmaduros + 2
		
		#loop for immature life stages
		while(dead != 1 || k <= sizeInmaduros){
			if(k > sizeInmaduros){
				break
			}
			if(dead ==1){
				break
			}
			
			#read the parameter for variation and for each specific stage
			slope = hfeno$slope_dv[[k]]
			
			#read value from Matrix "A"
			djx <- cuadro[Day,k*2-1]
			mjx <- cuadro[Day,k*2]
			
			#calculate physiological age
			if(Insect_age == 0){
				Insect_age = Insect_age + (djx / 2)
			}else{
				Insect_age <- Insect_age + djx
			}
			
			dev_p <- distrimodel(Insect_age,k)
			
			################################################
			# Se evalua si el insecto muere en cierto estado
			################################################
			
			#evaluate survival
			survival_p = survival_p * ((1 - mjx) ^ djx)  ## esta probabilidad es muy alta
			survival_rand = (1 - rand_m)
			
			if(survival_p < survival_rand || Day==365){
				dead = 1
				mat[Day,z] <- "dead"
				survival_p = 1
				k = 1
			}
			
			##################################################
			# Se evalua si el insecto pasa al siguiente estado
			##################################################
			
			#evaluate development to next stage
			if(dev_p > rand_d){   ### si dev_p es demasiado pequeño nunca va pasar al siguiente estado fenológico
				#read djx from Matrix "A" of the next stage
				djx <- cuadro[Day,(k+1)*2-1]
				Insect_age = djx / 2
				
				rand_d <- round(runif(1),3) #' random number for the next stage
				rand_m <- round(runif(1),3)
				
				survival_p = 1
				k = k + 1
			}
			
			#inmature dev stages
			if(k <= sizeInmaduros){
				kj = estadios[k]
			}
			
			if(k > sizeInmaduros){
				Insect_age = 0
			}
			if(dead == 1){
				kj = "dead"
			}
			mat[Day,z] <- kj
			Day = Day + 1
			dead
			
		}
		# end loop for inmatures
		
		#print(paste(Day,Rs[Day],sep="_")) ## Taza calculada en base ala temperatura que tiene según el dia cuando el insecto termina de ser inmaduro
		if(k == (sizeInmaduros+1) && Insect_age == 0 && Day!=nrow(cuadro) + 1){
			rand_sex <- round(runif(1),3)
			
			if(rand_sex < Rs[Day]){ ## Aquí se aplica
				k = numFemale
			}else{
				k = numMale
			}
		}
		
		#loop for adult males
		while(dead != 1 && k != numFemale && Day!=nrow(cuadro) + 1){
			if(dead ==1){
				break
			}
			
			#read the parameter for variation and for each specific stage
			slope = hfeno$slope_snm
			
			#read value from Matrix "A"
			djx <- cuadro[Day,k*2-2]
			
			#calculate physiological age
			if(Insect_age == 0){
				Insect_age = Insect_age + (djx / 2)
			}else{
				Insect_age <- Insect_age + djx
			}
			
			dev_p <- distriModelAdults(Insect_age,k, sizeInmaduros)
			#dev_p = log(Insect_age)
			#dev_p <- distrimodel(dev_p,k)
			
			#evaluate survival of males
			if(dev_p > rand_d){
				dead <- 1
				Insect_age = 0
				
				rand_d <- round(runif(1),3) #' random number for the next stage
				rand_m <- round(runif(1),3)
				
				survival_p = 1
			}
			
			kj = estadios[k]
			
			if(dead == 1){
				kj = "dead"
			}
			
			mat[Day,z] <- kj
			Day = Day + 1
		}# end loop for males
		#print(Day)
		#loop for adult females
		while(dead != 1 && k != numMale && Day!= nrow(cuadro) + 1){
			if(dead ==1){
				break
			}
			parametrosc <- hfeno$povih_h
			parametrosc<-as.list(parametrosc)
			for (i in names(parametrosc)){
				temp <- parametrosc[[i]]
				storage.mode(temp) <- "double"
				assign(i, temp)
			}
			formulac <- hfeno$fovih_h
			forexc <- formulac[[length(formulac)]]
			funcionc <- as.expression(forexc)
			
			slope = hfeno$slope_snh
			djx <- cuadro[Day,k*2-1]
			
			#calculate physiological age
			if(Insect_age == 0){
				djx <- djx + cuadro[Day-1,k*2-1]
			}
			
			Insect_age = Insect_age + djx
			dev_p <- distriModelAdults(Insect_age,k, sizeInmaduros)
			#dev_p = log(Insect_age)
			#dev_p <- distrimodel(dev_p,k)
			
			#evaluate survival of females
			if(dev_p > rand_d){
				dead = 1
				Insect_age = 0
				rand_d <- round(runif(1),3) #' random number for the next stage
				rand_m <- round(runif(1),3)
				
				survival_p = 1
			}
			
			#calculate oviposition
			#read parameters
			if(dead < 1){
				alfa = hfeno$povih_h[[1]]
				beta = hfeno$povih_h[[2]]
				cv = 0.3 # Averiguar acerca de este valor
				
				#read ovitot from Matrix "A"
				Ovitot = cuadro[Day-1,k*2+1]
				
				#calculate age-dependent oviposition frequency
				x  <-  Insect_age - djx
				OviFrec1  <- eval(funcionc)
				x  <-  Insect_age
				OviFrec2  <- eval(funcionc)
				OviFrec = OviFrec2 - OviFrec1
				
				#calculate oviposition
				Ovitot = qnorm(rand_ovi, Ovitot, Ovitot * cv);Ovitot[Ovitot<0]=0
				Oviposition = Ovitot * OviFrec
				Oviposition = round(Oviposition, 0)
				kj<- Oviposition
			}
			
			if(dead == 1){
				kj = "dead"
			}
			#print("error5");print(paste(Day,z))
			
			mat[Day,z] <- kj
			Day = Day + 1
		}#end loop females and oviposition
		
		Day = 1
	}# end loop individuos
	return(list(mat=data.frame(mat)))
}
##########################################
# Funcion generadora de Tasas y mortalidad
##########################################
RateI<-function(vec,Table2,Ki,parametrosc,parametrosm=NULL,funciont,funcionm=NULL,nmax,steps,J=NA)
{
	for (i in names(parametrosc)){temp <- parametrosc[i];storage.mode(temp) <- "double";assign(i, temp)}
	i=0:(steps-1)
	T1=((vec[3]-vec[2])/2)*cos(pi*(i+0.5)/steps) + (vec[3]+vec[2])/2
	if(vec[1]!=nmax){T2<-((vec[3]-Table2[vec[1]+1,2])/2)*cos(pi*(i+0.5)/steps) + (vec[3]+Table2[vec[1]+1,2])/2}else{
		T2<-((vec[3]-Table2[1,2])/2)*cos(pi*(i+0.5)/steps) + (vec[3]+Table2[1,2])/2}
	if((modelim[Ki]==1 || modelim[Ki]==2 || modelim[Ki]==3 ||
				modelim[Ki]==4 || modelim[Ki]==5 || modelim[Ki]==6 ||
				modelim[Ki]==7 || modelim[Ki]==8 || modelim[Ki]==9 ||
				modelim[Ki]==10 || modelim[Ki]==11 || modelim[Ki]==12 || modelim[Ki]==13 || modelim[Ki]==14) & is.na(J)){x = T1 + 273.15;x2 = T2 + 273.15}else{x = T1;x2 = T2}
	
	
	Rat=eval(funciont);if(!is.na(J)){Rat[Rat<=0]=0}
	Ratetot1<-(sum(Rat))/steps    ### aqui evalua la tasa de desarrollo de todas las divisiones
	x=x2;Rat=eval(funciont);if(!is.na(J)){Rat[Rat<=0]=0}
	Ratetot2<-(sum(Rat))/steps    ### aqui evalua la tasa de desarrollo de todas las divisiones del segundo dia
	Rate=(Ratetot1+Ratetot2)/2
	
	if(!is.null(funcionm))
	{
		for (i in names(parametrosm)){temp <- parametrosm[i];storage.mode(temp) <- "double";assign(i, temp)}
		x=T1;M1=eval(funcionm);M1[M1>1]=1; Mortality1 <- (sum(M1))/steps    ### aqui evalua la Mortalidad de todas las divisiones
		x=T2;M2=eval(funcionm);M2[M2>1]=1; Mortality2 <- (sum(M2))/steps    ### aqui evalua la Mortalidad de todas las divisiones del segundo dia
		Mortality=(Mortality1+Mortality2)/2
		return(c(Rate1=Ratetot1,Rate2=Rate,Mortality1=Mortality1,Mortality2=Mortality))
	}else{return(c(Rate1=Ratetot1,Rate2=Rate))}
}
##############################################


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
		funcionm <- as.expression(hfeno$mortal[[K]][[3]])
		
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
#####################################
AgeClases<- function(cuadro, estadios, hfeno){
	matrizA=cuadro
	sizeMatrizA = nrow(matrizA)
	inmad <-  estadios[-(length(estadios)-1):-(length(estadios))]
	numInmaduros = length(inmad)
	maduros   <-  estadios[(length(estadios)-1):(length(estadios))]
	numMaduros = length(maduros)
	estadiosAges<-list(1)
	ageClassMatriz = matrix(NA, nrow=sizeMatrizA, ncol=(numInmaduros+numMaduros))                                        
	
	Day = 1
	Ageclass = 1                     
	age = 0
	
	#primera linea  
	for(k in 1:(numInmaduros + numMaduros)){
		Ageclass = 1
		age = 0  
		Ageclass0 = c(rep(NA, 253))# porq 253??
		matrizAges = matrix(NA, nrow=sizeMatrizA, ncol=(length(Ageclass0)+1))
		#matrizAges = matrix(NA, nrow=366, ncol=(length(Ageclass0)+1))
		matrizAges[Day,1] = 0
		
		## if(k==1){
		##     Ageclass = Ageclass + 1
		## }
		
		Ageclass = Ageclass + 1
		#for(Ii in 2:length(Ageclass0)){
		for(Ii in 1:length(Ageclass0)){
			if(k == (numInmaduros + numMaduros)){
				Ageclass0[Ii] = matrizA[(sizeMatrizA+1)-Ii,k*2-2]              
			}else{
				Ageclass0[Ii] = matrizA[(sizeMatrizA+1)-Ii,k*2-1]
			}
			
			age = age + Ageclass0[Ii]
			matrizAges[Day,Ii+1] = age
		}
		if(k <= numInmaduros){
			ageClassMatriz[Day,k] = length(matrizAges[1,][matrizAges[1,]<2])+1
		}else{
			ageClassMatriz[Day,k] = length(matrizAges[1,][matrizAges[1,]<1.8])+1
		}
		estadiosAges[[k]] = matrizAges
	}
	
	# a partir de la 2da linea
	for(k in 1:(numInmaduros + numMaduros)){
		Day=2
		Ageclass = 1
		age = 0
		matrizAges = estadiosAges[[k]]
		while(Day <= sizeMatrizA){#porq 367 en excel?
			matrizAges[Day,1] = 0
			if(k == (numInmaduros + numMaduros)){        
				Rate = matrizA[Day,k*2-2]
			}else{
				Rate = matrizA[Day,k*2-1]            
			}
			
			if(k <= numInmaduros){
				while(age < 2){ # cambiar a valor del slop
					ageacc = matrizAges[Day-1, Ageclass]
					age = ageacc + Rate
					#matrizAges[Day-1, Ageclass+1] = age
					matrizAges[Day, Ageclass+1] = age
					Ageclass = Ageclass + 1
					
					if(Ageclass == 254){
						age = 3
					}
				}
				ageClassMatriz[Day,k] = Ageclass
			}else{
				while(age < 1.8){ # cambiar a valor del slop
					ageacc = matrizAges[Day-1, Ageclass]
					age = ageacc + Rate
					#matrizAges[Day-1, Ageclass+1] = age
					matrizAges[Day, Ageclass+1] = age
					Ageclass = Ageclass + 1
					
					
					if(Ageclass == 255){
						age = 3
					}
				}
				ageClassMatriz[Day,k] = Ageclass
			}
			age = 0
			Day = Day + 1
			Ageclass = 1
		}
		estadiosAges[[k]] = matrizAges
	}
	
	## Ovifreq ##
	matrizAges = estadiosAges[[numInmaduros+1]]
	oviFreq = matrix(NA, nrow=sizeMatrizA, ncol=length(Ageclass0))
	Day = 1
	Ageclass = 1
	age1 = 0
	age2 = 0
	
	parametrosc <- hfeno$povih_h
	parametrosc<-as.list(parametrosc)
	for (i in names(parametrosc)){
		temp <- parametrosc[[i]]
		storage.mode(temp) <- "double"
		assign(i, temp)
	}
	
	formulac <- hfeno$fovih_h
	forexc <- formulac[[length(formulac)]]
	funcionc <- as.expression(forexc)
	
	while(Day < sizeMatrizA){#porq 367 en excel?
		matrizAges[Day,1] = 0
		age1 = matrizAges[Day, Ageclass]
		x = age1
		OviFrec1 = eval(funcionc)
		
		while(age2 < 1.8){
			age2 = matrizAges[Day, Ageclass]
			x= age2
			OviFrec2 = eval(funcionc)
			OviFrec = OviFrec2 - OviFrec1
			OviFrec1 = OviFrec2
			if(Ageclass == 1){
				OviFrec = 0
			}
			oviFreq[Day,Ageclass] = OviFrec
			Ageclass = Ageclass + 1
		}
		age2 = 0
		Day = Day + 1
		Ageclass = 1
	}
	
	
	return(list(estadiosAges=estadiosAges, oviFreq=oviFreq, ageClassMatriz=ageClassMatriz))
}

#################################
# Tasas variables ###############

RateR<-function(vec,Table2,nmax,steps,ff,paramR)
{
	for (i in names(paramR)){temp <- paramR[i];storage.mode(temp) <- "double";assign(i, temp)}
	i=0:(steps-1) 
	T1=((vec[3]-vec[2])/2)*cos(pi*(i+0.5)/steps) + (vec[3]+vec[2])/2
	if(vec[1]!=nmax){T2<-((vec[3]-Table2[vec[1]+1,2])/2)*cos(pi*(i+0.5)/steps) + (vec[3]+Table2[vec[1]+1,2])/2}else{
		T2<-((vec[3]-Table2[1,2])/2)*cos(pi*(i+0.5)/steps) + (vec[3]+Table2[1,2])/2}
	
	x=T1;r1=eval(ff[[3]]);r1[r1 < 0]=0;r1[r1 > 1]=1
	R1=(sum(r1))/steps
	x=T2;r1=eval(ff[[3]]);r1[r1 < 0]=0;r1[r1 > 1]=1
	R2=(sum(r1))/steps
	Rt=(R1+R2)/2
	return(Rt)
}








#########################################################################3
###################### Metodos Antiguos #############################3##
#######################################################################
# En esta funcion ingresamos "xi" que es el valor de la edad normalizada hasta el 50% de la oviposición
Ratios<-function(Table,modelim,modelm,estadios,xi, steps){
	inmaduros <-  estadios[-(length(estadios)-1):-(length(estadios))]
	maduros   <-  estadios[(length(estadios)-1):(length(estadios))]
	matriz<-matrix(0,ncol=(length(inmaduros)*2+length(maduros)+1),nrow=nrow(Table))
	K<-1
	while(K<=length(inmaduros)){
		#  Desarrollo  # extraccion de funciones y parametros
		parametrosc <- hfeno$pdv_dv[[K]]
		parametrosc<-as.list(parametrosc)
		for (i in names(parametrosc)){
			temp <- parametrosc[[i]]
			storage.mode(temp) <- "double"
			assign(i, temp)
		}
		formulac <- hfeno$fdv_dv[[K]]
		forexc <- formulac[[length(formulac)]]
		funcionc <- as.expression(forexc)
		
		#   Mortalidad # extraccion de funciones y parametros
		parametrosm <- hfeno$pmortal[[K]]
		parametrosm<-as.list(parametrosm)
		for (i in names(parametrosm)){
			temp <- parametrosm[[i]]
			storage.mode(temp) <- "double"
			assign(i, temp)
		}
		formulam <- hfeno$mortal[[K]]
		forexm <- formulam[[length(formulam)]]
		funcionm <- as.expression(forexm)
		M<-R<-M2<-R2<-rep(0,length(Table[,1]))
		
		Ratesumme<-Rateacum<-Mortalitysumme<-Mortacum<-rep(0,length(Table[,1]))
		Ratetot1<-Mortality1<-Ratetot2<-Mortality2<-0
		
		for(i in 1:length(Table[,1])){#recorre cada valor de la tabla d temperaturas
			
			M<-(Table[,1][i]+Table[,2][i])/2
			R<-(Table[,2][i]-Table[,1][i])/2
			ifelse(i!=length(Table[,1]),M2<-(Table[,2][i]+Table[,1][i+1])/2,M2<-(Table[,2][i]+Table[,1][1])/2)
			ifelse(i!=length(Table[,1]),R2<-(Table[,2][i]-Table[,1][i+1])/2,R2<-(Table[,2][i]-Table[,1][1])/2)
			
			Ratetot1<-Mortality1<-Ratetot2<-Mortality2<-0
			#M <- round(M,1)
			#R <- round(R,1)
			#M2 <- round(M2,1)
			#R2 <- round(R2,1)
			
			
			Ratetot1<-Mortality1<-Ratetot2<-Mortality2<-0
			for(j in 0:(steps-1)){#por defecto 48 pasos
				T<-pi/steps*(j+.5)
				## eval1 ##
				DEGC<-R*cos(T)+M
				if(modelim[K]==1 || modelim[K]==2 || modelim[K]==3 ||
						modelim[K]==4 || modelim[K]==5 || modelim[K]==6 ||
						modelim[K]==7 || modelim[K]==8 || modelim[K]==9 ||
						modelim[K]==10 || modelim[K]==11 || modelim[K]==12 || modelim[K]==13 ||
						modelim[K]==14){
					DEGK<-DEGC+273.15
					x<-DEGK
				}else{ x<-DEGC }
				
				Rate<-eval(funcionc)
				Ratetot1<-Ratetot1+Rate
				x<-DEGC
				Mort<-eval(funcionm)
				if(Mort > 1){
					Mort = 1
				}else{
					Mort = Mort
				}
				Mortality1 <- Mortality1+ Mort
				
				## eval2 ##
				T<-pi/steps*(j+.5)#borrar
				DEGC<-R2*cos(T)+M2
				if(modelim[K]==1 || modelim[K]==2 || modelim[K]==3 ||
						modelim[K]==4 || modelim[K]==5 || modelim[K]==6 ||
						modelim[K]==7  || modelim[K]==8 || modelim[K]==9 ||
						modelim[K]==10 || modelim[K]==11 || modelim[K]==12 || modelim[K]==13 ||
						modelim[K]==14){
					DEGK<-DEGC+273.15
					x<-DEGK
				}else{x<-DEGC}
				
				Rate<-eval(funcionc)
				Ratetot2<-Ratetot2+Rate
				x<-DEGC
				Mort<-eval(funcionm)
				if(Mort > 1){
					Mort = 1
				}else{
					Mort = Mort
				}
				Mortality2 <- Mortality2+ Mort
				Mortality2
				
			}
			
			Ratetot1 <- Ratetot1/steps
			Ratetot2 <- Ratetot2/steps
			Ratetot<-(Ratetot1+Ratetot2)/2
			matriz[i,2*K-1] <- Ratetot
			
			Mortality1 <- Mortality1/steps
			Mortality2 <- Mortality2/steps
			Mortality<- (Mortality1+Mortality2)/2
			matriz[i,2*K] <- Mortality
			
		}
		K<-K+1
	}
	for(J in 1:(length(maduros)+1)){
		if(J==1)parametrosc <- hfeno$pfh_h
		if(J==2)parametrosc <- hfeno$pfm_m
		if(J==3)parametrosc <- hfeno$ptazaeh_h
		parametrosc<-as.list(parametrosc)
		for (i in names(parametrosc)){
			temp <- parametrosc[[i]]
			storage.mode(temp) <- "double"
			assign(i, temp)
		}
		if(J==1) formulac <- hfeno$fh_h
		if(J==2) formulac <- hfeno$fm_m
		if(J==3) formulac <- hfeno$ftazaeh_h
		forexc <- formulac[[length(formulac)]]
		funcionc <- as.expression(forexc)
		M<-R<-M2<-R2<-rep(0,length(Table[,1]))
		
		Ratesumme<-Rateacum<-rep(0,length(Table[,1]))
		ratetot1<-ratetot2<-ratetot<-0
		fectot1<-fectot2<-fectot<-fectotal<-0
		
		for(i in 1:length(Table[,1])){
			Ratetot1<-Ratetot2<-Ratetot<-0
			if(J==3){
				fectot1<-fectot2<-fectot<-0
			}
			M<-(Table[,1][i]+Table[,2][i])/2
			R<-(Table[,2][i]-Table[,1][i])/2
			ifelse(i!=length(Table[,1]),M2<-(Table[,2][i]+Table[,1][i+1])/2,M2<-(Table[,2][i]+Table[,1][1])/2)
			ifelse(i!=length(Table[,1]),R2<-(Table[,2][i]-Table[,1][i+1])/2,R2<-(Table[,2][i]-Table[,1][1])/2)
			
			for(j in 0:(steps-1)){
				T<-pi/steps*(j+.5)
				#### eval1 ###
				DEGC<-R*cos(T)+M
				x<-DEGC
				if(J==1 || J==2){
					if(modelm[J]==1 || modelm[J]==2 || modelm[J]==3 ||
							modelm[J]==4 || modelm[J]==5 || modelm[J]==6 ||
							modelm[J]==7 || modelm[J]==8 || modelm[J]==9 ||
							modelm[J]==10 || modelim[J]==11 || modelim[J]==12 || modelim[J]==13 ||
						modelim[J]==14){
						DEGK<-DEGC+273.15
						x<-DEGK
					}else{x<-DEGC}
				}
				
				Rate<-eval(funcionc)
				Ratetot1<-Ratetot1+Rate
				if(J==3){
					x<-DEGC
					fec1<-eval(funcionc)
					if(fec1>0){ fec1=fec1}else{ fec1=0}
					fectot1<-fectot1+fec1
				}
				
				#### eval2 ###
				DEGC<-R2*cos(T)+M2
				if(J==1 || J==2){
					if(modelm[J]==1 || modelm[J]==2 || modelm[J]==3 ||
							modelm[J]==4 || modelm[J]==5 || modelm[J]==6 ||
							modelm[J]==7 || modelm[J]==8 || modelm[J]==9 ||
							modelm[J]==10 || modelim[J]==11 || modelim[J]==12 || modelim[J]==13 ||
						modelim[J]==14){
						DEGK<-DEGC+273.15
						x<-DEGK
					}else{x<-DEGC}
				}
				#if(J==3) x<-DEGC
				Rate<-eval(funcionc)
				Ratetot2<-Ratetot2+Rate
				if(J==3){
					x<-DEGC
					fec2<-eval(funcionc)
					if(fec2>0){
						fec2=fec2
					}else{
						fec2=0
					}
					fectot2<-fectot2+fec2
				}
			}
			Ratetot1 <- Ratetot1/steps
			Ratetot2 <- Ratetot2/steps
			Ratetot<-(Ratetot1+Ratetot2)/2
			if(J==1){
				matriz[i, length(inmaduros)*2+1] <- Ratetot
			}
			if(J==2){
				matriz[i, length(inmaduros)*2+2] <- Ratetot
			}
			if(J==3){
				fectot1 = fectot1 / steps
				fectot2 = fectot2 / steps
				fectot = (fectot1 + fectot2) / 2
				if(fectot < 0){
					fectot = 0
				}
				matriz[i, length(inmaduros)*2+3] <- fectot
			}
		}
	}
	
	T<-1/matriz[,1]+1/matriz[,3]+1/matriz[,5]+xi/matriz[,7]
	Ro<-matriz[,9]*(1-matriz[,2])*(1-matriz[,4])*(1-matriz[,6])/2
	rm<-log(Ro)/T
	J<-exp(rm)
	Dt<-log(2)/rm
	E<-(1-matriz[,2])^matriz[,1]
	L<-(1-matriz[,4])^matriz[,3]
	P<-(1-matriz[,6])^matriz[,5]
	ERI<-(1-ifelse(length(E[E==0])==0,0,length(E[E==0])/12))*
			(1-ifelse(length(L[L==0])==0,0,length(L[L==0])/12))*
			(1-ifelse(length(P[P==0])==0,0,length(P[P==0])/12))
	
	return(list(MATRIZ=data.frame(matriz),ERI=ERI,J=J))
}
#######################################################3
parameters2 <- function (N, estad, ltb){
	estadios <- estad
	lifetable <- ltb
	s.x <- m.x <- l.x <- rep(0, nrow(lifetable))
	for (i in 1:nrow(lifetable)) l.x[i] <- sum(lifetable[i, 1:(length(estadios))])/N
	for (i in 1:nrow(lifetable)) if(l.x[i]==0) m.x[i] <- 0 else m.x[i] <- lifetable[,(length(estad)+1)][i]/N/l.x[i]
	s.x1 <- c(l.x[1], l.x[2])
	s.x2 <- rep(0, nrow(lifetable) - 2)
	for (i in 1:(nrow(lifetable) - 2)) s.x2[i] <- l.x[i + 2]/l.x[i + 1]
	s.x <- c(s.x1, s.x2)
	Ro <- sum(l.x * m.x)
	Day <- 0:(nrow(lifetable)-1)
	Tg <- sum(Day * l.x * m.x)/Ro
	GRR <- sum(m.x)
	r <- log(Ro)/Tg
	euler <- sum(l.x * m.x * exp(-r * (Day + 1)))
	r1 <- r - 1e-06
	for (k in 1:100) {
		euler1 <- sum(l.x * m.x * exp(-r1 * (Day + 1)))
		r2 <- ((1 - euler) * (r1 - r)/(euler1 - euler)) + r
		r1 <- r2
	}
	rf <- r1
	Dt <- log(2)/rf
	lamda <- exp(rf)
	T <- log(Ro)/rf
	parameters <- t(data.frame(r = rf, Ro = Ro, GRR = GRR, T = T, lambda = lamda, Dt = Dt))
	colnames(parameters) <- "Parameters"
	print(parameters)
	#return(list(m.x=m.x,s.x=s.x,l.x=as.matrix(l.x),parametro=parameters))
	return(list(r = rf, Ro = Ro, GRR = GRR, T = T, lambda = lamda, Dt = Dt))
}

