######################################################################
pr.c <- function(N,NH,Fec)
{
 ftp = y ~ 1 - exp(-(Fec*NH)/N)  ## creo que esta tasa de parasitacion es con respecto a una poblacion finita respecto a todos los parasitoides
 ptp = data.frame(Fec=Fec,NH=NH,N=N)
 modelp = 60
 return(list(ftp=ftp,ptp=ptp,modelp=modelp))
}
#######################################################################
# Calcula los parametros biologicos segun la temperatura en 2 especies
simultemp.2spec<-function(N,sexratio,isFixed,temp,xi,steps,ftp,ptp,modelp,NH,est.risk,..){
	n1=length(temp);n2=2*est.risk
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

		cuadro <- matrizA(estadios, hfeno, Table, steps)  ## solo se obtienen fits del dia en cuestion

		##esto lo inclui xq no existe el valor de M
		M<-nrow(Table)
		
		if(modelp==60)
		{
			PP=rep(pars.rate(temp[i],ftp,ptp,modelp),M) ## calcula directamente la probabilidad de parasitizacion
		}else{
			## ingreso la tasa de parasitacion
			tazpar=rep(pars.rate(temp[i],ftp,ptp,modelp),M) ## se calcula la tasa de parasitacion para la iesima temperatura

			## manipula la temperatura y la taza de sexo

			PP=(tazpar/cuadro[,n2-1])*(NH/N);PP[PP>1]=1  ## hallando el porcentaje de parasitacion
		}						#
		mort2=1-(1-cuadro[,n2])*(1-PP) ## nueva mortalidad modificada por la parasitacion
		cuadro[,n2]=mort2
		ageclases <- AgeClases(cuadro, estadios, hfeno)
		estadiosAges <- ageclases$estadiosAges
		oviFreq <- ageclases$oviFreq
		ageClassMatriz <- ageclases$ageClassMatriz
		Day=1
		severalYear=FALSE
		Steps=364
		simu2<-simulacionUnaGeneracion(Day,estadiosAges, oviFreq, ageClassMatriz, hfeno, cuadro, estadios,N,sexratio,Rs,Steps)
		matrizOut<-simu2$matrizOut
			
		pravida<-parameters(N, estadios, matrizOut)
		parametros<-pravida$parametro
		if(!is.nan(parametros[1])){if(parametros[1]<0){parametros[c(1,6)]=NA}}
		pars=rbind(pars,c(Temper=temp[i],t(parametros)))
	}
	pars=pars[-1,]
	return(pars)
}
###########################################################
# Predice las tasas de parasitacion para cierta temperatura
pars.rate<-function(x,ftp,ptp,modelp)
{
 x=mean(x) ## para asegurar que tambien se pueda usar vectores como imput cuando hay minimo y maximo
 Far=1:14;n1=length(Far[Far==modelp])
 if(n1!=0){x=x+273.15}
 parametros <- ptp
 for(i in names(parametros)){temp <- parametros[i];storage.mode(temp) <- "double";assign(i, temp)}
 funcion <- as.expression(ftp)
 pr=eval(funcion);pr[pr<=0]=0;pr=pr[1,1]
 return(pr)
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
				modelim[Ki]==10 || modelim[Ki]==11 || modelim[Ki]==12 || modelim[Ki]==13 ||
				modelim[Ki]==14) & is.na(J)){x = T1 + 273.15;x2 = T2 + 273.15}else{x = T1;x2 = T2}


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
#####################################################################
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
################################################################
simulacionUnaGeneracion <- function(Day, estadiosAges, oviFreq, ageClassMatriz, hfeno, cuadro, estadios,numIni,sexratio,Rs,Steps){
	Step=0

	matrizA = cuadro

	sizeMatriz = nrow(matrizA)

	inmad <-  estadios[-(length(estadios)-1):-(length(estadios))]

	numInmaduros = length(inmad)

	maduros   <-  estadios[(length(estadios)-1):(length(estadios))]

	numMaduros = length(maduros)

	NEWIND = 0

	dias<-254#porq?



	#matrizOut = matrix(NA, nrow=(nrow(matrizA)-1), (numInmaduros+numMaduros+1))

	matrizOut = matrix(NA, nrow=nrow(matrizA), (numInmaduros + numMaduros + 1))

	stagesN = c(rep(0,dias))

	stagesN1 = c(numIni, rep(0,(dias-1)))

	femalesN = c(rep(NA,dias))

	eggsN = stagesN1

	stagesDev = c(rep(0,dias))



	## vectorN <- list(1)

	## vectorNini <- list(1)



	#if(!severalYear){

	vectorN <- list(1)

	vectorNini <- list(1)

	for(k in 1:(numInmaduros + numMaduros)){

		v1 = c(100,rep(0,(ageClassMatriz[1,k]-1)))

		v2 = c(0,rep(0,(ageClassMatriz[1,k]-1)))

		v3 = c(0,rep(0,(ageClassMatriz[1,k]-1)))#



		if(k == 1){

			vectorN[[k]] = list(v1,v2,v3)#

			vectorNini[[k]] = v1

		}else{

			vectorN[[k]] = list(v2,v2,v3)#

			vectorNini[[k]] = v2

		}

	}



	matrizOut[1,] = c(numIni,rep(0,(numInmaduros + numMaduros)))# numero inicial de individuos

	#}

	## else{

	##     vectorNini <- list(1)

	##     for(k in 1:(numInmaduros + numMaduros)){

	##         vectorNini[[k]] <- vectorN[[k]][[1]]

	##         matrizOut[1,k] = round(numInd[[k]])

	##     }

	## }



	#while(Day < (sizeMatriz-1)){

	while(Step < 363){

		for(k in 1:(numInmaduros + numMaduros)){

			#Ageclass = ageClassMatriz[Day,k]

			Ageclass = ageClassMatriz[Day+1,k]



			if(k <= numInmaduros){

				stagesN = c(rep(0,Ageclass))

			}



			if(k == (numInmaduros+1)){

				femalesN = c(rep(0,Ageclass))

			}



			matrizAges = estadiosAges[[k]]

			NEWIND = 0



			E <- c(rep(0, Ageclass))

			EL <- c(rep(0, Ageclass))

			Eage <- c(rep(0, Ageclass))

			pacc <- c(rep(0, Ageclass))

			p1 <- c(rep(0, Ageclass))

			p2 <- c(rep(0, Ageclass))



			if(k <= numInmaduros){

				s1 <- c(rep(NA, Ageclass))

			}



			if(k == numInmaduros+1){

				vectorN[[k]][[2]] = vectorN[[k]][[3]]#

			}

			if(k == numInmaduros+2){

				vectorN[[k]][[2]] = vectorN[[k]][[3]]#

			}



			for(Ii in 1:(Ageclass-1)){

				Eage[Ii] = matrizAges[Day+1, Ii+1]



				if(k <= numInmaduros){

					survivalk = (1 - matrizA[Day,k*2]) ^ matrizA[Day,k*2-1]

				}



				if(k <= numInmaduros){

					if(Eage[Ii] > 1){

						s1[Ii] = 1

					}else{

						s1[Ii] = survivalk

					}

				}



				if(k <= numInmaduros){

					pacc[Ii] = distrimodeldeim(Eage[Ii], k)

				}else{

					pacc[Ii] = distrimodeldema(Eage[Ii], k - numInmaduros)

				}



				if(Ii == 1){

					p1[Ii] = pacc[Ii] - 0

				}else{

					p1[Ii] = pacc[Ii] - pacc[Ii - 1]

				}



				if(p1[Ii] == 0){

					p2[Ii] = 0

				}else{

					if(Ii == 1){

						p2[Ii] = p1[Ii] / (1 - 0)

					}else{

						p2[Ii] = p1[Ii] / (1 - pacc[Ii - 1])

					}

				}



				if(Day==1){

					E[Ii] = vectorNini[[k]][Ii]

				}else{

					E[Ii] = vectorN[[k]][[1]][Ii]#modifique el 1 por el 2

				}





				if(k <= numInmaduros){

					EL[Ii] = E[Ii] * p2[Ii] * s1[Ii]

				}else{

					EL[Ii] = E[Ii] * p2[Ii]

				}



				if(k <= numInmaduros){

					E1 = E[Ii] * s1[Ii] - EL[Ii]

				}else{

					E1 = E[Ii] - EL[Ii]

				}



				stagesN[Ii+1] = E1



				if(k <= numInmaduros){

					NEWIND = NEWIND + EL[Ii]

				}



			} # fin for clases



			vectorN[[k]][[2]] = c(vectorN[[k]][[2]][1], stagesN[2:length(stagesN)])

			vectorN[[k]][[1]] = vectorN[[k]][[2]]



			if(k <= (numInmaduros-1)){

				vectorN[[k]][[3]][1]=NEWIND

				vectorN[[k+1]][[2]] = vectorN[[k]][[3]]

			}



			if(k == numInmaduros){

				vectorN[[k+1]][[3]][1] = NEWIND * Rs[Day]

				vectorN[[k+2]][[3]][1] = NEWIND * (1-Rs[Day])



				vectorN[[k+1]][[2]] = vectorN[[k]][[3]]#

			}



			if(k == (numInmaduros+1)){

				femalesN = E

			}



			matrizOut[Day+1,k] = round(sum(vectorN[[k]][[2]]))

			#matrizOut[1:3,]



		}#fin for estadios





		#Oviposition

		Ageclass = ageClassMatriz[Day,numInmaduros+1]

		matrizAges = estadiosAges[[numInmaduros+1]]

		Ovi = c(rep(0,Ageclass))



		E <- c(rep(0, Ageclass))

		EL <- c(rep(0, Ageclass))

		Eage <- c(rep(0, Ageclass))

		pacc <- c(rep(0, Ageclass))

		p1 <- c(rep(0, Ageclass))

		p2 <- c(rep(0, Ageclass))



		for(Ii in 1:(Ageclass-1)){

			E[Ii] = vectorN[[k-1]][[2]][Ii]

			Eage[Ii] = oviFreq[Day, Ii]

			E1 = matrizA[Day,k*2-1]

			E1 = E1 * Eage[Ii] * E[Ii]

			Ovi[Ii+1] = E1

			NEWIND = NEWIND + E1

		}



		Ovi = na.omit(Ovi)

		matrizOut[Day+1,(numInmaduros + numMaduros +1)] = round(sum(Ovi),2)



		#vectorN[[1]][[2]][1] = NEWIND

		matrizOut[Day+1,1] = round(sum(vectorN[[1]][[2]]))

		vectorN[[1]][[1]] = vectorN[[1]][[2]]



#                      if(round(matrizOut[Day,(numInmaduros + numMaduros)]) > 0 && round(matrizOut[Day+1,(numInmaduros + numMaduros)]) == 0){

#                                  break;

#                      }

		if(matrizOut[Day,1]==0 && matrizOut[Day,2]==0 && matrizOut[Day,3]==0 && matrizOut[Day,4]==0 && matrizOut[Day,5]==0){

			#matrizOut[Day:nrow(matrizOut),]="Dead"
			break;

		}



		NEWIND=0

		Day=Day+1

		Step = Step + 1

	}#fin bucle steps





	namesMatriz = c(estadios, "New Egg");

	colnames(matrizOut) = namesMatriz

	matrizOutDead = matrizOut
	matrizOutDead[is.na(matrizOutDead)]="Dead"

	matrizOut[is.na(matrizOut)]=0

	return(list(matrizOut=matrizOut, vectorN=vectorN, matrizOutDead=matrizOutDead))







}
###############################################################################
###############################################################################

distrimodeldeim<-function(vec,sll) ## this is GML function
{
  # probit
  if(hfeno$distri_dv[[sll]]=="probit")
  {
    pdd<-pnorm(log(vec)*hfeno$slope_dv[[sll]])
  }
  # logit
  if(hfeno$distri_dv[[sll]]=="logit")
  {
    pdd<-1/(1+exp(-(log(vec)*hfeno$slope_dv[[sll]])))
  }
  # cloglog
  if(hfeno$distri_dv[[sll]]=="cloglog")
  {
    pdd<-1-exp(-exp(log(vec)*hfeno$slope_dv[[sll]]))
  }
return(pdd)
}
###############################################################################
###############################################################################
distrimodeldema<-function(vec,sll)
{
  if(sll==1)
  {
    #hembras
    # probit
    if(hfeno$distri_snh=="probit")
    {
      pdd<-pnorm(log(vec)*hfeno$slope_snh)
    }
    # logit
    if(hfeno$distri_snh=="logit")
    {
      pdd<-1/(1+exp(-(log(vec)*hfeno$slope_snh)))
    }
    # cloglog
    if(hfeno$distri_snh=="cloglog")
    {
      pdd<-1-exp(-exp(log(vec)*hfeno$slope_snh))
    }
  }
  if(sll==2)
  {
    #machos
    # probit
    if(hfeno$distri_snm=="probit")
    {
      pdd<-pnorm(log(vec)*hfeno$slope_snm)
    }
    # logit
    if(hfeno$distri_snm=="logit")
    {
      pdd<-1/(1+exp(-(log(vec)*hfeno$slope_snm)))
    }
    # cloglog
    if(hfeno$distri_snm=="cloglog")
    {
      pdd<-1-exp(-exp(log(vec)*hfeno$slope_snm))
    }
  }
return(pdd)
}
#######################################################3
parameters <- function (N, estad, ltb){
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
	return(list(m.x=m.x,s.x=s.x,l.x=as.matrix(l.x),parametro=parameters))
}
###########################
graphsimulp<-function(pars,modls=c(1,1,1,1,1,1),direc,ecua,nomec,ejex,ejey,tecu,grises,tit,ax,ay){

                if(grises==TRUE){ccol=c("gray20","gray50")}else{ccol=c("royalblue","red")}



                etip=c("Rm","Ro","GRR","GL","Lambda","Dt")

                nfile=etip[!is.na(modls)]

                NMS=paste(nfile,collapse="_")

                pathIm = paste(direc,NMS,".jpg",sep="")

                jpeg(pathIm, width = 720, height = 800,quality = 100)##

                # Validadando si se tiene mas observaciones que parametros
                ndatos<-function(vec){n1=length(vec)-length(vec[is.na(vec)]); if(length(vec[vec==0])==n1){n1=0}; return(n1)}
                ndts=apply(pars[,-1],2,ndatos)
                ind1=rep(nrow(pars),6)
                ind1[modls!=1]=3;ind1[modls==1]=4
                modls[ndts<=ind1]=NA

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





                                x=pars[,1]

                                y=pars[,k]



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
##############################################################################
# Calcula los parametros biologicos segun temperatura fluctuante en 2 especies

simultemp.2spec.fluc<-function(numIni,estadios,hfeno,params,sexratio,isFixed,Table,steps,ftp,ptp,modelp,NH,est.risk,multi.gen=FALSE)
{
	n2=2*est.risk
	if(isFixed){
		Rs=rep(sexratio,365)
	}else{
		Table2=cbind(id=1:nrow(Table),Table)
		paramR=params$paramR
		ff=params$ff
		Rs=apply(Table2,1,RateR,Table2,nmax=nrow(Table),steps=steps,ff,paramR)
	}
	cuadro <- matrizA(estadios, hfeno, Table, steps)

	if(modelp==60)
	{
		PP=apply(Table,1,pars.rate,ftp,ptp,modelp)

	}else{
		## ingreso la tasa de parasitacion
		tazpar=apply(Table,1,pars.rate,ftp,ptp,modelp)
		PP=(tazpar/cuadro[,n2-1])*(NH/numIni);PP[PP>1]=1
	}						#

	mort2=1-(1-cuadro[,n2])*(1-PP) ## nueva mortalidad modificada por la parasitacion
	cuadro[,n2]=mort2

	ageclases <- AgeClases(cuadro, estadios, hfeno)
	estadiosAges <- ageclases$estadiosAges
	oviFreq <- ageclases$oviFreq
	ageClassMatriz <- ageclases$ageClassMatriz
	Day <- 1
	severalYear=FALSE
	Steps=364
	if(multi.gen)
	{
		simu2<-simulation2(Day,estadiosAges, oviFreq, ageClassMatriz, hfeno, cuadro, estadios,numIni,sexratio,Rs,Steps,severalYear) ## en la funcion simulation2 se agrega el parametro "severalYear"
		matrizOut<-simu2$matrizOut
		matrizOut[is.na(matrizOut)]=0
		matrizOutDead<-simu2$matrizOutDead
		pravida<-parameters(numIni, estadios, matrizOut)
		parametros<-pravida$parametro
	 	return(list(matrizOut=matrizOut,matrizOutDead=matrizOutDead,cuadro=cuadro))
	}else{
		simu2<-simulacionUnaGeneracion(Day,estadiosAges, oviFreq, ageClassMatriz, hfeno, cuadro, estadios,numIni,sexratio,Rs,Steps)
		matrizOut<-simu2$matrizOut
		matrizOutDead<-simu2$matrizOutDead
		pravida<-parameters(numIni, estadios, matrizOut)
		parametros<-pravida$parametro

	 	return(list(parametros=parametros,matrizOut=matrizOut,matrizOutDead=matrizOutDead,cuadro=cuadro))
	}

}
#################################################################
simultemp.2spec.fluc.2<-function(matrizOut1,matrizOutP,numIni,estadios,hfeno,params,sexratio,isFixed,Table,steps,ftp,ptp,modelp,NH,est.risk,multi.gen=FALSE,FecunT=FALSE)
{
	n2=2*est.risk
	n3=(1:length(estadiosP))[length(estadiosP)-1]
	if(isFixed){
		Rs=rep(sexratio,365)
	}else{
		Table2=cbind(id=1:nrow(Table),Table)
		paramR=params$paramR
		ff=params$ff
		Rs=apply(Table2,1,RateR,Table2,nmax=nrow(Table),steps=steps,ff,paramR)
	}
	cuadro <- matrizA(estadios, hfeno, Table, steps)

	N=matrizOut1[,est.risk]+0.0001
	matrizOutP[1,n3]=NH
	NH=matrizOutP[,n3]

	if(modelp==60)
	{
		if(FecunT){PP=matrizOutP[,length(estadiosP)+1]/N;PP[PP>1]=1}else  ## si cumple la condicion se esta usando la oviposicion efectiva del parasitoide
		{
			Fec=ptp$Fec
			PP=eval(ftp);PP[PP>1]=1  ## la fecundidad "Fec" se queda constante mientras que los demas parametros varian por dia
		}



	}else{
		## ingreso la tasa de parasitacion
		tazpar=apply(Table,1,pars.rate,ftp,ptp,modelp) ## la fecundidad de la hembra parasitoide es variable segun la temperatura
		PP=tazpar*NH/N;PP[PP>1]=1
	}						

	mort2=1-(1-cuadro[,n2])*(1-PP) ## nueva mortalidad modificada por la parasitacion
	cuadro[,n2]=mort2

	ageclases <- AgeClases(cuadro, estadios, hfeno)
	estadiosAges <- ageclases$estadiosAges
	oviFreq <- ageclases$oviFreq
	ageClassMatriz <- ageclases$ageClassMatriz
	Day <- 1
	severalYear=FALSE
	Steps=364
	if(multi.gen)
	{
		simu2<-simulation2(Day,estadiosAges, oviFreq, ageClassMatriz, hfeno, cuadro, estadios,numIni,sexratio,Rs,Steps,severalYear) ## en la funcion simulation2 se agrega el parametro "severalYear"
		matrizOut<-simu2$matrizOut
		matrizOut[is.na(matrizOut)]=0
		matrizOutDead<-simu2$matrizOutDead
		pravida<-parameters(numIni, estadios, matrizOut)
		parametros<-pravida$parametro
	 	return(list(matrizOut=matrizOut,matrizOutDead=matrizOutDead,cuadro=cuadro))
	}else{
		simu2<-simulacionUnaGeneracion(Day,estadiosAges, oviFreq, ageClassMatriz, hfeno, cuadro, estadios,numIni,sexratio,Rs,Steps)
		matrizOut<-simu2$matrizOut
		matrizOutDead<-simu2$matrizOutDead
		pravida<-parameters(numIni, estadios, matrizOut)
		parametros<-pravida$parametro

	 	return(list(parametros=parametros,matrizOut=matrizOut,matrizOutDead=matrizOutDead,cuadro=cuadro))
	}

}
######################################################################################################
# grafica los conteos de la tabla de vida simulada por dia y por estado para una y varias generaciones

grafSimDete<-function(matrizOut, estadios, labx,laby,titulo,legx,legy,corrx1,corrx2,logscale=FALSE,see.male=FALSE){
	numEstadios = length(estadios)
	names1 = colnames(matrizOut)
	names1=names1[1:(length(names1)-2)]
	if(logscale){dat = log(matrizOut + 1)}else{dat = matrizOut}
	corrx=c(corrx1,corrx2)
	x1 = 1:nrow(dat)
	z1 = dat[,numEstadios-1]
	mod1 = lm(z1~x1)

	rang1=c(0.1,1,5,10,20,50,100,1000,10e+3,10e+3,10e+4,10e+5,10e+6,10e+7,10e+8,10e+9,10e+10,10e+11,10e+12,10e+13)
        if(is.null(corrx)){corrx=c(min(x1,0),max(x1)+max(x1)*0.2);corrx2=seq(0,10*round(corrx[2]/10,1),5)}else{corrx2=seq(corrx[1],corrx[2],5)} ## cambio
	corry=c(0,max(dat[,1])+max(dat[,1])*0.2);corry2=seq(0,round(max(corry),1),rang1[order(sqrt((rang1-round(max(corry),1)*0.1)^2))[1]])

	plot(x1,dat[,1],frame=F,type="n",pch=19,xlim=corrx,ylim=corry,xlab=labx,ylab=laby,axes=F,xaxt = "n", main=titulo) ## 4
	axis(1, corrx2,cex.axis=0.7)  ## cambio
	if(logscale){temp2=data.frame(exp(corry2));corry2=apply(temp2,1,val.dist.euc,rang1)
	axis(2, at=log(corry2),lab=corry2,las=2,cex.axis=0.6);abline(h = log(corry2), v = 0, lty = 3, lwd = .1, col = "gray78")}else{
	axis(2, corry2,las=2,cex.axis=0.6);abline(h = corry2, v = 0, lty = 3, lwd = .1, col = "gray78")}
	abline(mod1, lwd=2)

	if(see.male){n2=numEstadios}else{(n2=numEstadios-1)}
	for(i in 1:n2){
		points(x1,dat[,i],type="l",col=(i+1))
	}

	cols = c(2:(i+1),1)
	legend(legx,legy,c(names1,paste("Expon.(",names1[length(names1)],")")),col =cols,lty = 1)
}
#############################################################
simulation2<-function(Day, estadiosAges, oviFreq, ageClassMatriz, hfeno, cuadro, estadios,numIni,sexratio,Rs, Steps,severalYear){
                #Step=0

                matrizA = cuadro

                sizeMatriz = nrow(matrizA)

                inmad <-  estadios[-(length(estadios)-1):-(length(estadios))]

                numInmaduros = length(inmad)

                maduros   <-  estadios[(length(estadios)-1):(length(estadios))]

                numMaduros = length(maduros)

                NEWIND = 0

                dias<-254#porq?



                matrizOut = matrix(NA, nrow=(nrow(matrizA)), (numInmaduros+numMaduros+1))

                #matrizOut[1,1] = numIni# numero inicial de individuos

                stagesN = c(rep(0,dias))

                stagesN1 = c(numIni, rep(0,(dias-1)))

                femalesN = c(rep(NA,dias))

                eggsN = stagesN1

                stagesDev = c(rep(0,dias))

                Step=0



                ## vectorN <- list(1)

                ## vectorNini <- list(1)



                if(!severalYear){

                                vectorN <- list(1)

                                vectorNini <- list(1)

                                for(k in 1:(numInmaduros + numMaduros)){

                                                #v1 = c(100,rep(0,(ageClassMatriz[1,k]-1)))
                                                v1 = c(numIni,rep(0,(ageClassMatriz[1,k]-1)))####################

                                                v2 = c(0,rep(0,(ageClassMatriz[1,k]-1)))

                                                v3 = c(0,rep(0,(ageClassMatriz[1,k]-1)))#



                                                if(k == 1){

                                                                vectorN[[k]] = list(v1,v2,v3)#

                                                                vectorNini[[k]] = v1

                                                }else{

                                                                vectorN[[k]] = list(v2,v2,v3)#

                                                                vectorNini[[k]] = v2

                                                }

                                }



                                matrizOut[1,1] = numIni# numero inicial de individuos
                                matrizOut[1,2:ncol(matrizOut)] = 0#####################

                }else{

                                vectorNini <- list(1)

                                for(k in 1:(numInmaduros + numMaduros)){

                                                vectorNini[[k]] <- vectorN[[k]][[1]]

                                                matrizOut[1,k] = round(numInd[[k]])

                                }

                }



                #while(Day < (sizeMatriz-1)){

                while(Step < Steps){

                                for(k in 1:(numInmaduros + numMaduros)){

                                                Ageclass = ageClassMatriz[Day+1,k]



                                                if(k <= numInmaduros){

                                                                stagesN = c(rep(0,Ageclass))

                                                }



                                                if(k == (numInmaduros+1)){

                                                                femalesN = c(rep(0,Ageclass))

                                                }



                                                matrizAges = estadiosAges[[k]]

                                                NEWIND = 0



                                                E <- c(rep(0, Ageclass))

                                                EL <- c(rep(0, Ageclass))

                                                Eage <- c(rep(0, Ageclass))

                                                pacc <- c(rep(0, Ageclass))

                                                p1 <- c(rep(0, Ageclass))

                                                p2 <- c(rep(0, Ageclass))



                                                if(k <= numInmaduros){

                                                                s1 <- c(rep(NA, Ageclass))

                                                }



                                                if(k == numInmaduros+1){

                                                                vectorN[[k]][[2]] = vectorN[[k]][[3]]#

                                                }

                                                if(k == numInmaduros+2){

                                                                vectorN[[k]][[2]] = vectorN[[k]][[3]]#

                                                }



                                                for(Ii in 1:(Ageclass-1)){

                                                                Eage[Ii] = matrizAges[Day+1, Ii+1]



                                                                if(k <= numInmaduros){

                                                                                survivalk = (1 - matrizA[Day,k*2]) ^ matrizA[Day,k*2-1]

                                                                }



                                                                if(k <= numInmaduros){

                                                                                if(Eage[Ii] > 1){

                                                                                                s1[Ii] = 1

                                                                                }else{

                                                                                                s1[Ii] = survivalk

                                                                                }

                                                                }



                                                                if(k <= numInmaduros){

                                                                                pacc[Ii] = distrimodeldeim(Eage[Ii], k)

                                                                }else{

                                                                                pacc[Ii] = distrimodeldema(Eage[Ii], k - numInmaduros)

                                                                }



                                                                if(Ii == 1){

                                                                                p1[Ii] = pacc[Ii] - 0

                                                                }else{

                                                                                p1[Ii] = pacc[Ii] - pacc[Ii - 1]

                                                                }



                                                                if(p1[Ii] == 0){

                                                                                p2[Ii] = 0

                                                                }else{

                                                                                if(Ii == 1){

                                                                                                p2[Ii] = p1[Ii] / (1 - 0)

                                                                                }else{

                                                                                                p2[Ii] = p1[Ii] / (1 - pacc[Ii - 1])

                                                                                }

                                                                }



                                                                if(Day==1){

                                                                                E[Ii] = vectorNini[[k]][Ii]

                                                                }else{

                                                                                E[Ii] = vectorN[[k]][[1]][Ii]#modifique el 1 por el 2

                                                                }





                                                                if(k <= numInmaduros){

                                                                                EL[Ii] = E[Ii] * p2[Ii] * s1[Ii]

                                                                }else{

                                                                                EL[Ii] = E[Ii] * p2[Ii]

                                                                }



                                                                if(k <= numInmaduros){

                                                                                E1 = E[Ii] * s1[Ii] - EL[Ii]

                                                                }else{

                                                                                E1 = E[Ii] - EL[Ii]

                                                                }



                                                                stagesN[Ii+1] = E1



                                                                if(k <= numInmaduros){

                                                                                NEWIND = NEWIND + EL[Ii]

                                                                }



                                                } # fin for clases



                                                vectorN[[k]][[2]] = c(vectorN[[k]][[2]][1], stagesN[2:length(stagesN)])

                                                vectorN[[k]][[1]] = vectorN[[k]][[2]]



                                                if(k <= (numInmaduros-1)){

                                                                vectorN[[k]][[3]][1]=NEWIND

                                                                vectorN[[k+1]][[2]] = vectorN[[k]][[3]]

                                                }



                                                if(k == numInmaduros){

                                                                vectorN[[k+1]][[3]][1] = NEWIND * Rs[Day]

                                                                vectorN[[k+2]][[3]][1] = NEWIND * (1-Rs[Day])



                                                                vectorN[[k+1]][[2]] = vectorN[[k]][[3]]#

                                                }



                                                if(k == (numInmaduros+1)){

                                                                femalesN = E

                                                }



                                                matrizOut[Day+1,k] = round(sum(vectorN[[k]][[2]]))

                                                #matrizOut[1:3,]



                                }#fin for estadios



                                #Oviposition

                                Ageclass = ageClassMatriz[Day,numInmaduros+1]

                                matrizAges = estadiosAges[[numInmaduros+1]]

                                Ovi = c(rep(0,Ageclass))



                                E <- c(rep(0, Ageclass))

                                EL <- c(rep(0, Ageclass))

                                Eage <- c(rep(0, Ageclass))

                                pacc <- c(rep(0, Ageclass))

                                p1 <- c(rep(0, Ageclass))

                                p2 <- c(rep(0, Ageclass))



                                for(Ii in 1:(Ageclass-1)){

                                                E[Ii] = vectorN[[k-1]][[2]][Ii]

                                                Eage[Ii] = oviFreq[Day, Ii]

                                                E1 = matrizA[Day,k*2-1]

                                                E1 = E1 * Eage[Ii] * E[Ii]

                                                Ovi[Ii+1] = E1

                                                NEWIND = NEWIND + E1

                                }



                                Ovi = na.omit(Ovi)

                                matrizOut[Day+1,(numInmaduros + numMaduros +1)] = round(sum(Ovi),2)



                                vectorN[[1]][[2]][1] = NEWIND

                                matrizOut[Day+1,1] = round(sum(vectorN[[1]][[2]]))

                                vectorN[[1]][[1]] = vectorN[[1]][[2]]



                                NEWIND=0

                                Day=Day+1

                                Step = Step + 1

                }#fin bucle steps



                namesMatriz = c(estadios, "New Egg");

                colnames(matrizOut) = namesMatriz

                matrizOutDead = matrizOut
                #matrizOutDead[is.infinite(matrizOutDead)]="Dead"
                matrizOutDead[is.na(matrizOutDead)]="Dead"
                matrizOutDead[matrizOutDead == "Inf"]="Dead"
                #matrizOutDead[is.nan(matrizOutDead)]="Dead"

                return(list(matrizOut=matrizOut, vectorN=vectorN, matrizOutDead=matrizOutDead))




}
######################
# distancio euclidiana

val.dist.euc<-function(x,rang){rang[order(sqrt((rang-x)^2))[1]]}

#################################################################################################################
# grafica los conteos de la tabla de vida simulada por dia de varias generaciones del Hospedero y del Parasitoide

grafSimDete.2esp<-function(matrizOut,matrizOutP,estadios,estadiosP,labx,laby,titulo,lgx1,lgy1,lgx2,lgy2,corrx1,corrx2,logscale=FALSE,see.host,see.pars){

	numEstadios = length(estadios);numEstadiosP = length(estadiosP)
	if(logscale){dat = log(matrizOut + 1);dat2 = log(matrizOutP + 1)}else{dat = matrizOut;dat2 = matrizOutP}
	corrx0=c(corrx1,corrx2)
	x1 = 1:nrow(dat)
	z1 = dat[,numEstadios-1]
	z2 = dat2[,numEstadiosP-1]

	mod1 = lm(z1~x1)
	mod2 = lm(z2~x1)
	if(max(dat[,1]) > max(dat2[,1])){ y1= dat[,1]}else{y1= dat2[,1]}

	rang1=c(0.1,1,5,10,20,50,100,1000,10e+3,10e+3,10e+4,10e+5,10e+6,10e+7,10e+8,10e+9,10e+10,10e+11,10e+12,10e+13)
        if(is.null(corrx0)){corrx=c(min(x1,0),max(x1)+max(x1)*0.2);corrx2=seq(0,10*round(corrx[2]/10,1),5)}else{corrx2=seq(corrx0[1],corrx0[2],5)} ## cambio
	#corry=c(0,max(dat[,1])+max(dat[,1])*0.2);corry2=seq(0,round(max(corry),1),rang1[order(sqrt((rang1-round(max(corry),1)*0.1)^2))[1]])
	corry=c(0,c(max(dat[,1],dat2[,1])+max(dat[,1],dat2[,1])*0.2));corry2=seq(0,round(max(corry),1),rang1[order(sqrt((rang1-round(max(corry),1)*0.1)^2))[1]])
	plot(x1,y1,frame=F,type="n",pch=19,xlim=corrx0,ylim=corry,xlab=labx,ylab=laby,axes=F,xaxt = "n", main=titulo)
	axis(1, corrx2,cex.axis=0.7)  ## cambio
	if(logscale){temp2=data.frame(exp(corry2));corry2=apply(temp2,1,val.dist.euc,rang1)
	axis(2, at=log(corry2),lab=corry2,las=2,cex.axis=0.6);abline(h = log(corry2), v = 0, lty = 3, lwd = .1, col = "gray78")}else{
	axis(2, corry2,las=2,cex.axis=0.6);abline(h = corry2, v = 0, lty = 3, lwd = .1, col = "gray78")}
	abline(mod1, lwd=2)
	abline(mod2, lwd=2,col="gray50")

	# puntos del Host	
	estadH=data.frame(1:length(estadios));rownames(estadH)=estadios
	n2=estadH[see.host,1]
	col1=brewer.pal(length(estadios)+2,"Greens");col1=col1[-(1:2)]

	for(i in n2){
		points(x1,dat[,i],type="l",col=col1[i])
	}

	# puntos del Parasitoide
	estadP=data.frame(1:length(estadiosP));rownames(estadP)=estadiosP
	n3=estadP[see.pars,1]
	col2=brewer.pal(length(estadiosP)+2,"Reds");col2=col2[-(1:2)]

	for(i in n3){
		points(x1,dat2[,i],type="l",col=col2[i])
	}
	
	# leyenda del Host
	legend(lgx1,lgy1,c(estadios[n2],paste("Expon.(",estadios[length(estadios)-1],")")),col =c(col1[n2],"gray20"),lty = 1,lwd=2.5,cex=0.65,title="HOST")

	# leyenda del Parasitoid
	legend(lgx2,lgy2,c(estadiosP[n3],paste("Expon.(",estadiosP[length(estadiosP)-1],")")),col =c(col2[n3],"gray50"),lty = 1,lwd=2.5,cex=0.65,title="PARASITOID")
}


##################################################################################################################################
# grafica los conteos de la tabla de vida simulada por dia de varias generaciones del Hospedero con varias temperaturas constantes

grafSimDete.2esp.Temps<-function(matrizOut,estadios,labx,laby,titulo,lgx,lgy,corrx1,corrx2,logscale=FALSE,see.host,temprs,see.temprs){

	numEstadios = length(estadios)
	names1 = colnames(matrizOut)
	names1=names1[1:(length(names1)-2)]
	if(logscale){dat = log(matrizOut + 1)}else{dat = matrizOut}
	corrx=c(corrx1,corrx2)
	ntt=(nrow(dat)/length(temprs))
	x1 = 1:ntt
	#z1 = dat[,numEstadios-1]
	#mod1 = lm(z1~x1)
	Temprs=rep(temprs,rep(ntt,length(temprs)))
	dat=data.frame(Temprs,dat)

	rang1=c(0.1,1,5,10,20,50,100,1000,10e+3,10e+3,10e+4,10e+5,10e+6,10e+7,10e+8,10e+9,10e+10,10e+11,10e+12,10e+13)
        if(is.null(corrx)){corrx=c(min(x1,0),max(x1)+max(x1)*0.2);corrx2=seq(0,10*round(corrx[2]/10,1),5)}else{corrx2=seq(corrx[1],corrx[2],5)} ## cambio
	corry=c(0,max(dat[,2])+max(dat[,2])*0.1);corry2=seq(0,round(max(corry),1),rang1[order(sqrt((rang1-round(max(corry),1)*0.1)^2))[1]])

	plot(c(x1[1],x1[length(x1)]),c(min(dat[,2]),max(dat[,2])),frame=F,type="n",pch=19,xlim=corrx,ylim=corry,xlab=labx,ylab=laby,axes=F,xaxt = "n", main=titulo) ## 4
	axis(1, corrx2,cex.axis=0.7)  ## cambio
	if(logscale){temp2=data.frame(exp(corry2));corry2=apply(temp2,1,val.dist.euc,rang1)
	axis(2, at=log(corry2),lab=corry2,las=2,cex.axis=0.6);abline(h = log(corry2), v = 0, lty = 3, lwd = .1, col = "gray78")}else{
	axis(2, corry2,las=2,cex.axis=0.6);abline(h = corry2, v = 0, lty = 3, lwd = .1, col = "gray78")}
	#abline(mod1, lwd=2)


	# estados para ver
	estadH=data.frame(1:length(estadios));rownames(estadH)=estadios
	n2=estadH[see.host,1]

	# estados para ver
	temprsH=data.frame(1:length(temprs));rownames(temprsH)=temprs
	n3=temprsH[as.character(see.temprs),1]

	orden.col=c(35, 19, 29, 23, 20, 27, 21, 22, 30, 24, 33, 31, 25, 18, 28, 32, 26, 34)
	colores=rownames(brewer.pal.info[orden.col,])  ## este objeto ya existe cuando se activa la libreria gráfica

	col.leg=rep(NA,length(estadios))
	for(j in n3)
	{
		col1=brewer.pal(length(estadios)+2,colores[j]);col1=col1[-(1:2)];col.leg[j]=col1[length(col1)]
		for(i in n2){
			points(x1,dat[dat[,1]==temprs[j],i+1],type="l",col=col1[i],lwd = 1)
		}
	}

	legend(lgx,lgy,temprs[n3],col =col.leg[n3],lty = 1,lwd=2.5,cex=0.65)
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


