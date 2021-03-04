histlife<-function(estadios,fenologia){
  inmaduros <-  estadios[-(length(estadios)-1):-(length(estadios))]
  maduros   <-  estadios[(length(estadios)-1):(length(estadios))]
  finma_dtdv<-finma_spdv<-finma_fdv<-finma_fpdv<-finma_mort<-finma_pmort<-list(1)
  finma_sedt<-finma_sedr<-finma_sem<-list(1)

  # Inmaduros
    for(j in 1:length(inmaduros))
      {
        finma_dtdv[[j]]   <-  fenologia[[1]][[j]]
        finma_spdv[[j]]   <-  fenologia[[2]][[j]]
	finma_sedt[[j]]	  <-  fenologia[[3]][[j]]

        finma_fdv[[j]]    <-  fenologia[[4]][[j]]
        finma_fpdv[[j]]   <-  fenologia[[5]][[j]]
	finma_sedr[[j]]	  <-  fenologia[[6]][[j]]

        finma_mort[[j]]   <-  fenologia[[7]][[j]]
        finma_pmort[[j]]  <-  fenologia[[8]][[j]]
	finma_sem[[j]]	  <-  fenologia[[9]][[j]]
      }
  # Maduros

    fmadu_hfdv<-fmadu_hspdv<-fmadu_hff<-fmadu_hspf<-fmadu_htz<-fmadu_hovi<-fmadu_htzp<-fmadu_hovip<-list(1)
    fmadu_hsedt<-fmadu_hsedr<-fmadu_hseto<-fmadu_hsero<-list(1)   

    # Hembras

        fmadu_hfdv[[1]]   <-  fenologia[[1]][[j+1]]
        fmadu_hspdv[[1]]  <-  fenologia[[2]][[j+1]]
	fmadu_hsedt[[1]]  <-  fenologia[[3]][[j+1]]

        fmadu_hff[[1]]    <-  fenologia[[4]][[j+1]]
        fmadu_hspf[[1]]   <-  fenologia[[5]][[j+1]]
	fmadu_hsedr[[1]]  <-  fenologia[[6]][[j+1]]

        fmadu_htz[[1]]    <-  fenologia[[10]][[1]]
        fmadu_htzp[[1]]   <-  fenologia[[11]][[1]]
	fmadu_hseto[[1]]  <-  fenologia[[12]][[1]]

        fmadu_hovi[[1]]   <-  fenologia[[13]][[1]]
        fmadu_hovip[[1]]  <-  fenologia[[14]][[1]]
	fmadu_hsero[[1]]  <-  fenologia[[15]][[1]]

    fmadu_mfdv<-fmadu_mspdv<-fmadu_mff<-fmadu_mspf<-list(1)
    fmadu_msedt<-fmadu_msedr<-list(1)

    # Machos

    if(length(fenologia[[1]]) +1 == length(estadios))
    {
        fmadu_mfdv[[1]]   <-  "logit" ## debe tener un valor muy alto para que el tiempo medio sea pequeñisimo
        fmadu_mspdv[[1]]  <-  100000000  ## Cosiderando este modelo se acomoda a nuestra condicion
	fmadu_msedt[[1]]  <-  c(slope=7777777)
        fmadu_mff[[1]]    <-  y~1+a*1e-012*x   ## obligamos que la tasa sea mayor a 1, lo cual indicaria que viviría menos de 1 dia
        fmadu_mspf[[1]]   <-  c(a=1)
	fmadu_msedr[[1]]  <-  c(a=0.01)
    }else
    {
        fmadu_mfdv[[1]]   <-  fenologia[[1]][[j+2]]
        fmadu_mspdv[[1]]  <-  fenologia[[2]][[j+2]]
	fmadu_msedt[[1]]  <-  fenologia[[3]][[j+2]]
        fmadu_mff[[1]]    <-  fenologia[[4]][[j+2]]
        fmadu_mspf[[1]]   <-  fenologia[[5]][[j+2]]
	fmadu_msedr[[1]]  <-  fenologia[[6]][[j+2]]
    }

 salidas<-list(distri_dv=finma_dtdv,slope_dv=finma_spdv,StEr_DTime=finma_sedt,fdv_dv=finma_fdv,pdv_dv=finma_fpdv,StEr_DRate=finma_sedr,mortal=finma_mort,pmortal=finma_pmort,StEr_Mort=finma_sem,
                distri_snh=fmadu_hfdv[[1]],slope_snh=fmadu_hspdv[[1]],StEr_DTimeH=fmadu_hsedt[[1]],fh_h=fmadu_hff[[1]],pfh_h=fmadu_hspf[[1]],StEr_DRateH=fmadu_hsedr[[1]],
		ftazaeh_h=fmadu_htz[[1]],ptazaeh_h=fmadu_htzp[[1]],StEr_ToviH=fmadu_hseto[[1]],fovih_h=fmadu_hovi[[1]],povih_h=fmadu_hovip[[1]],StEr_RoviH=fmadu_hsero[[1]],
                distri_snm=fmadu_mfdv[[1]],slope_snm=fmadu_mspdv[[1]],StEr_DTimeM=fmadu_msedt[[1]],fm_m=fmadu_mff[[1]],pfm_m=fmadu_mspf[[1]],StEr_DRateM=fmadu_msedr[[1]])
  return(salidas)
}
###########################################################
histlife_ant<-function(estadios,fenologia)
{
  inmaduros <-  estadios[-(length(estadios)-1):-(length(estadios))]
  maduros   <-  estadios[(length(estadios)-1):(length(estadios))]
  finma_dtdv<-finma_spdv<-finma_fdv<-finma_fpdv<-finma_mort<-finma_pmort<-list(1)
  # Inmaduros
    for(j in 1:length(inmaduros))
      {
        finma_dtdv[[j]]   <-  fenologia[[1]][[j]]
        finma_spdv[[j]]   <-  fenologia[[2]][[j]]
        finma_fdv[[j]]    <-  fenologia[[3]][[j]]
        finma_fpdv[[j]]   <-  fenologia[[4]][[j]]
        finma_mort[[j]]   <-  fenologia[[5]][[j]]
        finma_pmort[[j]]  <-  fenologia[[6]][[j]]
      }
  # Maduros
    fmadu_hfdv<-fmadu_hspdv<-fmadu_hff<-fmadu_hspf<-fmadu_htz<-fmadu_hovi<-fmadu_htzp<-fmadu_hovip<-list(1)
    # Hembras
        fmadu_hfdv[[1]]   <-  fenologia[[1]][[j+1]]
        fmadu_hspdv[[1]]  <-  fenologia[[2]][[j+1]]
        fmadu_htz[[1]]    <-  fenologia[[7]][[1]]
        fmadu_hovi[[1]]   <-  fenologia[[9]][[1]]
        fmadu_hovip[[1]]  <-  fenologia[[10]][[1]]
        fmadu_htzp[[1]]   <-  fenologia[[8]][[1]]
        fmadu_hff[[1]]    <-  fenologia[[3]][[j+1]]
        fmadu_hspf[[1]]   <-  fenologia[[4]][[j+1]]
    fmadu_mfdv<-fmadu_mspdv<-fmadu_mff<-fmadu_mspf<-list(1)
    # Machos
    if(length(fenologia[[1]]) +1 == length(estadios))
    {
        fmadu_mfdv[[1]]   <-  "logit" ## debe tener un valor muy alto para que el tiempo medio sea pequeñisimo
        fmadu_mspdv[[1]]  <-  100000000  ## Cosiderando este modelo se acomoda a nuestra condicion
        fmadu_mff[[1]]    <-  y~1+a*1e-012*x   ## obligamos que la tasa sea mayor a 1, lo cual indicaria que viviría menos de 1 dia
        fmadu_mspf[[1]]   <-  c(a=1)
    }else
    {
        fmadu_mfdv[[1]]   <-  fenologia[[1]][[j+2]]
        fmadu_mspdv[[1]]  <-  fenologia[[2]][[j+2]]
        fmadu_mff[[1]]    <-  fenologia[[3]][[j+2]]
        fmadu_mspf[[1]]   <-  fenologia[[4]][[j+2]]
    }

 salidas<-list(distri_dv=finma_dtdv,slope_dv=finma_spdv,fdv_dv=finma_fdv,pdv_dv=finma_fpdv,moratl=finma_mort,pmortal=finma_pmort,
                distri_snh=fmadu_hfdv[[1]],slope_snh=fmadu_hspdv[[1]],fh_h=fmadu_hff[[1]],pfh_h=fmadu_hspf[[1]],ftazaeh_h=fmadu_htz[[1]],ptazaeh_h=fmadu_htzp[[1]],fovih_h=fmadu_hovi[[1]],povih_h=fmadu_hovip[[1]],
                distri_snm=fmadu_mfdv[[1]],slope_snm=fmadu_mspdv[[1]],fm_m=fmadu_mff[[1]],pfm_m=fmadu_mspf[[1]])
  return(salidas)
}
#############################################################################################################
#############################################################################################################
tosimulateinsim<-function(modelim,modelm,fenologia, estadios, hfeno, xi){
return(list(modelim=modelim,modelm=modelm,fenologia=fenologia, estadios=estadios, hfeno=hfeno, xi=xi))
}
#############################################################################################################
#############################################################################################################
tosimulateinsimVariable<-function(modelim,modelm,fenologia, estadios, hfeno, xi, ff, paramR){
	return(list(modelim=modelim,modelm=modelm,fenologia=fenologia, estadios=estadios, hfeno=hfeno, xi=xi, ff=ff, paramR=paramR))
}
#############################################################################################################
#############################################################################################################
tosimulate<-function(modelim,modelm,GImat,AImat,fenologia, estadios, hfeno, xi){
return(list(modelim=modelim,modelm=modelm,Gimat=GImat,Aimat=AImat,fenologia=fenologia, estadios=estadios, hfeno=hfeno, xi=xi))
}
#############################################################################################################
#############################################################################################################
indicesNew<-function(method,modelim,modelm,estadios,Index)
{
Tmin <- seq(-57, 33,0.5)
Tmax <- seq(-51, 48, 0.5)

gis2<-gis1<-matrix(0,nrow=length(Tmin), ncol=length(Tmax))
colnames(gis1)<-Tmax
rownames(gis1)<-Tmin 
colnames(gis2)<-Tmax
rownames(gis2)<-Tmin
if(Index=="GI") for(i in 1:nrow(gis1)) for(j in 1:ncol(gis1)) if(Tmin[i] < Tmax[j]) gis1[i,j]<-GenActIndex("Dete",data.frame(Tmin[i],Tmax[j]),modelim,modelm,estadios)$GI
if(Index=="AI") for(i in 1:nrow(gis2)) for(j in 1:ncol(gis2)) if(Tmin[i] < Tmax[j]) gis2[i,j]<-GenActIndex("Dete",data.frame(Tmin[i],Tmax[j]),modelim,modelm,estadios)$AI

if(Index=="GI") return(list(GI_gis=gis1))
if(Index=="AI")return(list(AI_gis=gis2))
}