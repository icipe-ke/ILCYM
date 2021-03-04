################### Control calidad datos cohorte ##############################
control.cohorte<-function(datos)
{
 for(i in 1:length(datos))
 {
  ## Problemas de observaciones 
  d1=tapply(datos[[i]][,4],datos[[i]][,1],length)
  d2= d1<3 
  #if(length(d2[d2==TRUE])>=1){cat("\n","Para la data ",i," en la temperatura de",paste(names(d2)[d2],"°C",sep="")," hay menos de 3 observaciones","\n")}
   if(length(d2[d2==TRUE])>=1){cat("\n","In temperature file ", i," at ", paste(names(d2)[d2],"°C",sep="") ," has less than 3 observations","\n")}

  ## Problemas de repeticion de dias 
  f1=function(vec){dif=length(vec)-length(unique(vec));return(dif)}
  e1=tapply(datos[[i]][,2],datos[[i]][,1],f1)
  e2= e1>0 
#  if(length(e2[e2==TRUE])>=1){cat("\n","Para la data ",i," en la temperatura de",paste(names(e2)[e2],"°C",sep="")," hay dias que se repiten","\n")}
  if(length(e2[e2==TRUE])>=1){cat("\n","In temperature file ",i," at ",paste(names(e2)[e2],"°C",sep="")," there are repeated days","\n")}

  ## Problemas de repeticion de dias 
  f2=function(vec){dif2=rep(vec[1],length(vec))==vec;return(length(dif2[dif2==FALSE]))}
  g1=tapply(datos[[i]][,3],datos[[i]][,1],f2)
  g2= g1>0 
#  if(length(g2[g2==TRUE])>=1){cat("\n","Para la data ",i," en la temperatura de",paste(names(g2)[g2],"°C",sep="")," el número de muestras no es el mismo","\n")}
  if(length(g2[g2==TRUE])>=1){cat("\n","In temperature file ",i," at ",paste(names(g2)[g2],"°C",sep="")," the number of samples are not the same","\n")}

 }
}

##################### Verifica si hay espacios en blanco #######################

control.esp<-function(datos)
{
  f3=function(vec){dif3= is.na(vec);return(length(dif3[dif3==TRUE]))}
  h1=apply(datos,2,f3)
  h2= h1>0
#  if(length(h2[h2==TRUE])>=1){cat("\n","Para la data ",i," en la temperatura de",paste(names(h2)[h2],"°C",sep="")," hay obervaciones con espacios en blanco o vacios","\n")} 
  if(length(h2[h2==TRUE])>=1){cat("\n","File ",i," in the temperature of ",paste(names(h2)[h2],"°C",sep="")," there are observations with empty spaces","\n")}
  }


######################  Retorna valores unicos en archivo ######################

findu=function(dat)
{
 n1=ncol(dat)
 k1=as.character(dat[1,1])
 for(i in 1:n1)
 {
  p1=as.character(unique(dat[,i]))
  k1=as.character(union(k1,p1))
 }

 p1=p2=k1
 p1=as.character(p1)
 p1=as.numeric(p1)
 k2<-c(1);s=1
 for(j in 1:length(p1))
 {
  pp=p1[j]+1
  if(is.na(pp)){k2[s]=p2[j];s=s+1}
 }
 return(k2)
}







































