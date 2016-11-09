#************************************************************************
#************************************************************************
#Calculo del área dentro de la línea crítica cuando conocemos
#los puntos (x,y) de la línea crítica
#************************************************************************
#************************************************************************

rm(list=ls())

#************************************************************************
#Uso
#************************************************************************

# >Rscript file_input
# >Rscript "ce.dat"

#Read command arguments

args<-commandArgs(TRUE)
arg_input<-(args[1])

if(file.exists(arg_input)){
    print("Input file exists")
}else{
    print("ERROR: Input file not found")
    q(save='no')
}

#************************************************************************
#Es necesario cargar las siguientes librerías en R
#Quitar la # para instalar librerías (sólo en la primera ocasión)
#************************************************************************

#install.packages('magic')
#install.packages('abind')
#install.packages('geometry')
#install.packages('plotrix')

library('magic')
library('abind')
library('geometry')

#Leo el archivo de entrada
#Ordeno los puntos para generar un polígono cerrado
#------------------------------------------------------------------------

ce_path<-arg_input
ce<-read.table(paste(ce_path))

index<-c(1:length(ce[,2]))
flag<-c(1:length(ce[,2]))*0

orden<-c(1:500)
counter<-1

ind<-1
flag[ind]<-1
orden[counter]<-ind

while(length(which(flag==0))>0){

    dist<-sqrt( (ce[,2]-ce[ind,2])^2 + (ce[,3]-ce[ind,3])^2 )
    dist[which(flag==1)]<-0
    distnonzero<-dist[which(dist>0)]
    ind_nonzero<-which(dist==min(distnonzero))

    ind_included<-c(1:length(ind_nonzero))*0

    for (rr in 1:length(ind_nonzero)){
        
        if(flag[ind_nonzero[rr]]==0){
            ind_included[rr]<-1
        }
        
    }

    if(length(which(ind_included==1))>1){
        ind<-min(ind_nonzero)
    }
    
    if(length(which(ind_included==1))==1){
        ind<-ind_nonzero[which(ind_included==1)]
    }

    flag[ind_nonzero]<-1
    #print(ind)
    counter<-counter+1
    #print(counter)
    orden[counter]<-ind

}

orden<-orden[1:counter]
orden<-c(orden,orden[1])

#------------------------------------------------------------------------
#Área dentro de la línea crítica
#------------------------------------------------------------------------

area<-polyarea(ce[orden,2],ce[orden,3])
print(paste('AREA = ',area,sep=""))

#------------------------------------------------------------------------
#Gráficos
#------------------------------------------------------------------------

#library('plotrix')

#Para ce.dat
#plot(ce[,2],ce[,3],type='b',xlim=c(-16,16),ylim=c(-16,16),xlab='x',ylab='y')
#lines(ce[orden,2],ce[orden,3],col='red',pch=20,type='b')
#draw.circle(0,0,sqrt(area/pi),border='blue')

#Para ci.dat
#plot(ce[,2],ce[,3],type='b',xlim=c(-4,4),ylim=c(-4,4),xlab='x',ylab='y')
#lines(ce[orden,2],ce[orden,3],col='red',pch=20,type='b')
#draw.circle(0,0,sqrt(area/pi),border='blue')

q(save='no')
#-------------------------------------------------------------------------
