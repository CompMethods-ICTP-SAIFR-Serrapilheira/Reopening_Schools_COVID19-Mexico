
rm(list=ls())

library(openxlsx)
library(stringi) #Package to remove accents
library(fitdistrplus)
library(xtable)
library(plotly)
library(EnvStats) #For pdfPlot function
library(rsample)
library(dplyr) #Sample from dataframe
library(boot)
set.seed(3)

start.date = as.Date("2021-03-01")
end.date   = as.Date("2021-07-18")
database ="210801COVID19MEXICO.csv"

setwd("C:/Users/ruthc/Google Drive/Tesis/CODIGOS_FINALES_Julio8_ComparModels/Datos")
R = read.csv(paste("./",database, sep=""), header=T) # raw COVID-19 data
Pos = R[R$CLASIFICACION_FINAL %in% c(1, 2, 3),]

esc="universidad"
min_alumn=18
max_alumn=30
min_prof=31
max_prof=59

alumnos=Pos[Pos$EDAD>=min_alumn& Pos$EDAD<=max_alumn,]
maestros=Pos[Pos$EDAD>min_prof & Pos$EDAD<max_prof,]

#dim(alumnos) #475036 casos positivos en todo el país
#dim(maestros)# 1242936 casoso positivos en todo el país

setwd("C:/Users/ruthc/Google Drive/Tesis/CODIGOS_FINALES_Julio8_ComparModels")
dir.create(paste("Distribution_adjustment/Distributions_",format(end.date+14, "%Y-%m-%d"),"_",esc, sep=""))
setwd(paste("C:/Users/ruthc/Google Drive/Tesis/CODIGOS_FINALES_Julio8_ComparModels/Distribution_adjustment/Distributions_",format(end.date+14, "%Y-%m-%d"),"_",esc, sep=""))



#----Function to fit a distribution----
fit.dist<-function(n, vac){
    pop.tot=n
    
    if(vac=="SI"){
      rho=0.9 #rho: porcentaje de eficacia de vacuna cansino de no presentar enfermedad grave.
    }else{
      rho=0
    }
    
    prop.alumnos=0.8
    prop.maestros=0.2
    print(paste("Población total= ",pop.tot, " Alumnos=",pop.tot*prop.alumnos ))
    
    boot.size=1000
    prob.sev=c(NA, numeric(boot.size))
    
    for (i in 1:boot.size){
      x=sample_n(alumnos, size=pop.tot*prop.alumnos, replace=F)
      y=sample_n(maestros, size=pop.tot*prop.maestros, replace=F)
      y.sever=y[y$INTUBADO == 1 | y$FECHA_DEF!="9999-99-99",]
      x.sever=x[x$INTUBADO == 1 | x$FECHA_DEF!="9999-99-99",]
      muestra<- rbind(x, y)
      Sev=muestra[muestra$INTUBADO == 1 | muestra$FECHA_DEF!="9999-99-99",]
      p=(dim(Sev)[1]-rho*(dim(y.sever)[1]+dim(x.sever)[1]))/pop.tot
      
      prob.sev[i]=p
    }
    
    #---Fit distribution----
    x=prob.sev[prob.sev>0]
    x=x[!is.na(x)]
    boxplot(x, main='Probabilidad')
    cuantil=0.98
    quant.p= quantile(x, cuantil)
    
    fit.beta     = fitdist(x[x<=quant.p], distr = "beta") 
    Dist.beta    = summary(fit.beta)
    a            = Dist.beta$estimate["shape1"]
    b            = Dist.beta$estimate["shape2"]
    mean.p         = a/(a + b)
    plot(fit.beta) 
    
    setwd(paste("C:/Users/ruthc/Google Drive/Tesis/CODIGOS_FINALES_Julio8_ComparModels/Distribution_adjustment/Distributions_",format(end.date+14, "%Y-%m-%d"),"_",esc, sep=""))
    if(vac=="SI"){
      pdf(paste("./VacSev_popTot",pop.tot,"_al",prop.maestros, "_alAge",min_alumn,"-",max_alumn, "_prfAge",min_prof,"-",max_prof,"_q",cuantil,"_B",sprintf(a, fmt = '%#.4f'),"_",sprintf(b, fmt = '%#.4f'),".pdf", sep=""))
    }else{
      pdf(paste("./Sev_popTot",pop.tot,"_al",prop.maestros, "_alAge",min_alumn,"-",max_alumn, "_prfAge",min_prof,"-",max_prof,"_q",cuantil,"_B",sprintf(a, fmt = '%#.4f'),"_",sprintf(b, fmt = '%#.4f'),".pdf", sep=""))
    }
    plot(fit.beta) 
    dev.off()
    
    print(c( paste("p.sev.school", pop.tot,"_al",prop.alumnos,"_alAge",min_alumn,"-",max_alumn, "_prfAge",min_prof,"-",max_prof, quant.p, a, b, mean.p)))
}


#----Execution----
pop=c(300,500,800,1000)
vacuna=c("SI", "NO")
for(k in pop){
  for(v in vacuna){
    fit.dist(k, v)
  }
}


