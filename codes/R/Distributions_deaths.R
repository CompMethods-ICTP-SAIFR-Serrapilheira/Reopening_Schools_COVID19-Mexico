#Aqu? ya se considera vacunaci?n a >18

library(dplyr) #(sample_n) Sample from dataframe
library(fitdistrplus) #fit.dist
source("codes/R/fun_fitdist_deaths.R")


start.date = as.Date("2021-03-01")
end.date   = as.Date("2021-07-18")
database ="210801COVID19MEXICO.csv"
esc="universidad"
min_alumn=18
max_alumn=30
min_prof=31
max_prof=59
prop.alumnos=0.8
prop.maestros=0.2


if(!dir.exists("figs/Distribution_adjustment")){
  dir.create("figs/Distribution_adjustment")
}

if(!dir.exists(paste("figs/Distribution_adjustment/",format(end.date+14, "%Y-%m-%d"),"_",esc, sep=""))){
  dir.create(paste("figs/Distribution_adjustment/",format(end.date+14, "%Y-%m-%d"),"_",esc, sep=""))
}

if(!dir.exists(paste("figs/Distribution_adjustment/",format(end.date+14, "%Y-%m-%d"),"_",esc,"/Deaths", sep=""))){
  dir.create(paste("figs/Distribution_adjustment/",format(end.date+14, "%Y-%m-%d"),"_",esc,"/Deaths", sep=""))
}


R = read.csv(paste("data/",database, sep=""), header=T) # raw COVID-19 data
Pos = R[R$CLASIFICACION_FINAL %in% c(1, 2, 3),]

alumnos=Pos[Pos$EDAD>=min_alumn& Pos$EDAD<=max_alumn,]
maestros=Pos[Pos$EDAD>min_prof & Pos$EDAD<max_prof,]

alumnos_deaths=alumnos[alumnos$FECHA_DEF!="9999-99-99",]
maestros_deaths=maestros[maestros$FECHA_DEF!="9999-99-99",]


#----Execution----
pop=c(300) #,500,800,1000)
vac=c("SI") #, "NO")

for (v in vac){
  for(k in pop){
    fit.dist(k, v, paste(format(end.date+14, "%Y-%m-%d"),"_",esc, sep="")) #"SI" o "NO" se refiere a vacunaci?n
  }
}

