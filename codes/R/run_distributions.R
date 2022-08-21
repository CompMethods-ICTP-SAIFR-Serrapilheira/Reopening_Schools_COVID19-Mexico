#' --------------------------------------
#' title: run_distributions.R
#' author: Ruth Corona-Moreno
#' description: 
#' The script performed the adjustment of beta distribution probability function for severe cases and deaths.
#' In addition, a sample of the probabilities Ps, Pr, Pr2, Pr3, sigma, Psev and Pd are performed 
#'
#' input parameters:
#' --------------------------------------
#' --database[string]          Database name according to official website https://www.gob.mx/salud/documentos/datos-abiertos-bases-historicas-direccion-general-de-epidemiologia 
#' --start.date[string]        First date considered in data set to adjust the probability distribution (Date string "Y-m-d")
#' --end.date[string]          Last date considered in data set to adjust the probability distribution (Date string "Y-m-d")
#' --esc[string]               Name assigned to save simulations
#' --min_stud[int]             Minimum age for students
#' --max_stud[int]     		     Maximum age for students
#' --min_prof[int]             Minimum age for teachers
#' --max_prof[int]             Maximum age for teachers
#' --prop.students[float]      Proportion of students in each sample
#' --vac[vect]                 c("YES", "NO" )
#'
#'
#' Output
#' --------------------------------------
#' 1) File: paste("output/DistProb_Vac",vac,".csv",sep=""))
#' 
#' Dependencies from base R:
#' --------------------------------------
#' 1) dplyr
#' 2) fitdistrplus
#' 
#' Databases required:
#' --------------------------------------
#' data/210801COVID19MEXICO.csv   Downloaded from https://www.gob.mx/salud/documentos/datos-abiertos-bases-historicas-direccion-general-de-epidemiologia
#'
#' Dependencies from codes in this repository:
#' --------------------------------------
#' 1) codes/R/functions/fun_fitdist.R
#'


library(dplyr) #required for "sample_n" function
library(fitdistrplus) #required for "fitdist" function
source("codes/R/functions/fun_fitdist.R") 

#----------Parameters----------
database ="210801COVID19MEXICO.csv"

#Time window to sample the population
start.date = as.Date("2021-03-01") 
end.date   = as.Date("2021-07-18")
n=abs(as.numeric(as.Date(start.date, format="%Y-%m-%d")-as.Date(end.date, format="%Y-%m-%d"), unit="days"))

esc="university" 
min_stud=18 
max_stud=30 
min_prof=31 
max_prof=59 
prop.students=0.8 
prop.teachers=1-prop.students 
pop=c(1000) 
vac=c("YES") 

#----------Creating directory to save figures----------
if(!dir.exists("figs/Distribution_adjustment")){
  dir.create("figs/Distribution_adjustment")
}

if(!dir.exists(paste("figs/Distribution_adjustment/",format(end.date+14, "%Y-%m-%d"),"_",esc, sep=""))){
  dir.create(paste("figs/Distribution_adjustment/",format(end.date+14, "%Y-%m-%d"),"_",esc, sep=""))
}

if(!dir.exists(paste("figs/Distribution_adjustment/",format(end.date+14, "%Y-%m-%d"),"_",esc,"/Deaths", sep=""))){
  dir.create(paste("figs/Distribution_adjustment/",format(end.date+14, "%Y-%m-%d"),"_",esc,"/Deaths", sep=""))
}

if(!dir.exists(paste("figs/Distribution_adjustment/",format(end.date+14, "%Y-%m-%d"),"_",esc,"/Severe", sep=""))){
  dir.create(paste("figs/Distribution_adjustment/",format(end.date+14, "%Y-%m-%d"),"_",esc,"/Severe", sep=""))
}

#---------Importing data----------
R = read.csv(paste("data/",database, sep=""), header=T) # raw COVID-19 data

#Filtering positive cases
Pos = R[R$CLASIFICACION_FINAL %in% c(1, 2, 3),]

#Filtering students and teachers data according to the age.
students=Pos[Pos$EDAD>=min_stud& Pos$EDAD<=max_stud,]
teachers=Pos[Pos$EDAD>min_prof & Pos$EDAD<max_prof,]

students_deaths=students[students$FECHA_DEF!="9999-99-99",]
teachers_deaths=teachers[teachers$FECHA_DEF!="9999-99-99",]

#----Execution----
#Density probability adjustment for each population size (pop) and vaccination status (vac)
for (v in vac){
  for(k in pop){
    #Run "fit.dist.deaths" to adjust a beta distribution to deaths cases
    adjust.deaths=fit.dist.deaths(k, v, paste(format(end.date+14, "%Y-%m-%d"),"_",esc, sep="")) 
    #First parameter (a) of beta distribution adjusted
    a_death=adjust.deaths[[1]]
    #Second parameter (b) of beta distribution adjusted
    b_death=adjust.deaths[[2]]
    
    #Run "fit.dist.severe" to adjust a beta distribution to severe cases
    adjust.severe=fit.dist.severe(k, v, paste(format(end.date+14, "%Y-%m-%d"),"_",esc, sep="")) 
    #First parameter (a) of beta distribution adjusted
    a_severe=adjust.severe[[1]]
    #Second parameter (b) of beta distribution adjusted
    b_severe=adjust.severe[[2]]
    
    
    #Sampling probability distributions
    x=distributions(n, v, a_death, b_death, a_severe, b_severe)  
    print(head(x))
    write.csv(x, paste("output/DistProb_Vac",vac,".csv",sep=""))
  }
}

