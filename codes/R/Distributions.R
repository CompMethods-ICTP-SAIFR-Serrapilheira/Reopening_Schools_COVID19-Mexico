library(dplyr) #required for "sample_n" function
library(fitdistrplus) #required for "fitdist" function
source("codes/R/fun_fitdist.R") #Import "fit.dist" function

#----------Parameters----------
database ="210801COVID19MEXICO.csv"

#Time window to sample the population
start.date = as.Date("2021-03-01") 
end.date   = as.Date("2021-07-18")
n=abs(as.numeric(as.Date(start.date, format="%Y-%m-%d")-as.Date(end.date, format="%Y-%m-%d"), unit="days"))

esc="university" #Name of the school
min_stud=18 #Minimum age for students
max_stud=30 #Maximum age for students
min_prof=31 #Minimum age for teachers
max_prof=59 #Maximum age for teachers
prop.students=0.8 #Proportion of students in each sample
prop.teachers=1-prop.students #Proportion of teacher in each sample
pop=c(1000) #population in school (size sample)
vac=c("YES") #, "NO") #vaccination status for population in the school

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
    adjust.deaths=fit.dist.deaths(k, v, paste(format(end.date+14, "%Y-%m-%d"),"_",esc, sep="")) 
    a_death=adjust.deaths[[1]]
    b_death=adjust.deaths[[2]]
    
    adjust.severe=fit.dist.severe(k, v, paste(format(end.date+14, "%Y-%m-%d"),"_",esc, sep="")) 
    a_severe=adjust.severe[[1]]
    b_severe=adjust.severe[[2]]
    
    x=distributions(n, v, a_death, b_death, a_severe, b_severe)  
    print(head(x))
    write.csv(x, paste("data/DistProb_","_Vac",vac,".csv",sep=""))
  }
}

