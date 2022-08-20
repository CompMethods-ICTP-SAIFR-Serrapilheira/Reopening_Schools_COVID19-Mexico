#Purpose:Simulate data of COVID-19 epidemic dynamic in a university in Mexico
#Author: Ruth Corona-Moreno

library(dplyr) 
library(interp)
library(tidyr) 
library(plotly)
source("codes/R/functions/fun_simulations.R", local = TRUE)


#*************************************************************************
#------------------------Set parameters values----------------------------

#Location name and total population
lugar.name <- 'QUERETARO_CAP'
pop_tot <- 1.05*(10^6)

#School population initial conditions
vac <- 'YES'
esc.size <- 1000
TP0inic <-  0 
FP0inic <- 0
D0inic <- 0
S0inic <- 0 
Rec0inic <- 8
ci <- 3 #E0inic=A0inic=ci

#Estimations from covidestim
estimations_date <- '2021-09-26'
TimeWindow <- 3

#Vaccination conditions
type.test <- 3
tau <- 0.05
tau.ext <- 0
Sp <- 0.97
Se <- 0.80

#Parameters of the model
b.esc <- 0.1
w <- 0
phi <- 0.3
TauMax=esc.size/2 #Required to define the linear screening strategy


#*************************************************************************
#------------------------Time windows definition---------------------------
if(TimeWindow==1){
  inic_esc='2020-08-02'
  fin_esc= '2020-12-17'
}else if(TimeWindow==2){
  inic_esc='2021-03-05'
  fin_esc= '2021-07-20'
}else if(TimeWindow==3){
  inic_esc='2021-05-01'
  fin_esc= '2021-09-15'
}else{
  print("Error")
}


#*************************************************************************
#--------------------------Making simulations-----------------------------
x <- modelo_func(
  vacuna=vac,
  b.esc=b.esc,
  lugar=lugar.name,
  esc.size =esc.size,
  w=w,
  type.test=type.test,
  tau=tau,
  A0inic=ci, 
  E0inic=ci 
)

#Saving simulations
df_simulations <- data.frame(x$fechas_simul,x$x.simul.Simul.U, x$simul.E, x$simul.A, x$simul.TP, x$simul.FP, x$simul.S, x$simul.R, x$simul.D, x$simul.beta.v, x$simul.rt_esc, x$simul.tests)
write.csv(df_simulations, paste("output/",lugar.name,"_",inic_esc,"_",fin_esc,"_Vac",vac,"_bEsc",b.esc,"_Tot",esc.size,"_w",w,"_Test",type.test,"_tau",tau,"_A0E0",ci,".csv", sep=""))


#Plot simulations
plot_func(lugar.name, inic_esc, fin_esc, vac, b.esc, esc.size, w, type.test, tau, ci)
