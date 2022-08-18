library(dplyr) # >%>
library(interp)
library(tidyr) #expand function
library(plotly)
source("codes/R/fun_simulations.R", local = TRUE)

phi <- 0.3
Sp <- 0.97
Se <- 0.80
TP0inic <-  0 
FP0inic <- 0
estimations_date <- '2021-09-26'
TimeWindow <- 3
pop_tot <- 1.05*(10^6)
tau.ext <- 0
vac <- 'YES'
lugar.name <- 'QUERETARO_CAP'
esc.size <- 1000
D0inic <- 0
S0inic <- 0 
Rec0inic <- 8
b.esc <- 0.1
w <- 0
type.test <- 3
tau <- 0.05
ci <- 3


#########
#tau_Func
TauMax=esc.size/2

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

df_simulations <- data.frame(x$fechas_simul,x$x.simul.Simul.U, x$simul.E, x$simul.A, x$simul.TP, x$simul.FP, x$simul.S, x$simul.R, x$simul.D, x$simul.beta.v, x$simul.rt_esc, x$simul.tests)
write.csv(df_simulations, paste("output/",lugar.name,"_",inic_esc,"_",fin_esc,"_Vac",vac,"_bEsc",b.esc,"_Tot",esc.size,"_w",w,"_Test",type.test,"_tau",tau,"_A0E0",ci,".csv", sep=""))


############   Plot   ##############################
plot_func(lugar.name, inic_esc, fin_esc, vac, b.esc, esc.size, w, type.test, tau, ci)
