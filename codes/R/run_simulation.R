library(dplyr) # >%>
library(interp)
library(tidyr) #expand function
source("codes/R/fun_simulations.R", local = TRUE)

phi <- 0.3
Sp <- 0.97
Se <- 0.80
TP0inic <-  0 
FP0inic <- 0
estimations_date <- '2021-09-26'
inic_esc='2021-05-01'
fin_esc= '2021-09-15'
pop_tot <- 1.05*(10^6)
tau.ext <- 0
esc <-'universidad'
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


x <- modelo_func(
  esc=esc,
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

df_simulations <- data.frame(x$fechas_simul,x$simul.U, x$simul.E, x$simul.A, x$simul.TP, x$simul.FP, x$simul.S, x$simul.R, x$simul.D, x$simul.beta.v, x$simul.rt_esc, x$simul.tests)
write.csv(df_simulations, paste("output/",lugar.name,"_",esc,"_Vac",vac,"_bEsc",b.esc,"_Tot",esc.size,"_w",w,"_Test",type.test,"_tau",tau,"_A0E0",ci,".csv", sep=""))

