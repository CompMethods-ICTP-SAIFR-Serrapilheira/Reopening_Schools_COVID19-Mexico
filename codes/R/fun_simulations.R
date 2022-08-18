modelo_func<-function(
  vacuna,
  b.esc,
  lugar,
  esc.size, 
  w,
  type.test,
  tau, 
  A0inic, 
  E0inic 
){
  U0inic= esc.size-(D0inic+S0inic+Rec0inic+A0inic+E0inic+TP0inic+FP0inic)

  cltr=read.csv(paste("data/Estimations_",estimations_date,"/covidestim_",lugar,".csv", sep=""))
  cltr=transform(cltr, date= as.Date(date, origin = "1900-01-01"))
  
  idx_inic_cltr=match(as.Date(inic_esc, "%Y-%m-%d"),cltr$date)
  idx_fin_cltr=match(as.Date(fin_esc, "%Y-%m-%d"), cltr$date)
  
  #--Rt
  RT=cltr$Rt
  RT=RT[idx_inic_cltr: idx_fin_cltr]
  n=length(RT)
  
  pop_inf=cltr$pop.infectiousness
  activos<-c(numeric(length(RT)))
  
  for (t in (idx_inic_cltr:(idx_fin_cltr+1))){#(idx_fin_cltr+1))){
    activos_vect=pop_inf[(t-13):t]
    activos[t-idx_inic_cltr+1]=sum(activos_vect)
  }
  
  DistProb=read.csv(paste("output/DistProb","_Vac",vacuna,".csv", sep=""))
  generations=abs(as.numeric(as.Date(format(as.Date(inic_esc, "%Y-%m-%d"), "%Y%m%d"), format="%Y%m%d")-as.Date(format(as.Date(fin_esc, "%Y-%m-%d"), "%Y%m%d"), format="%Y%m%d"), unit="days"))
  
  #------------------------Modelo--------------------------
  U<-c(U0inic, numeric(generations-1))
  E<-c(E0inic, numeric(generations-1))
  A<-c(A0inic, numeric(generations-1))
  TP<-c(TP0inic,numeric(generations-1))
  FP<-c(FP0inic,numeric(generations-1))
  S<-c(S0inic,numeric(generations-1))
  R<-c(Rec0inic,numeric(generations-1))
  D<-c(D0inic,numeric(generations-1))
  beta.v<-c(0,numeric(generations-1))
  rt_esc<-c(0,numeric(generations-1))
  
  if(type.test==1){
    tau=esc.size*tau
    tests<-c(tau/(U0inic+A0inic+E0inic+Rec0inic),numeric(generations-1))  
  }else if(type.test==2){
    tests<-c(((2*tau*esc.size-esc.size)/esc.size)*(U0inic+E0inic+A0inic+Rec0inic)+esc.size-tau*esc.size,numeric(generations-1))
    tau=tests[1]
  }else if(type.test==3){
    tests<- c(tau, numeric(generations-1))
    tau=tests[1]
  }
  
  
  for (t in (1:(generations-1))){
    R0=RT[t]
    Ps=DistProb$Ps_output[t]
    Pr=DistProb$Pr_output[t]
    Pr2=DistProb$Pr2_output[t]
    Pr3=DistProb$Pr3_output[t]
    sigma=DistProb$sigma_output[t]
    Psev=DistProb$Psev_output[t]
    Pd=DistProb$Pd_output[t]
    
    beta=R0*(Pr*(esc.size-tau.ext)+Se*tau.ext)*(Ps*(esc.size-tau.ext)+Se*tau.ext)/(esc.size*(Pr*sigma*(esc.size-tau.ext)+Ps*(esc.size-tau.ext)*(1-sigma)+Se*tau.ext))
    beta.v[t+1]=beta
    rt_esc[t+1]=b.esc*(esc.size*(Pr*sigma*(esc.size-tau)+Ps*(esc.size-tau)*(1-sigma)+Se*tau))/((Pr*(esc.size-tau)+Se*tau)*(Ps*(esc.size-tau)+Se*tau))
    
    if(vacuna=="YES")
      alph=0.27
    else if(vacuna=="NO"){
      alph=1
    }
      
    
    x_pop=1-alph*beta*(activos[t]/pop_tot)
    mu=(1-alph*b.esc*(A[t]+E[t])/(U[t]+A[t]+E[t]+R[t]))*x_pop
    g=R[t]/(R[t]+U[t])
    
    
    if (tau/(U[t]+A[t]+E[t]+R[t])<=1 & anyNA(c(U[t],E[t], A[t]))==FALSE){
      if(type.test==1){
        tests[t]=tau/(U[t]+A[t]+E[t]+R[t])  
      }else if(type.test==2){
        tests[t+1]=trunc(((2*tau*esc.size-esc.size)/esc.size)*(U[t]+E[t]+A[t]+R[t])+esc.size-(tau*esc.size))
      }else if(type.test==3){
        if(((TP[t]+S[t])/tests[t])<0.5){
          tests[t+1]=TP[t]+S[t]+1
        }else if(((TP[t]+S[t])/tests[t])>=0.5 & ((TP[t]+S[t])/tests[t])<1){
          tests[t+1]=tests[t]
        }else if(((TP[t]+S[t])/tests[t])>1){
          tests[t+1]=TP[t]+S[t]+1
        }
      }
      
      
      U[t+1]=U[t]*(mu)*(1-(tau/(U[t]+A[t]+E[t]+R[t]))*(1-Sp))+ (1-g)*FP[t]*phi + R[t]*(1-w)
      E[t+1]=U[t]*sigma*(1-mu)+E[t]*((tau/(U[t]+A[t]+E[t]+R[t]))*(1-Se)+(1-(tau/(U[t]+A[t]+E[t]+R[t])))*(1-Ps))
      A[t+1]=U[t]*((1-sigma)*(1-mu))+A[t]*((tau/(U[t]+A[t]+E[t]+R[t]))*(1-Se)+(1-(tau/(U[t]+A[t]+E[t]+R[t])))*(1-Pr))
      TP[t+1]=E[t]*(tau/(U[t]+A[t]+E[t]+R[t]))*Se+A[t]*(tau/(U[t]+A[t]+E[t]+R[t]))*Se+TP[t]*((1-sigma)*(1-Pr)+sigma*(1-Ps))
      FP[t+1]=U[t]*(tau/(U[t]+A[t]+E[t]+R[t]))*(1-Sp)*(mu)+(1-phi)*FP[t]+R[t]*w*(tau/(U[t]+A[t]+E[t]+R[t]))*(1-Sp)
      S[t+1]=E[t]*(1-(tau/(U[t]+A[t]+E[t]+R[t])))*Ps+TP[t]*sigma*Ps+S[t]*((1-Psev)*(1-Pr2)+Psev*(1-Pd)*(1-Pr3))
      R[t+1]=TP[t]*(1-sigma)*Pr+A[t]*(1-(tau/(U[t]+A[t]+E[t]+R[t])))*Pr+S[t]*((1-Psev)*Pr2+Psev*(1-Pd)*Pr3)+R[t]*w*((tau/(U[t]+A[t]+E[t]+R[t]))*Sp+(1-(tau/(U[t]+A[t]+E[t]+R[t]))))+g*phi*FP[t] 
      D[t+1]=D[t]+S[t]*Psev*Pd
    }else{
      U[t+1]=NA
      E[t+1]=NA
      A[t+1]=NA
      TP[t+1]=NA
      FP[t+1]=NA
      S[t+1]=NA
      R[t+1]=NA
      D[t+1]=NA
      tests[t]=NA
    }
    
  }
  
  ciclo<-seq(0,generations-1)
  suma_nat=U+E+A+TP+FP+S+R+D
  simul=data.frame(ciclo, U, E, A, TP, FP, S,  R, D, suma_nat, beta.v, rt_esc, tests)
  
  
  fechas_simul=seq((as.Date("2021-08-02", "%Y-%m-%d")),(as.Date("2021-08-02", "%Y-%m-%d")+generations-1),by=1)
  df <- data.frame(fechas_simul, simul$ciclo, simul$U, simul$E, simul$A, simul$TP, simul$FP, simul$S, simul$R,  simul$D, simul$beta.v, simul$rt_esc, simul$tests)
  return(df)
}



#############################################################################################

plot_func<-function(lugar.name, inic_esc, fin_esc, vac, b.esc, esc.size, w, type.test, tau, ci){
  csv_simu=read.csv(paste("output/",lugar.name,"_",inic_esc,"_",fin_esc,"_Vac",vac,"_bEsc",b.esc,"_Tot",esc.size,"_w",w,"_Test",type.test,"_tau",tau,"_A0E0",ci,".csv",sep=""))
  
  
  plot1 <- plot_ly(x = csv_simu$x.fechas_simul[2:length((csv_simu$x.fechas_simul))], y = csv_simu$x.simul.E[2:length((csv_simu$x.fechas_simul))], type = 'scatter', mode = 'lines',name="Presymptomatic", line = list(color = "rgb(214, 39, 40)", width = 2, dash = FALSE)) 
  plot1 <- plot1%>%add_trace(x =csv_simu$x.fechas_simul[2:length((csv_simu$x.fechas_simul))], y = csv_simu$x.simul.A[2:length((csv_simu$x.fechas_simul))], type = 'scatter', mode = 'lines',name = "Aymptomatic", line = list(color = "rgb(214, 39, 40)", width = 2, dash = 'dot')
  )
  
  
  plot2 <- plot_ly(x = csv_simu$x.fechas_simul[2:length((csv_simu$x.fechas_simul))], y = csv_simu$x.simul.TP[2:length((csv_simu$x.fechas_simul))], type = 'scatter', mode = 'lines',name="True Positive", line = list(color = "rgb(252, 189, 1)", width = 2, dash = 'dash')) 
  plot2 <- plot2%>%add_trace(x =csv_simu$x.fechas_simul[2:length((csv_simu$x.fechas_simul))], y = csv_simu$x.simul.S[2:length((csv_simu$x.fechas_simul))], type = 'scatter', mode = 'lines',name = "Symptomatic", line = list(color = "rgb(252, 189, 1)", width = 2, dash = FALSE)
  )
  
  
  plot3 <- plot_ly(x = csv_simu$x.fechas_simul[2:length((csv_simu$x.fechas_simul))], y = csv_simu$x.simul.U[2:length((csv_simu$x.fechas_simul))], type = 'scatter', mode = 'lines',name="Susceptible", line = list(color = "rgb(153, 197, 49)", width = 2, dash = 'dash')) 
  plot3 <- plot3%>%add_trace(x =csv_simu$x.fechas_simul[2:length((csv_simu$x.fechas_simul))], y = (csv_simu$x.simul.U+csv_simu$x.simul.E+csv_simu$x.simul.A+csv_simu$R)[2:length((csv_simu$x.fechas_simul))], type = 'scatter', mode = 'lines',name = "In school", line = list(color = "rgb(153, 197, 49)", width = 2, dash = FALSE)
  )
  
  plot4 <- plot_ly(x = csv_simu$x.fechas_simul[2:length((csv_simu$x.fechas_simul))], y = csv_simu$x.simul.tests[2:length((csv_simu$x.fechas_simul))], type = 'scatter', mode = 'lines',name="Number of Test", line = list(color = "rgb(206, 37,201)", width = 2, dash = 'dot')) 
  
  
  mrg <- list(l = 50, r = 50,
              b = 0, t = 50,
              pad = 20)
  
  annotations = list( 
    list( 
      x = 0.5,  
      y = 1.0,  
      text = paste("Simulations for a school. Vaccination:",vac, " / Type of Test:",tau, sep=""),  font=list(size=20),  
      xref = "paper",  
      yref = "paper",  
      xanchor = "center",  
      yanchor = "bottom",  
      showarrow = FALSE
    )
  )
  
  
  fig <- subplot(plot1, plot2, plot3, plot4, nrows = 4, shareX = TRUE,  margin = c(0.05,0.05,0.05,0.05)) %>%
    layout(#plot_bgcolor='#e5ecf6', 
      xaxis = list(
        tickfont = list(size = 20),
        tickangle=-90, #Angulo de labels
        autotick=F,  #Elegir xlabels values
        dtick = 30, #Intervalo entre labels
        zerolinecolor = '#000000',#'#cccc', 
        zerolinewidth = 2, 
        gridcolor = '#cccc'), 
      yaxis = list( 
        tickfont = list(size = 19),
        zerolinecolor = '#cccc', 
        zerolinewidth = 2, 
        gridcolor = '#cccc',
        title = 'Population',
        titlefont=list(size=17)),
      yaxis2 = list( 
        tickfont = list(size = 19),
        zerolinecolor = '#000000', 
        zerolinewidth = 2, 
        gridcolor = '#cccc',
        title = 'Population',
        titlefont=list(size=17)),
      yaxis3 = list( 
        tickfont = list(size = 19),
        zerolinecolor = '#000000', 
        zerolinewidth = 2, 
        gridcolor = '#cccc',
        title = 'Population',
        titlefont=list(size=17)),
      yaxis4 = list( 
        tickfont = list(size = 19),
        zerolinecolor = '#000000', 
        zerolinewidth = 2, 
        gridcolor = '#cccc',
        title = 'Tests',
        titlefont=list(size=17)),
      margin = mrg,
      annotations = annotations,
      legend=list(font=list(size=20)))
  
  return(fig)
}


