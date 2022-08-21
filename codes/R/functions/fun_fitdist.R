#' --------------------------------------
#' title: fun_fitdist.R
#' author: Ruth Corona-Moreno
#' description: The script has the functions "fit.dist.deaths", "fit.dist.severe" and "distributions" called by the "run_distributions.R" code.
#' 
#' "fit.dist.deaths"
#' --------------
#' This function perform the estimation of the density probability distribution of deaths 
#' 
#' "fit.dist.severe"
#' --------------
#' This function perform the estimation of the density probability distribution of deaths 
#' 
#' "distributions"
#' --------------
#' This function perform the sampling of all density probability distributions, which are required to perform the simulations of the model. 



#----Function to fit a distribution-----
fit.dist.deaths<-function(n, vac, directory){
  pop.tot=n

  if(vac=="YES"){
    rho=0.9 #rho: porcentaje de eficacia de vacuna cansino de no presentar enfermedad grave.
  }else{
    rho=0
  }

  boot.size=1000
  prob.die=c(NA, numeric(boot.size))

  for (i in 1:boot.size){
    x=sample_n(students, size=pop.tot*prop.students, replace=F)
    y=sample_n(teachers, size=pop.tot*prop.teachers, replace=F)
    y.deaths=dim(y)[1]-table(y$FECHA_DEF)[length(table(y$FECHA_DEF))][[1]]
    x.deaths=dim(x)[1]-table(x$FECHA_DEF)[length(table(x$FECHA_DEF))][[1]]
    muestra<- rbind(x, y)
    deaths=muestra[muestra$FECHA_DEF!="9999-99-99",]
    p=(dim(deaths)[1]-rho*(y.deaths+x.deaths))/pop.tot

    prob.die[i]=p
  }

  #---Fit distribution----
  x=prob.die[prob.die>0]
  x=x[!is.na(x)]
  boxplot(x, main='Probability')
  cuantil=0.98
  quant.p= quantile(x, cuantil)

  fit.beta     = fitdist(x[x<=quant.p], distr = "beta")
  Dist.beta    = summary(fit.beta)
  a            = Dist.beta$estimate["shape1"]
  b            = Dist.beta$estimate["shape2"]
  mean.p         = a/(a + b)
  plot(fit.beta)


  if(vac=="YES"){
    pdf(paste("figs/Distribution_adjustment/",directory,"/Deaths/VacDth_popTot",pop.tot,"_al",prop.teachers, "_alAge",min_stud,"-",max_stud, "_prfAge",min_prof,"-",max_prof,"_q",cuantil,"_B",sprintf(a, fmt = '%#.4f'),"_",sprintf(b, fmt = '%#.4f'),".pdf", sep=""))
  }else{
    pdf(paste("figs/Distribution_adjustment/",directory,"/Deaths/Dth_popTot",pop.tot,"_al",prop.teachers, "_alAge",min_stud,"-",max_stud, "_prfAge",min_prof,"-",max_prof,"_q",cuantil,"_B",sprintf(a, fmt = '%#.4f'),"_",sprintf(b, fmt = '%#.4f'),".pdf", sep=""))
  }
  plot(fit.beta)
  dev.off()

  return(list(a,b))
}


fit.dist.severe<-function(n, vac, directory){
  pop.tot=n
  
  if(vac=="YES"){
    rho=0.9 #rho: porcentaje de eficacia de vacuna cansino de no presentar enfermedad grave.
  }else{
    rho=0
  }
  
  boot.size=1000
  prob.sev=c(NA, numeric(boot.size))
  
  for (i in 1:boot.size){
    x=sample_n(students, size=pop.tot*prop.students, replace=F)
    y=sample_n(teachers, size=pop.tot*prop.teachers, replace=F)
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
  boxplot(x, main='Probability')
  cuantil=0.98
  quant.p= quantile(x, cuantil)
  
  fit.beta     = fitdist(x[x<=quant.p], distr = "beta") 
  Dist.beta    = summary(fit.beta)
  a            = Dist.beta$estimate["shape1"]
  b            = Dist.beta$estimate["shape2"]
  mean.p         = a/(a + b)
  plot(fit.beta) 
  

  if(vac=="YES"){
    pdf(paste("figs/Distribution_adjustment/",directory,"/Severe/VacSev_popTot",pop.tot,"_al",prop.teachers, "_alAge",min_stud,"-",max_stud, "_prfAge",min_prof,"-",max_prof,"_q",cuantil,"_B",sprintf(a, fmt = '%#.4f'),"_",sprintf(b, fmt = '%#.4f'),".pdf", sep=""))
  }else{
    pdf(paste("figs/Distribution_adjustment/",directory,"/Severe/Sev_popTot",pop.tot,"_al",prop.teachers, "_alAge",min_stud,"-",max_stud, "_prfAge",min_prof,"-",max_prof,"_q",cuantil,"_B",sprintf(a, fmt = '%#.4f'),"_",sprintf(b, fmt = '%#.4f'),".pdf", sep=""))
  }
  plot(fit.beta) 
  dev.off()
  
  return(list(a,b))
}


distributions<- function(n, vacuna, a_death, b_death, a_severe, b_severe){
  nrows=1000
  ncols=n
  
  Ps<-NULL
  Pr<-NULL
  Pr2<-NULL
  Pr3<-NULL
  sigma<-NULL
  Psev<-NULL
  Pd<-NULL
  
  
  for (i in 1:nrows){
    Ps_sample=as.numeric(dgamma(rgamma(n, shape=3.413, rate = 0.6051),shape=3.413, rate = 0.6051 ))
    #Prob recuperarse asintom?ticos
    Pr_sample=as.numeric(dgamma(rgamma(n, shape=14, rate= 2), shape=14, rate= 2) )
    #Prob recuperarse sintom?ticos NO severos
    Pr2_sample=as.numeric(dgamma(rgamma(n, shape=8.163265, scale= 3.02575), shape=8.163265, scale= 3.02575) )#no severa
    #Prob recuperarse sintom?ticos Severos
    Pr3_sample=as.numeric(dgamma(rgamma(n, shape=8.163265, scale= 3.02575), shape=8.163265, scale= 3.02575)) #severa
    sigma_sample=as.numeric(rbeta(n, shape1 = 5.143, shape2=3.536))
    if(vacuna=="NO"){
      Pd_sample=as.numeric(rbeta(n, shape1=a_death, shape2=b_death))
      Psev_sample=as.numeric(rbeta(n, shape1=a_severe, shape2=b_severe))
    }else{
      Pd_sample=as.numeric(rbeta(n, shape1=a_death, shape2=b_death))
      Psev_sample=as.numeric(rbeta(n, shape1=a_severe, shape2=b_severe))
    }
    
    
    Ps<-rbind(Ps,Ps_sample)
    Pr<-rbind(Pr,Pr_sample)
    Pr2<-rbind(Pr2,Pr2_sample)
    Pr3<-rbind(Pr3,Pr3_sample)
    sigma<-rbind(sigma,sigma_sample)
    Psev<-rbind(Psev,Psev_sample)
    Pd<-rbind(Pd,Pd_sample)
  }
  
  Ps_output<-c()
  Pr_output<-c()
  Pr2_output<-c()
  Pr3_output<-c()
  sigma_output<-c()
  Psev_output<-c()
  Pd_output<-c()
  
  for (i in 1:ncols){
    Ps_output=append(Ps_output, mean(Ps[,i]))
    Pr_output=append(Pr_output, mean(Pr[,i]))
    Pr2_output=append(Pr2_output, mean(Pr2[,i]))
    Pr3_output=append(Pr3_output, mean(Pr3[,i]))
    sigma_output=append(sigma_output, mean(sigma[,i]))
    Psev_output=append(Psev_output, mean(Psev[,i]))
    Pd_output=append(Pd_output, mean(Pd[,i]))
  }
  
  df=data.frame(Ps_output, Pr_output, Pr2_output, Pr3_output, sigma_output, Psev_output, Pd_output)
  return(df)
}