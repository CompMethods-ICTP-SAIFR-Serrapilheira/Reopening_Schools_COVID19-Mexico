#----Function to fit a distribution-----
fit.dist<-function(n, vac, directory){
  pop.tot=n

  if(vac=="SI"){
    rho=0.9 #rho: porcentaje de eficacia de vacuna cansino de no presentar enfermedad grave.
  }else{
    rho=0
  }

  ##Sample according school population
  #prop.alumnos=0.8
  #prop.maestros=0.2
  ##print(paste("Poblaci?n total= ",pop.tot, " Alumnos=",pop.tot*prop.alumnos ))

  boot.size=1000
  prob.die=c(NA, numeric(boot.size))

  for (i in 1:boot.size){
    x=sample_n(alumnos, size=pop.tot*prop.alumnos, replace=F)
    y=sample_n(maestros, size=pop.tot*prop.maestros, replace=F)
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
  boxplot(x, main='Probabilidad')
  cuantil=0.98
  quant.p= quantile(x, cuantil)

  fit.beta     = fitdist(x[x<=quant.p], distr = "beta")
  Dist.beta    = summary(fit.beta)
  a            = Dist.beta$estimate["shape1"]
  b            = Dist.beta$estimate["shape2"]
  mean.p         = a/(a + b)
  plot(fit.beta)


  if(vac=="SI"){
    pdf(paste("figs/Distribution_adjustment/",directory,"/Deaths/VacDth_popTot",pop.tot,"_al",prop.maestros, "_alAge",min_alumn,"-",max_alumn, "_prfAge",min_prof,"-",max_prof,"_q",cuantil,"_B",sprintf(a, fmt = '%#.4f'),"_",sprintf(b, fmt = '%#.4f'),".pdf", sep=""))
  }else{
    pdf(paste("figs/Distribution_adjustment/",directory,"/Deaths/Dth_popTot",pop.tot,"_al",prop.maestros, "_alAge",min_alumn,"-",max_alumn, "_prfAge",min_prof,"-",max_prof,"_q",cuantil,"_B",sprintf(a, fmt = '%#.4f'),"_",sprintf(b, fmt = '%#.4f'),".pdf", sep=""))
  }
  plot(fit.beta)
  dev.off()

}
