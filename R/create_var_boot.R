create_var_boot=function(nb, dt, confounders=c(), intv=3, reNAME=NULL, x0, x1){
  var.boot=matrix(NA, nr=nb, nc=12)
  for(boot.count in 1:nb){
    dt.boot=as.data.frame(dt[sample(1:nrow(dt), replace=TRUE),])
    GT=get_theta(dt.boot,reNAME)
    if(intv==3){
      p000=omega(GT$theta_hat, c(x0,x0,x0), confounders)
      p100=omega(GT$theta_hat, c(x1,x0,x0), confounders)
      p110=omega(GT$theta_hat, c(x1,x1,x0), confounders)
      p111=omega(GT$theta_hat, c(x1,x1,x1), confounders)
      RD1=p100[1]-p000[1]
      RD2=p110[1]-p100[1]
      RD3=p111[1]-p110[1]
      RR1=p100[1]/p000[1]
      RR2=p110[1]/p100[1]
      RR3=p111[1]/p110[1]
      OR1=(p100[1]/(1-p100[1]))/(p000[1]/(1-p000[1]))
      OR2=(p110[1]/(1-p110[1]))/(p100[1]/(1-p100[1]))
      OR3=(p111[1]/(1-p111[1]))/(p110[1]/(1-p110[1]))
      var.boot[boot.count,]=c(RD1,RD2,RD3,RR1,RR2,RR3,OR1,OR2,OR3,NA,NA,NA)
    }else if(intv==4){
      p0000=omega(GT$theta_hat, c(x0,x0,x0,x0), confounders)
      p1000=omega(GT$theta_hat, c(x1,x0,x0,x0), confounders) #first part of difference
      p1100=omega(GT$theta_hat, c(x1,x1,x0,x0), confounders) #first part of difference
      p1110=omega(GT$theta_hat, c(x1,x1,x1,x0), confounders) #first part of difference
      p1111=omega(GT$theta_hat, c(x1,x1,x1,x1), confounders) #first part of difference
      RD1=p1000[1]-p0000[1]
      RD2=p1100[1]-p1000[1]
      RD3=p1110[1]-p1100[1]
      RD4=p1111[1]-p1110[1]
      RR1=p1000[1]/p0000[1]
      RR2=p1100[1]/p1000[1]
      RR3=p1110[1]/p1100[1]
      RR4=p1111[1]/p1110[1]
      OR1=(p1000[1]/(1-p1000[1]))/(p0000[1]/(1-p0000[1]))
      OR2=(p1100[1]/(1-p1000[1]))/(p1000[1]/(1-p1000[1]))
      OR3=(p1110[1]/(1-p1110[1]))/(p1100[1]/(1-p1100[1]))
      OR4=(p1111[1]/(1-p1111[1]))/(p1110[1]/(1-p1110[1]))
      var.boot[boot.count,]=c(RD1,RD2,RD3,RR1,RR2,RR3,OR1,OR2,OR3,RD4,RR4,OR4)
    }
  }
  return(var.boot)
}
