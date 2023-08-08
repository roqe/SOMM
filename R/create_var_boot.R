#' @import foreach
create_var_boot=function(nb, dt, confounders=c(), intv=3, reNAME=NULL, x0, x1, stra=1){
  if(stra>1){
    dt.BB=c()
    for (i in 1:stra){
      dt.Q=dt[dt$Q<=i & dt$Q>(i-1),]
      dt.B=c()
      for (j in 1:stra){
        w_range=which(dt.Q$W<=j & dt.Q$W>(j-1))
        dt.b=as.data.frame(dt.Q[w_range,][sample(1:sum(w_range), replace=T),])
        dt.B=rbind(dt.B, dt.b)
      }
      dt.BB=rbind(dt.BB,dt.B)
    }
    if (!is.null(reNAME)) {
      dt.BB=dt.BB[-is.na(dt.BB$id),]
      a=table(dt.BB$id)
      dt.boot=c()
      for (AA in min(a):max(a)){
        b=dt[which(dt$id %in% names(which(a>=AA))),]
        dt.boot=rbind(dt.boot,b)
      }
    }
    if (is.null(reNAME)) dt.boot=dt.BB[!is.na(dt.BB$Y),]
  }else{
    dt.boot=as.data.frame(dt[sample(1:nrow(dt), replace=TRUE),])
  }

    GT=get_theta(dt.boot,reNAME)
    BF=(summary(GT$reg$y.reg)$family=="binomial")
    if(intv==3){
      p000=omega(GT$theta_hat, c(x0,x0,x0), confounders, BF)
      p100=omega(GT$theta_hat, c(x1,x0,x0), confounders, BF)
      p110=omega(GT$theta_hat, c(x1,x1,x0), confounders, BF)
      p111=omega(GT$theta_hat, c(x1,x1,x1), confounders, BF)
      RD1=p100[1]-p000[1]
      RD2=p110[1]-p100[1]
      RD3=p111[1]-p110[1]
      RR1=p100[1]/p000[1]
      RR2=p110[1]/p100[1]
      RR3=p111[1]/p110[1]
      OR1=(p100[1]/(1-p100[1]))/(p000[1]/(1-p000[1]))
      OR2=(p110[1]/(1-p110[1]))/(p100[1]/(1-p100[1]))
      OR3=(p111[1]/(1-p111[1]))/(p110[1]/(1-p110[1]))
      RDT=p111[1]-p000[1]
      RRT=p111[1]/p000[1]
      ORT=(p111[1]/(1-p111[1]))/(p000[1]/(1-p000[1]))
      return(list(RD1,RD2,RD3,RR1,RR2,RR3,OR1,OR2,OR3,NA,NA,NA,RDT,RRT,ORT))
    }else if(intv==4){
      p0000=omega(GT$theta_hat, c(x0,x0,x0,x0), confounders, BF)
      p1000=omega(GT$theta_hat, c(x1,x0,x0,x0), confounders, BF) #first part of difference
      p1100=omega(GT$theta_hat, c(x1,x1,x0,x0), confounders, BF) #first part of difference
      p1110=omega(GT$theta_hat, c(x1,x1,x1,x0), confounders, BF) #first part of difference
      p1111=omega(GT$theta_hat, c(x1,x1,x1,x1), confounders, BF) #first part of difference
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
      RDT=p1111[1]-p0000[1]
      RRT=p1111[1]/p0000[1]
      ORT=(p1111[1]/(1-p1111[1]))/(p0000[1]/(1-p0000[1]))
      return(list(RD1,RD2,RD3,RR1,RR2,RR3,OR1,OR2,OR3,RD4,RR4,OR4,RDT,RRT,ORT))
    }
 }
