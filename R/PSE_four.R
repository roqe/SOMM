#' @import data.table
PSE_four=function(GT,x0,x1,confounders,V.matrix){
  BF=(summary(GT$reg$y.reg)$family$family=="binomial")
  p0000=omega(GT$theta_hat, c(x0,x0,x0,x0), confounders, BF)
  p1000=omega(GT$theta_hat, c(x1,x0,x0,x0), confounders, BF) #first part of difference
  p1100=omega(GT$theta_hat, c(x1,x1,x0,x0), confounders, BF) #first part of difference
  p1110=omega(GT$theta_hat, c(x1,x1,x1,x0), confounders, BF) #first part of difference
  p1111=omega(GT$theta_hat, c(x1,x1,x1,x1), confounders, BF) #first part of difference
  RD1=rd(p1000,p0000,V.matrix)
  RD2=rd(p1100,p1000,V.matrix)
  RD3=rd(p1110,p1100,V.matrix)
  RD4=rd(p1111,p1110,V.matrix)
  RDT=rd(p1111,p0000,V.matrix)
  RR1=rr(p1000,p0000,V.matrix)
  RR2=rr(p1100,p1000,V.matrix)
  RR3=rr(p1110,p1100,V.matrix)
  RR4=rr(p1111,p1110,V.matrix)
  RRT=rr(p1111,p0000,V.matrix)
  OR1=or(p1000,p0000,V.matrix)
  OR2=or(p1100,p1000,V.matrix)
  OR3=or(p1110,p1100,V.matrix)
  OR4=or(p1111,p1110,V.matrix)
  ORT=or(p1111,p0000,V.matrix)
  omega_values=data.table::data.table(p0000=p0000[1],p1000=p1000[1],p1100=p1100[1],p1110=p1110[1],p1111=p1111[1])
  return(list(omega=omega_values,RD=data.table::data.table(rbind(RD1,RD2,RD3,RD4,RDT)),
              RR=data.table::data.table(rbind(RR1,RR2,RR3,RR4,RRT)),OR=data.table::data.table(rbind(OR1,OR2,OR3,OR4,ORT))))
}
