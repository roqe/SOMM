#' @import data.table
PSE_four=function(GT,x0,x1,cnfd,cnfd2,cnfd3,V.matrix,BF){
  p0000=omega(GT$theta_hat, c(x0,x0,x0,x0), cnfd,cnfd2,cnfd3, BF)
  p1000=omega(GT$theta_hat, c(x1,x0,x0,x0), cnfd,cnfd2,cnfd3, BF) #first part of difference
  p1100=omega(GT$theta_hat, c(x1,x1,x0,x0), cnfd,cnfd2,cnfd3, BF) #first part of difference
  p1110=omega(GT$theta_hat, c(x1,x1,x1,x0), cnfd,cnfd2,cnfd3, BF) #first part of difference
  p1111=omega(GT$theta_hat, c(x1,x1,x1,x1), cnfd,cnfd2,cnfd3, BF) #first part of difference
  omega_values=data.table::data.table(p0000=p0000[1],p1000=p1000[1],p1100=p1100[1],p1110=p1110[1],p1111=p1111[1])

  RD1=rd(p1000,p0000,V.matrix); RD2=rd(p1100,p1000,V.matrix); RD3=rd(p1110,p1100,V.matrix); RD4=rd(p1111,p1110,V.matrix); RDT=rd(p1111,p0000,V.matrix)
  RDs=ifelse(abs(sum(sign(c(RD1[1],RD2[1],RD3[1],RD4[1]))))==sum(c(RD1[1],RD2[1],RD3[1],RD4[1])!=0),T,F)
  RD=data.table::data.table(rbind(RD1,RD2,RD3,RD4,RDT),
                            PM=c(ifelse(RDs,RD1[1]/RDT[1],NA),ifelse(RDs,RD2[1]/RDT[1],NA),
                                 ifelse(RDs,RD3[1]/RDT[1],NA),ifelse(RDs,RD4[1]/RDT[1],NA),1))
  RR=OR=NULL
  if(BF){
    RR1=rr(p1000,p0000,V.matrix); RR2=rr(p1100,p1000,V.matrix); RR3=rr(p1110,p1100,V.matrix); RR4=rr(p1111,p1110,V.matrix); RRT=rr(p1111,p0000,V.matrix)
    OR1=or(p1000,p0000,V.matrix); OR2=or(p1100,p1000,V.matrix); OR3=or(p1110,p1100,V.matrix); OR4=or(p1111,p1110,V.matrix); ORT=or(p1111,p0000,V.matrix)
    RRs=ifelse(abs(sum(sign(log(c(RR1[1],RR2[1],RR3[1],RR4[1])))))==sum(log(c(RR1[1],RR2[1],RR3[1],RR4[1]))!=0),T,F)
    ORs=ifelse(abs(sum(sign(log(c(OR1[1],OR2[1],OR3[1],OR4[1])))))==sum(log(c(OR1[1],OR2[1],OR3[1],OR4[1]))!=0),T,F)
    RR=data.table::data.table(rbind(RR1,RR2,RR3,RR4,RRT),
                              PM=c(ifelse(RRs,log(RR1[1])/log(RRT[1]),NA),ifelse(RRs,log(RR2[1])/log(RRT[1]),NA),
                                   ifelse(RRs,log(RR3[1])/log(RRT[1]),NA),ifelse(RRs,log(RR4[1])/log(RRT[1]),NA),
                                   ifelse(RRs,log(RRT[1])/log(RRT[1]),NA)))
    OR=data.table::data.table(rbind(OR1,OR2,OR3,OR4,ORT),
                              PM=c(ifelse(ORs,log(OR1[1])/log(ORT[1]),NA),ifelse(ORs,log(OR2[1])/log(ORT[1]),NA),
                                   ifelse(ORs,log(OR3[1])/log(ORT[1]),NA),ifelse(ORs,log(OR4[1])/log(ORT[1]),NA),
                                   ifelse(ORs,log(ORT[1])/log(ORT[1]),NA)))
  }

  return(list(omega=omega_values,RD=RD,RR=RR,OR=OR))
}
