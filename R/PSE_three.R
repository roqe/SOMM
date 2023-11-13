#' @import data.table
PSE_three=function(GT,x0,x1,cnfd,cnfd2,cnfd3,V.matrix,BF){
  p000=omega(GT$theta_hat, c(x0,x0,x0), cnfd,cnfd2,cnfd3, BF)
  p100=omega(GT$theta_hat, c(x1,x0,x0), cnfd,cnfd2,cnfd3, BF)
  p110=omega(GT$theta_hat, c(x1,x1,x0), cnfd,cnfd2,cnfd3, BF)
  p111=omega(GT$theta_hat, c(x1,x1,x1), cnfd,cnfd2,cnfd3, BF)
  omega_values=data.table::data.table(p000=p000[1],p100=p100[1],p110=p110[1],p111=p111[1])

  RD1=rd(p100,p000,V.matrix); RD2=rd(p110,p100,V.matrix); RD3=rd(p111,p110,V.matrix); RDT=rd(p111,p000,V.matrix)
  RDs=ifelse(abs(sum(sign(c(RD1[1],RD2[1],RD3[1]))))==sum(c(RD1[1],RD2[1],RD3[1])!=0),T,F)
  RD=data.table::data.table(rbind(RD1,RD2,RD3,RDT),
                            PM=c(ifelse(RDs,RD1[1]/RDT[1],NA),ifelse(RDs,RD2[1]/RDT[1],NA),ifelse(RDs,RD3[1]/RDT[1],NA),1))
  RR=OR=NULL
  if(BF){
    RR1=rr(p100,p000,V.matrix); RR2=rr(p110,p100,V.matrix); RR3=rr(p111,p110,V.matrix); RRT=rr(p111,p000,V.matrix)
    OR1=or(p100,p000,V.matrix); OR2=or(p110,p100,V.matrix); OR3=or(p111,p110,V.matrix); ORT=or(p111,p000,V.matrix)
    RRs=ifelse(abs(sum(sign(log(c(RR1[1],RR2[1],RR3[1])))))==sum(log(c(RR1[1],RR2[1],RR3[1]))!=0),T,F)
    ORs=ifelse(abs(sum(sign(log(c(OR1[1],OR2[1],OR3[1])))))==sum(log(c(OR1[1],OR2[1],OR3[1]))!=0),T,F)
    RR=data.table::data.table(rbind(RR1,RR2,RR3,RRT),
                              PM=c(ifelse(RRs,log(RR1[1])/log(RRT[1]),NA),ifelse(RRs,log(RR2[1])/log(RRT[1]),NA),
                                   ifelse(RRs,log(RR3[1])/log(RRT[1]),NA),ifelse(RRs,log(RRT[1])/log(RRT[1]),NA)))
    OR=data.table::data.table(rbind(OR1,OR2,OR3,ORT),
                              PM=c(ifelse(ORs,log(OR1[1])/log(ORT[1]),NA),ifelse(ORs,log(OR2[1])/log(ORT[1]),NA),
                                   ifelse(ORs,log(OR3[1])/log(ORT[1]),NA),ifelse(ORs,log(ORT[1])/log(ORT[1]),NA)))
  }
  return(list(omega=omega_values,RD=RD,RR=RR,OR=OR))
}
