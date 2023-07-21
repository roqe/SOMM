PSE_three=function(GT,x0,x1,confounders,V.matrix){
  p000=omega(GT$theta_hat, c(x0,x0,x0), confounders)
  p100=omega(GT$theta_hat, c(x1,x0,x0), confounders)
  p110=omega(GT$theta_hat, c(x1,x1,x0), confounders)
  p111=omega(GT$theta_hat, c(x1,x1,x1), confounders)
  RD1=rd(p100,p000,V.matrix)
  RD2=rd(p110,p100,V.matrix)
  RD3=rd(p111,p110,V.matrix)
  RR1=rr(p100,p000,V.matrix)
  RR2=rr(p110,p100,V.matrix)
  RR3=rr(p111,p110,V.matrix)
  OR1=rr(p100/(1-p100),p000/(1-p000),V.matrix)
  OR2=rr(p110/(1-p110),p100/(1-p100),V.matrix)
  OR3=rr(p111/(1-p111),p110/(1-p110),V.matrix)
  omega_values=data.table(p000=p000[1],p100=p100[1],p110=p110[1],p111=p111[1])

  return(list(omega=omega_values,RD=data.table(rbind(RD1,RD2,RD3)),
              RR=data.table(rbind(RR1,RR2,RR3)),OR=data.table(rbind(OR1,OR2,OR3))))
}
