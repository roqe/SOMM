create_vmatrix=function(GT){
  y.reg=GT$reg$y.reg
  s.reg=GT$reg$s.reg
  q.reg=GT$reg$q.reg
  beta=GT$theta_hat$bc
  alpha=GT$theta_hat$ac
  delta=GT$theta_hat$dc
  lb=length(beta); la=length(alpha); ld=length(delta)
  dimensions=sum(lb+la+ld)
  V.matrix=matrix(0, nr=dimensions, nc=dimensions)
  V.matrix[c(1:lb),c(1:lb)]=check_vmatrix(y.reg,beta)
  V.matrix[c((1+lb):(lb+la)),c((1+lb):(lb+la))]=check_vmatrix(s.reg,alpha)
  V.matrix[c((1+lb+la):(lb+la+ld)),c((1+lb+la):(lb+la+ld))]=check_vmatrix(q.reg,delta)
  return(V.matrix)
}

check_vmatrix=function(reg,theta){
  ll=length(theta)
  names(theta)[1]="(Intercept)"
  tt=matrix(0, nr=ll, nc=ll, dimnames=list(names(theta),names(theta)))
  if(sum(theta)!=0){
    vv=as.matrix(vcov(reg))
    on=dimnames(vv)[[1]]
    on=on[on%in%names(theta)]
    tt[on,on]=vv[on,on]
  }
  return(tt)
}
