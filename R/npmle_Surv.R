predit.haz=function(x,R,Y){
  Ys=sort(Y)
  ph=apply(matrix(x),1,function(t){
    cid=sum(t>=Ys)
    if(cid<1) ans=0
    if(1<=cid & cid<length(Ys)) ans=R[cid]
    if(length(Ys)<=cid) ans=R[length(Ys)]
    return(ans)
  })
  return(ph)
}

tran.npmle_est=function(Y,D,X,p.time){
  p=length(X)/length(Y)
  n=length(Y)

  B=survival::coxph(Surv(Y,D)~X,method="breslow")$coef
  names(B)=colnames(X)

  Y=sort(Y,index.return=TRUE)
  D=D[Y$ix]
  if(p==1){X=matrix(X[Y$ix],n,p)}
  if(p>1){X=matrix(X[Y$ix,],n,p)}
  Y=Y$x

  Yp=unique(Y)
  tiefun_1=function(cid){ sum(Y[D==1]==cid) }
  tiefun_2=function(cid){ sum(Y==cid) }

  nYjp=apply(matrix(Yp),1,tiefun_1)
  Y.times=apply(matrix(Yp),1,tiefun_2)
  lc.jump=which(diff(c(0,Y))!=0)

  np=length(Yp)
  R=((1:np)/(np))*max(Y)

  ###########################  transformation function G ##################################
  #########################################################################################
  G=function(t){ return(-log(plnorm(t,meanlog = 0,sdlog = 1,lower.tail=FALSE))) }
  g=function(t){ return(dlnorm(t,sdlog = 1)/(plnorm(t,meanlog = 0,sdlog = 1,lower.tail=FALSE))) }

  dg=function(t){
    temp=t
    sig2=1^2
    sig = sqrt(sig2)
    Phi = pnorm(log(temp),0,sig)
    temp2 = dnorm(log(temp),0,sig)
    phi = temp2/temp
    dphi = -temp2/(temp^2)-((log(temp))/(temp*sig2))*temp2/temp
    d2phi = 2*(temp^-3)*temp2+2*(temp^-2)*temp2*((log(temp))/(temp*sig2))+(temp^-1)*temp2*((log(temp))/(temp*sig2))^2-(temp^-1)*temp2*(1/(sig2*(temp^2)))+(temp^-1)*temp2*(log(temp))/((temp^2)*sig2)
    d2G = dphi/(1-Phi)+(phi/(1-Phi))^2
    return(d2G)
  }

  kpa=function(t){ return(dg(t)/g(t)) }
  #############################################################################
  # eat R
  #############################################################################

  estR=function(B,R){
    dRk=diff(c(0,R))
    if(p==1){BX=c(X)*c(B)}                             # \beta^{'} X_i, i=1,2,3,...,n
    if(p>1){BX=apply(matrix(t(t(X)*B),n,p),1,sum)}      # \beta^{'} X_i, i=1,2,3,...,n
    eBX=exp(BX)
    BXM=matrix(BX,n,np)
    eBXM=exp(BXM)
    up=nYjp
    adjg=g(R*eBX)
    adjK_0=kpa(R*eBX)
    adjK_1=ifelse(adjK_0=="NaN",0,adjK_0)
    adjK_2=ifelse(adjK_1>adjg,0,adjK_1)
    down=cumsum((eBX*(adjg-adjK_2*D))[n:1])[n:1]
    equ=cumsum(up/down)
    return(equ)
  }

  estB=function(B,R){
    if(p==1){BX=X*c(B)}
    if(p>1){BX=apply(matrix(t(t(X)*B),n,p),1,sum)}
    eBX=exp(BX)
    BXM=matrix(BX,n,n)
    eBXM=exp(BXM)
    RM=matrix(R,n,n,byrow=TRUE)
    adjk_b0=kpa(R*eBX)
    adjk_b1=ifelse(adjk_b0=="NaN",0,adjk_b0)
    adjk_b1=ifelse(adjk_b1<0,adjk_b1,0)
    part2=-g(R*eBX)*R*eBX
    part1=(adjk_b1*R*eBX+1)*D
    partall=part1+part2
    equ=apply(X,2,function(xx){ return(sum(partall*xx)) })
    return(equ)
  }

  newton=function(B,R){
    eva=1/n
    goal=estB(B,R)
    JM=matrix(NA,p,p)
    for(j in 1:p){
      Be=B
      Be[j]=Be[j]+eva
      JM[j,]=(estB(Be,R)-goal)/eva
    }
    rse<-list(JM=JM,goal=goal)
    return(rse)
  }
  ###############################################################################
  R0=0
  RRR=0

  while(abs(sum(R0)-sum(R))>10^-3&RRR<50){
    R0=R
    R=estR(B,R0)
    RRR=RRR+1
  }
  goal=estB(B,R)

  if(sum(goal^2)^0.5>10^(-3)){
    BBnew=B
    BBold=rep(Inf,p)
  }
  if(sum(goal^2)^0.5<10^(-3)){
    BBold=BBnew=B
  }

  ###############################################################################
  kk=0
  while((max(abs(BBold-BBnew))>10^-3|sum(goal^2)^0.5>10^(-6))&kk<10){
    BBold=BBnew
    rrr=0
    while(sum(goal^2)^0.5>10^(-6)&rrr<50){
      resJM=newton(BBnew,R)
      if(p==1){BBnew=BBnew-c(1/(resJM$JM)*resJM$goal)}
      if(p>1){BBnew=BBnew-c(solve(resJM$JM)%*%resJM$goal)/10}
      BBnew=BBnew-c(solve(resJM$JM)%*%resJM$goal)
      goal=resJM$goal
      rrr=rrr+1
    }

    RRR=0
    while(abs(sum(R0)-sum(R))>10^-3&RRR<10){
      R0=R
      R=estR(BBnew,R0)
      RRR=RRR+1
    }
    goal=estB(BBnew,R)
    kk=kk+1
  }

  ptime.Lm=predit.haz(p.time,R,Y)
  rse<-list(coef=BBnew,bas.haz=R,iter=kk,event.time=Y[D==1],
            ptime.Lm=ptime.Lm,ptime.Ht=log(ptime.Lm))
  return(rse)
}
