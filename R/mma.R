mma=function(dt,exposure,outcome,mediators,cnfd=NULL){
  Z=dt[[exposure]]
  Y=dt[[outcome]]
  M=dt[mediators]
  X=dt[!colnames(dt)%in%c(exposure,outcome,mediators)]
  print(paste0("Note: Estimating PSE on ",exposure," > ",paste0(mediators,collapse = " > ")," > ",outcome))
  total.effect=summary(glm(Y~.,data = cbind(Z,X)))$coef[,1]

  regR=lapply(0:ncol(M),function(i){
      if(i==ncol(M)){
        return(summary(lm(Y~., data=cbind(Z,M,X)))$coef[,1])
      }else if(i==0){
        return(summary(lm(M[[i+1]]~., data=cbind(Z,X)))$coef[,1])
      }else{
        return(summary(lm(M[[i+1]]~., data=cbind(Z,M[1:i],X)))$coef[,1])
      }
    })

  invR=lapply(0:(2^ncol(M)),function(s){
    zz=c(rep(1,s),rep(0,2^(ncol(M))-s))
    mu=muu(regR,zz,cnfd)
    names(mu)=paste0("z",paste0(zz,collapse = ""))
    return(mu)
  })

}

muu=function(RR,zz,cnfd,w=0){
  inK=length(RR)
  theta=RR[[inK]]
  muj=theta["(Intercept)"]+theta[["Z"]]*zz[w+1]+sum(theta[names(X)]*cnfd)
  if(inK>1){
    cmu=sapply(1:(inK-1),function(k){
      return(muu(RR[1:k],zz,cnfd,w=w+2^(k-1)))
    })
    return(muj+sum(theta[mediators][1:length(cmu)]*cmu))
  }else{
    return(muj)
  }
}
