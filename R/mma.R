#' @export
mma=function(dt,exposure,outcome,mediators,cnfd=NULL,nb=500,mc=5,seed=217){
  covariates=names(dt)[!names(dt)%in%c(exposure,outcome,mediators)]
  regR=thetaHAT(dt,exposure,outcome,mediators,covariates)
  BF=ifelse(all(dt$Y%in%c(0,1)),T,F)
  pathR=PSE(regR,cnfd,BF,exposure,outcome,mediators,covariates)

  if(nb>0){
    cl = snow::makeCluster(mc)
    doSNOW::registerDoSNOW(cl)
    pb = txtProgressBar(max = nb, style = 3)
    progress = function(n) setTxtProgressBar(pb, n)
    opts = list(progress = progress)
    var.boot=foreach::foreach(boot.count=1:nb, .options.snow = opts, .combine=rbind,
                              .export = c("PSE","thetaHAT","muu","vqq","bn","regR")) %dopar% {
      seed=seed*boot.count
      dt.boot=as.data.frame(dt[sample(1:nrow(dt), replace=TRUE),])
      regR.boot=thetaHAT(dt.boot,exposure,outcome,mediators,covariates)
      pathR.boot=PSE(regR.boot,cnfd,BF,exposure,outcome,mediators,covariates)
      return(pathR.boot$pse)
    }
    snow::stopCluster(cl)
  }
  if(BF){
    pseR=rbind(pvbd(stat="RD",path=colnames(pathR$pse),effect=pathR$pse["pse_RD",],var.boot=var.boot[rownames(var.boot)=="pse_RD",]),
               pvbd(stat="RR",path=colnames(pathR$pse),effect=pathR$pse["pse_RD",],var.boot=var.boot[rownames(var.boot)=="pse_RD",]),
               pvbd(stat="OR",path=colnames(pathR$pse),effect=pathR$pse["pse_RD",],var.boot=var.boot[rownames(var.boot)=="pse_RD",]))
  }else{
    pseR=pvbd(stat="mean",path=names(pathR$pse),effect=pathR$pse,var.boot=var.boot)
  }
  return(list(OMEGA=pathR$invZ,PSE=pseR))
}

pvbd=function(stat,path,effect,var.boot){
  ulbd=apply(var.boot,2,btbd)
  pseR=data.table(stat=stat,path=path,effect=effect,
                  `lower(b)`=ulbd[1,],`upper(b)`=ulbd[2,],`pv(b)`=ulbd[3,],
                  `lower(n)`=ulbd[4,],`upper(n)`=ulbd[5,],`pv(n)`=ulbd[6,])
  return(pseR)
}

btbd=function(bb){
  vv=mean(bb>0)
  pb=min(vv,(1-vv))*2
  exact=c(quantile(bb,c(.05,.95)),pb)
  cc=rnorm(1e+6, mean=mean(bb,na.rm = T), sd=sd(bb,na.rm = T))
  vv=mean(cc>0);
  pa=min(vv,(1-vv))*2
  appxn=c(quantile(cc,c(.05,.95)),pa)
  return(c(exact,appxn))
}

thetaHAT=function(dt,exposure,outcome,mediators,covariates){
  Z=dt[[exposure]]
  Y=dt[[outcome]]
  M=dt[mediators]
  X=dt[covariates]
  print(paste0("Note: Estimating PSE on ",exposure," > ",paste0(mediators,collapse = " > ")," > ",outcome))

  BF=ifelse(all(Y%in%c(0,1)),T,F)
  if(BF){
    crude=glm(Y~.,data = cbind(Z,X),family = binomial(link="probit"))
  }else{
    crude=lm(Y~.,data = cbind(Z,X))
  }

  regR=lapply(0:ncol(M),function(i){
    if(i==ncol(M)){
      if(BF){
        return(glm(Y~.,data = cbind(Z,M,X),family = binomial(link="probit")))
      }else{
        return(lm(Y~., data=cbind(Z,M,X)))
      }
    }else if(i==0){
      return(lm(M[[i+1]]~., data=cbind(Z,X)))
    }else{
      return(lm(M[[i+1]]~., data=cbind(Z,M[1:i],X)))
    }
  })

  return(regR)
}

bn=function(j,i){
  bnv=rev(as.numeric(intToBits(i)))
  return(bnv[(length(bnv)-j+1):length(bnv)])
}

PSE=function(regR,cnfd,BF,exposure,outcome,mediators,covariates){
  if(BF){ varYY=1 }else{ varYY=mean(resid(regR[[length(regR)]])^2) }
  varMY=sapply(1:(length(regR)-1),function(k){
    varM=mean(resid(regR[[k]])^2)*vqq(k,regR,mediators)
  })
  sdY=sqrt(varYY+sum(varMY))
  invZ=sapply(0:(2^length(mediators)),function(s){
    zz=c(rep(1,s),rep(0,2^(length(mediators))-s))
    mu=muu(regR,zz,cnfd,w=0,mediators,covariates)
    names(mu)=paste0(ifelse(BF,"p","z"),paste0(zz,collapse = ""))
    return(mu/sdY)
  })

  nnn=c(paste0(exposure,">",outcome),unlist(sapply(1:(2^length(mediators)-1),function(i){
    nn=paste0(exposure,">",paste(mediators[rev(as.logical(bn(length(mediators),i)))],collapse = ">"),">",outcome)
    return(nn)
  })),"total")
  if(BF){
    invZ=pnorm(invZ)
    p1=invZ[[length(invZ)]]
    pse_RD=c(diff(invZ),p1-invZ[[1]])
    frn=invZ[2:length(invZ)]
    bck=invZ[1:(length(invZ)-1)]
    pse_RR=c(frn/bck,p1/invZ[[1]])
    pse_OR=c((frn/(1-frn))/(bck/(1-bck)),(p1/(1-p1))/(invZ[[1]]/(1-invZ[[1]])))
    pse=rbind(pse_RD,pse_RR,pse_OR)
    colnames(pse)=nnn
  }else{
    pse=c(diff(invZ),invZ[[length(invZ)]]-invZ[[1]])
    names(pse)=nnn
  }

  return(list(invZ=invZ,pse=pse))
}

muu=function(RR,zz,cnfd,w=0,mediators,covariates){
  inK=length(RR)
  theta=summary(RR[[inK]])$coef[,1]
  muJ=theta["(Intercept)"]+theta[["Z"]]*zz[w+1]+sum(theta[covariates]*cnfd)
  if(inK>1){
    cmu=sapply(1:(inK-1),function(k){
      return(muu(RR[1:k],zz,cnfd,w=w+2^(k-1),mediators,covariates))
    })
    sumM=sum(theta[mediators][1:length(cmu)]*cmu)
  }else{
    sumM=0
  }
  return(muJ+sumM)
}

vqq=function(k,regR,mediators){
  thetaJK=summary(regR[[length(regR)]])$coef[,1][mediators][k]
  if((k+1)>(length(regR)-1)){
    sumKI=0
  }else{
    sumKI=sapply((k+1):(length(regR)-1),function(i){
      thetaIK=summary(regR[[i]])$coef[,1][names(thetaJK)]
      return(thetaIK^2*vqq(i,regR,mediators))
    })
  }
  return(thetaJK^2+sum(sumKI))
}
