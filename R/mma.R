#' @export
mma=function(dt,exposure,outcome,mediators,cnfd=NULL,nb=500,mc=5,seed=217,cid=NULL,ntp=100,oti=c(.25,.5,.75)){
  covariates=names(dt)[!names(dt)%in%c(exposure,outcome,mediators,cid)]
  regR=thetaHAT(dt,exposure,outcome,mediators,covariates,cid,ntp)
  BF=ifelse(all(dt$Y%in%c(0,1)),T,F)
  pathR=PSE(regR,cnfd,BF,exposure,outcome,mediators,covariates,cid)

  if(nb>0){
    cl = snow::makeCluster(mc)
    doSNOW::registerDoSNOW(cl)
    pb = txtProgressBar(max = nb, style = 3)
    progress = function(n) setTxtProgressBar(pb, n)
    opts = list(progress = progress)
    var.boot=foreach::foreach(boot.count=1:nb, .options.snow = opts, .combine=rbind,
                              .packages = c("survival","data.table"),
                              .export = c("PSE","thetaHAT","muu","vqq","bn","regR","stat",
                                          "tran.npmle_est","predit.haz")) %dopar% {
      seed=seed*boot.count
      dt.boot=as.data.frame(dt[sample(1:nrow(dt), replace=TRUE),])
      regR.boot=thetaHAT(dt.boot,exposure,outcome,mediators,covariates,cid,ntp)
      pathR.boot=PSE(regR.boot,cnfd,BF,exposure,outcome,mediators,covariates,cid)
      return(pathR.boot$PSE)
    }
    cat("\n")
    if(BF){
      ll=ifelse(is.list(pathR$PSE),length(pathR$PSE),1)
      pb = txtProgressBar(max = ll, style = 3)
      progress = function(n) setTxtProgressBar(pb, n)
      pseR=foreach::foreach(t=1:ll, .options.snow = opts, .export = c("pvbd","btbd")) %dopar% {
        if(ll==1){
          pse=pathR$PSE
          VBt=var.boot
        }else{
          pse=pathR$PSE[[t]]
          VBt=do.call(rbind,var.boot[,t])
        }
        path=colnames(pse)
        pseRt=rbind(pvbd(stat="RD",path=path,effect=pse["pse_RD",],var.boot=VBt[rownames(VBt)=="pse_RD",],nv=0),
                    pvbd(stat="RR",path=path,effect=pse["pse_RR",],var.boot=VBt[rownames(VBt)=="pse_RR",],nv=1),
                    pvbd(stat="OR",path=path,effect=pse["pse_OR",],var.boot=VBt[rownames(VBt)=="pse_OR",],nv=1))
        return(pseRt)
      }
    }else{
      pseR=list(pvbd(stat="std_mean",path=names(pathR$PSE),effect=pathR$PSE,var.boot=var.boot,nv=0))
    }
    snow::stopCluster(cl)
    if(!is.null(cid)){
      mm=quantile(1:nrow(pathR$OMEGA),oti)
      ps=pseR[mm]
      names(ps)=paste0("pse_",round(mm))
      return(list(OMEGA=cbind(quantile=mm,pathR$OMEGA[mm,]),PSE=ps))
    }else{
      return(list(OMEGA=pathR$OMEGA,PSE=pseR[[1]]))
    }
  }
  return(pathR)
}

pvbd=function(stat,path,effect,var.boot,nv){
  ulbd=apply(var.boot,2,btbd,nv)
  pseR=data.table(stat=stat,path=path,effect=effect,
                  `lower(b)`=ulbd[1,],`upper(b)`=ulbd[2,],`pv(b)`=ulbd[3,],
                  `lower(n)`=ulbd[4,],`upper(n)`=ulbd[5,],`pv(n)`=ulbd[6,])
  return(pseR)
}

btbd=function(bb,nv){
  vv=mean(bb>nv)
  pb=min(vv,(1-vv))*2
  exact=c(quantile(bb,c(.05,.95)),pb)
  cc=rnorm(1e+6, mean=mean(bb,na.rm = T), sd=sd(bb,na.rm = T))
  vv=mean(cc>nv);
  pa=min(vv,(1-vv))*2
  appxn=c(quantile(cc,c(.05,.95)),pa)
  return(c(exact,appxn))
}

thetaHAT=function(dt,exposure,outcome,mediators,covariates,cid,ntp){
  Z=dt[[exposure]]
  Y=dt[[outcome]]
  if(length(mediators)>1){
    M=dt[mediators]
  }else{
    M=data.frame(dt[[mediators]])
    colnames(M)=mediators
  }
  X=dt[covariates]
  if(!is.null(cid)){ D=dt[[cid]] }
  print(paste0("Note: Forward Estimating PSE on ",exposure," > ",paste0(mediators,collapse = " > ")," > ",outcome))

  regR=lapply(0:ncol(M),function(i){
    if(i==ncol(M)){
      if(!is.null(cid)){
        print("Survival outcome, probit link.")
        return(tran.npmle_est(Y,D,X=as.matrix(cbind(Z,M,X)),p.time = seq(min(Y),max(Y),length.out=ntp)))
      }else if(all(Y%in%c(0,1))){
        print("Binary outcome, probit link.")
        return(glm(Y~.,data = cbind(Z,M,X),family = binomial(link="probit")))
      }else{
        print("Continuous outcome, identity link.")
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

PSE=function(regR,cnfd,BF,exposure,outcome,mediators,covariates,cid){
  if(BF){ varYY=1 }else{ varYY=mean(resid(regR[[length(regR)]])^2) }
  varMY=sapply(1:(length(regR)-1),function(k){
    varM=mean(resid(regR[[k]])^2)*vqq(k,regR,mediators)
  })
  sdY=sqrt(varYY+sum(varMY))
  invZ=sapply(0:(2^length(mediators)),function(s){
    zz=c(rep(1,s),rep(0,2^(length(mediators))-s))
    mu=muu(regR,zz,cnfd,w=0,mediators,covariates)
    if(!is.null(cid)){
      mu=data.table(mu+regR[[length(regR)]]$ptime.Ht)
      names(mu)=paste0("Ht",paste0(zz,collapse = ""))
    }else{
      names(mu)=paste0(ifelse(BF,"p","z"),paste0(zz,collapse = ""))
    }
    return(mu/sdY)
  })
  if(is.list(invZ)){ invZ=rbindlist(list(invZ)) }

  nnn=c(paste0(exposure,">",outcome),unlist(sapply(1:(2^length(mediators)-1),function(i){
    nn=paste0(exposure,">",paste(mediators[rev(as.logical(bn(length(mediators),i)))],collapse = ">"),">",outcome)
    return(nn)
  })),"total")
  if(BF){
    if(!is.data.table(invZ)){ pse=stat(invZ,nnn) }else{ pse=apply(invZ,1,stat,nnn,simplify = F) }
  }else{
    pse=c(diff(invZ),invZ[[length(invZ)]]-invZ[[1]])
    names(pse)=nnn
  }

  return(list(OMEGA=invZ,PSE=pse))
}

stat=function(invZ,nnn){
  invZ=pnorm(invZ)
  p1=invZ[[length(invZ)]]
  pse_RD=c(diff(invZ),p1-invZ[[1]])
  frn=invZ[2:length(invZ)]
  bck=invZ[1:(length(invZ)-1)]
  pse_RR=c(frn/bck,p1/invZ[[1]])
  pse_OR=c((frn/(1-frn))/(bck/(1-bck)),(p1/(1-p1))/(invZ[[1]]/(1-invZ[[1]])))
  pse=rbind(pse_RD,pse_RR,pse_OR)
  colnames(pse)=nnn
  return(pse)
}

muu=function(RR,zz,cnfd,w=0,mediators,covariates){
  inK=length(RR)
  if(length(RR[[inK]])==6){
    theta=RR[[inK]]$coef
    muJ=theta[["Z"]]*zz[w+1]+sum(theta[covariates]*cnfd)
  }else{
    theta=summary(RR[[inK]])$coef[,1]
    muJ=theta["(Intercept)"]+theta[["Z"]]*zz[w+1]+sum(theta[covariates]*cnfd)
  }
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
  CM=regR[[length(regR)]]
  if(length(CM)==6){
    thetaJK=CM$coef[mediators][k]
  }else{
    thetaJK=summary(CM)$coef[,1][mediators][k]
  }
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
