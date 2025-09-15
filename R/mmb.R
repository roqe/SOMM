#' @export
mmb=function(dt,exposure,outcome,mediators,cnfd=NULL,nb=500,mc=5,seed=217,intv="weak",MS=NULL){
  covariates=names(dt)[!names(dt)%in%c(exposure,outcome,mediators)]
  pathR=SIM(dt,exposure,outcome,mediators,covariates,cnfd,z0=0,z1=1,intv,MS)

  if(nb>0){
    cl = snow::makeCluster(mc)
    doSNOW::registerDoSNOW(cl)
    # pb = txtProgressBar(max = nb, style = 3)
    # progress = function(n) setTxtProgressBar(pb, n)
    # opts = list(progress = progress)
    opts = list()

    var.boot=foreach::foreach(boot.count=1:nb, .options.snow = opts, .combine=rbind,
                              .packages = c("data.table"), .export = c("SIM","BE","BEs","PSEnames","bn")) %dopar% {
      seed=seed*boot.count
      ind=sample(1:nrow(dt), replace=TRUE)
      dt.boot=as.data.frame(dt[ind,])
      if(is.null(MS)){ MS.boot=MS }else{ MS.boot=data.frame(MS[ind,]); names(MS.boot)=names(MS) }
      pathR.boot=SIM(dt.boot,exposure,outcome,mediators,covariates,cnfd,z0=0,z1=1,intv,MS.boot)
      return(pathR.boot$PSE)
    }
    snow::stopCluster(cl)
    #cat("\n")
    pathR$PSE=pvbd(stat="mean",path=names(pathR$PSE),effect=pathR$PSE,var.boot=var.boot,nv=0)
  }

  return(pathR)
}

SIM=function(dt,exposure,outcome,mediators,covariates,cnfd,z0,z1,intv,MS){
  Z=dt[[exposure]]
  Y=dt[[outcome]]
  if(length(mediators)>1){
    M=dt[mediators]
  }else{
    M=data.frame(dt[[mediators]])
    colnames(M)=mediators
  }
  X=dt[covariates]

  #print(paste0("Note: Backward Estimating PSE on ",exposure," > ",paste0(mediators,collapse = " > ")," > ",outcome))
  Xmean=X
  Xmean[names(cnfd)]=cnfd
  if(intv=="weak"){
    RR=unlist(BE(Z,Y,M,X,Xmean,z0,z1))
    ind=2^(0:(length(mediators)+1))
  }else{
    RR=unlist(BEs(Z,Y,M,1:ncol(M),X,Xmean,z0,z1,MS))
    ind=2^(0:(2^length(mediators)))
  }
  names(RR)=paste0("y",gsub(".y","",substr(names(RR),2,nchar(names(RR)[1])-2)))
  nnn=PSEnames(exposure,mediators,outcome,intv,re=T)

  MD=sapply(1:length(nnn),function(i){
    if(i==length(nnn)){
      return(RR[ind[1]]-RR[[ind[i]]])
    }else{
      return(RR[ind[i]]-RR[[ind[i+1]]])
    }
  })
  names(MD)=nnn
  nnn=PSEnames(exposure,mediators,outcome,intv,re=F)

  return(list(CO=RR[ind],PSE=MD[nnn]))
}

PSEnames=function(exposure,mediators,outcome,intv,re=F){
  if(intv=="weak"){
    nnn=c(unlist(sapply(1:length(mediators),function(i){
      nn=paste0(exposure,">",paste(mediators[i],collapse = ">"),
                if(i<length(mediators)){ paste0(">",paste(sapply((i+1):length(mediators),function(j){
                  return(paste0("(",mediators[j],")")) }),collapse = ">")) },
                ">",outcome)
      return(nn)
    })),paste0(exposure,">",outcome),"total")
  }else{
    if(re){
      nnn=c(unlist(sapply((2^length(mediators)-1):1,function(i){
        nn=paste0(exposure,">",paste(mediators[rev(as.logical(bn(length(mediators),i)))],collapse = ">"),">",outcome)
        return(nn)
      })),paste0(exposure,">",outcome),"total")
    }else{
      nnn=c(unlist(sapply((2^length(mediators)-1):1,function(i){
        nn=paste0(exposure,">",paste(mediators[(as.logical(bn(length(mediators),i)))],collapse = ">"),">",outcome)
        return(nn)
      })),paste0(exposure,">",outcome),"total")
    }
  }
  return(nnn)
}

### problematic for three and above mediators (finding sequential estimation order)
BE=function(Z,Y,M,X,Xmean,z0,z1){
  if(is.null(M)){
    mY=lm(Y~., data=cbind(Z,X))
    MM_z1=predict(mY,cbind(Z=rep(z1,length(Y)),Xmean)) #Y(za,M1,M2)
    MM_z0=predict(mY,cbind(Z=rep(z0,length(Y)),Xmean))
  }else{
    mY=lm(Y~., data=cbind(Z,M,X))
    MM_z1=predict(mY,cbind(Z=rep(z1,length(Y)),M,Xmean)) #Y(za,M1,M2)
    MM_z0=predict(mY,cbind(Z=rep(z0,length(Y)),M,Xmean))
  }
  if(is.null(M)){
    return(list(y1=MM_z1[1],y0=MM_z0[1]))
  }else if(ncol(M)>1){
    return(list(y1=BE(Z,MM_z1,M[1:(ncol(M)-1)],X,Xmean,z0,z1),y0=BE(Z,MM_z0,M[1:(ncol(M)-1)],X,Xmean,z0,z1)))
  }else{
    return(list(y1=BE(Z,MM_z1,NULL,X,Xmean,z0,z1),y0=BE(Z,MM_z0,NULL,X,Xmean,z0,z1)))
  }
}

BEs=function(Z,Y,M,mind,X,Xmean,z0,z1,MS){
  if(is.null(mind)){
    mY=lm(Y~., data=cbind(Z,X))
    MM_z1=predict(mY,cbind(Z=rep(z1,length(Y)),Xmean)) #Y(za,M1,M2)
    MM_z0=predict(mY,cbind(Z=rep(z0,length(Y)),Xmean))
  }else{
    mY=lm(Y~., data=cbind(Z,M[mind],X))
    MM_z1=predict(mY,cbind(Z=rep(z1,length(Y)),M,Xmean)) #Y(za,M1,M2)
    MM_z0=predict(mY,cbind(Z=rep(z0,length(Y)),M,Xmean))
  }
  if(is.null(mind)){
    return(list(y1=MM_z1[1],y0=MM_z0[1]))
  }else if(length(mind)>1){
    return(list(y1=BEs(Z,MM_z1,M,mind[-1],X,Xmean,z0,z1,MS),y0=BEs(Z,MM_z0,M,mind[-1],X,Xmean,z0,z1,MS)))
  }else if(length(mind)==1&mind>1){
    if(!is.null(MS)){
      return(list(y1=BEs(Z,MM_z1,MS,1:(mind-1),X,Xmean,z0,z1,MS),y0=BEs(Z,MM_z0,MS,1:(mind-1),X,Xmean,z0,z1,MS)))
    }else{
      return(list(y1=BEs(Z,MM_z1,M,1:(mind-1),X,Xmean,z0,z1,MS),y0=BEs(Z,MM_z0,M,1:(mind-1),X,Xmean,z0,z1,MS)))
    }
  }else{
    return(list(y1=BEs(Z,MM_z1,M,NULL,X,Xmean,z0,z1,MS),y0=BEs(Z,MM_z0,M,NULL,X,Xmean,z0,z1,MS)))
  }
}

stTB=function(BB,TVf,TVb){
  BWE=lapply(BB,function(bb){
    return(as.data.table(as.list(setNames(bb$BwE$effect,bb$BwE$path))))
  }) %>% rbindlist()
  FWE=lapply(BB,function(bb){
    return(as.data.table(as.list(setNames(bb$FwE$effect,bb$FwE$path))))
  }) %>% rbindlist()
  BWC=lapply(BB,function(bb){
    return(as.data.table(as.list(setNames(bb$BwE[,4]<TVb&TVb<bb$BwE[,5],bb$BwE$path))))
  }) %>% rbindlist()
  FWC=lapply(BB,function(bb){
    return(as.data.table(as.list(setNames(bb$FwE[,4]<TVf&TVf<bb$FwE[,5],bb$FwE$path))))
  }) %>% rbindlist()
  statFE=lapply(1:ncol(FWE),ss,EE=FWE,CI=FWC,TV=TVf) %>% rbindlist()
  statBE=lapply(1:ncol(BWE),ss,EE=BWE,CI=BWC,TV=TVb) %>% rbindlist()
  apply(statFE,1,function(l){
    cat(paste0(paste0(l,collapse = " & ")," \\\\"))
    cat("\n")
  })
  apply(statBE,1,function(l){
    cat(paste0(paste0(l,collapse = " & ")," \\\\"))
    cat("\n")
  })
  return(list(statBE=statBE,statFE=statFE))
}

ss=function(c,EE,CI,TV){
  options(scipen = -2)
  PSEhat=EE[[c]]
  return(list(
    effect=names(EE)[c],
    true=TV[c],
    RB=signif(mean(PSEhat-TV[c])/TV[c],3),
    SD=signif(sd(PSEhat),3),
    MSE=signif(sum((PSEhat-TV[c])^2)/nrow(EE),3),
    CP=mean(CI[[c]])
  ))
}
