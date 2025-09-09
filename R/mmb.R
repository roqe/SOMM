#' @export
mmb=function(dt,exposure,outcome,mediators,cnfd=NULL,nb=500,mc=5,seed=217){
  covariates=names(dt)[!names(dt)%in%c(exposure,outcome,mediators)]
  pathR=SIM(dt,exposure,outcome,mediators,covariates,cnfd,z0=0,z1=1)

  if(nb>0){
    cl = snow::makeCluster(mc)
    doSNOW::registerDoSNOW(cl)
    pb = txtProgressBar(max = nb, style = 3)
    progress = function(n) setTxtProgressBar(pb, n)
    opts = list(progress = progress)

    var.boot=foreach::foreach(boot.count=1:nb, .options.snow = opts, .combine=rbind,
                              .packages = c("data.table"),
                              .export = c("SIM","BE")) %dopar% {
                                seed=seed*boot.count
                                dt.boot=as.data.frame(dt[sample(1:nrow(dt), replace=TRUE),])
                                pathR.boot=SIM(dt.boot,exposure,outcome,mediators,covariates,cnfd,z0=0,z1=1)
                                return(pathR.boot$PSE)
                              }
    snow::stopCluster(cl)
    cat("\n")
    pathR$PSE=pvbd(stat="mean",path=names(pathR$PSE),effect=pathR$PSE,var.boot=var.boot,nv=0)
  }

  return(pathR)
}

SIM=function(dt,exposure,outcome,mediators,covariates,cnfd=NULL,z0=0,z1=1){
  Z=dt[[exposure]]
  Y=dt[[outcome]]
  if(length(mediators)>1){
    M=dt[mediators]
  }else{
    M=data.frame(dt[[mediators]])
    colnames(M)=mediators
  }
  X=dt[covariates]

  print(paste0("Note: Backward Estimating PSE on ",exposure," > ",paste0(mediators,collapse = " > ")," > ",outcome))
  Xmean=X
  Xmean[names(cnfd)]=cnfd
  RR=unlist(BE(Z,Y,M,X,Xmean,z0,z1))
  names(RR)=paste0("y",gsub(".z","",substr(names(RR),2,nchar(names(RR)[1])-2)))

  nnn=c(unlist(sapply(1:length(mediators),function(i){
    nn=paste0(exposure,">",paste(mediators[i],collapse = ">"),
              if(i<length(mediators)){ paste0(">",paste(sapply((i+1):length(mediators),function(j){
                return(paste0("(",mediators[j],")")) }),collapse = ">")) },
              ">",outcome)
    return(nn)
  })),paste0(exposure,">",outcome),"total")

  ind=2^(0:(length(mediators)+1))
  MD=sapply(1:length(nnn),function(i){
    if(i==length(nnn)){
      (RR[ind[1]]-RR[[ind[i]]])
    }else{
      return(RR[ind[i]]-RR[[ind[i+1]]])
    }
  })
  names(MD)=nnn
  return(list(DELTA=RR[ind],PSE=MD))
}

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
    return(list(z1=MM_z1[1],z0=MM_z0[1]))
  }else if(ncol(M)>1){
    return(list(z1=BE(Z,MM_z1,M[1:(ncol(M)-1)],X,Xmean,z0,z1),z0=BE(Z,MM_z0,M[1:(ncol(M)-1)],X,Xmean,z0,z1)))
  }else{
    return(list(z1=BE(Z,MM_z1,NULL,X,Xmean,z0,z1),z0=BE(Z,MM_z0,NULL,X,Xmean,z0,z1)))
  }
}


