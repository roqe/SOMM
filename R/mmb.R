#' @export
mmb=function(dt,exposure,outcome,mediators,cnfd=NULL,nb=500,mc=5,seed=217,intv="weak",MS=NULL,model="linear",std_data=T,z0=0,z1=1,bwr=F){
  covariates=names(dt)[!names(dt)%in%c(exposure,outcome,mediators)]
  if(std_data){
    dt_ori=dt
    dt=as.data.frame(scale(dt))
  }
  pathR=npPSE(dt,exposure,outcome,mediators,covariates,cnfd,z0,z1,intv,MS,model,bw=NULL)

  if(nb>0){
    cl = snow::makeCluster(mc)
    doSNOW::registerDoSNOW(cl)
    # pb = txtProgressBar(max = nb, style = 3)
    # progress = function(n) setTxtProgressBar(pb, n)
    # opts = list(progress = progress)
    # opts = list()

    var.boot=foreach::foreach(boot.count=1:nb, # .options.snow = opts,
                              .packages = c("data.table","mgcv","xgboost","grf"),
                              .export = c("npPSE","BE","BExg","BErf","BWselect","BWselect2","PSEnames","bn")) %dopar% {
      seed=seed*boot.count
      ind=sample(1:nrow(dt), replace=TRUE)
      dt.boot=as.data.frame(dt[ind,])
      if(is.null(MS)){ MS.boot=MS }else{ MS.boot=data.frame(MS[ind,]); names(MS.boot)=names(MS) }
      if(std_data){ dt.boot=as.data.frame(scale(dt.boot)) }
      if(model!="linear"&bwr){
        #pathR.boot=tryCatch({npPSE(dt.boot,exposure,outcome,mediators,covariates,cnfd,z0,z1,intv,MS.boot,model,bw=pathR$bwL)},error = function(e){ return(NULL) })
        pathR.boot=npPSE(dt.boot,exposure,outcome,mediators,covariates,cnfd,z0,z1,intv,MS.boot,model,bw=pathR$bwL)
      }else{
        pathR.boot=npPSE(dt.boot,exposure,outcome,mediators,covariates,cnfd,z0,z1,intv,MS.boot,model,bw=NULL)
      }
      return(list(pseR=pathR.boot$PSE,cfoR=pathR.boot$CFO))
    }
    snow::stopCluster(cl)
    #cat("\n")
    cfo.boot=lapply(var.boot,function(vb){ return(as.list(vb$cfoR)) }) %>% rbindlist(use.names = F)
    var.boot=lapply(var.boot,function(vb){ return(as.list(vb$pseR)) }) %>% rbindlist()
    pathR$CFO=pvbd(stat="cfoc",path=names(pathR$CFO),effect=pathR$CFO,var.boot=cfo.boot,nv=NA)
    pathR$PSE=pvbd(stat="mean",path=names(pathR$PSE),effect=pathR$PSE,var.boot=var.boot,nv=0)
  }

  return(list(CFO=pathR$CFO[,c(-6,-9)],PSE=pathR$PSE))
}

npPSE=function(dt,exposure,outcome,mediators,covariates,cnfd,z0,z1,intv,MS,model,bw){
  Z=dt[[exposure]]
  Y=dt[[outcome]]
  if(length(mediators)>1){
    M=dt[mediators]
  }else{
    M=data.frame(dt[[mediators]])
    colnames(M)=mediators
  }
  if(length(covariates)>1){
    X=dt[covariates]
  }else{
    X=data.frame(dt[[covariates]])
    colnames(X)=covariates
  }
  Xmean=X
  Xmean[names(cnfd)]=cnfd

  #print(paste0("Note: Backward Estimating PSE on ",exposure," > ",paste0(mediators,collapse = " > ")," > ",outcome))
  if(intv=="weak"){
    if(model=="linear"){
      RR=unlist(BE(Z,Y,M,X,Xmean,z0,z1))
    }else if(model=="sim"){
      RR=unlist(BEss(Z,Y,M,X,Xmean,z0,z1))
      # RR=unlist(BEiw(Z,Y,M,X,Xmean,z0,z1,bw),recursive = F)
      # RR=BWselect(RR)
      # nR=names(RR)
      # bwi=substr(nR,nchar(nR),nchar(nR))
      # bwL=RR[bwi=="w"]
      # RR=unlist(RR[bwi!="w"])
    }else if(model=="gam"){
      RR=unlist(BEgm(Z,Y,M,X,Xmean,z0,z1,bw),recursive = F)
      RR=BWselect2(RR)
      nR=names(RR)
      bwi=substr(nR,nchar(nR),nchar(nR))
      bwL=RR[bwi=="w"]
      RR=unlist(RR[bwi!="w"])
    }else if(model=="grf"){
      RR=unlist(BErf(Z,Y,M,X,Xmean,z0,z1))
    }else if(model=="xgb"){
      RR=unlist(BExg(Z,Y,M,X,Xmean,z0,z1))
    }else{
      RR=unlist(BEn(Z,Y,M,X,Xmean,z0,z1))
    }
    ind=2^(0:(length(mediators)+1))
  }else{
    RR=unlist(BEs(Z,Y,M,1:ncol(M),X,Xmean,z0,z1,MS))
    ind=2^(0:(2^length(mediators)))
  }
  ii=nchar(names(RR)[1])-1
  names(RR)=paste0("y",gsub(".y","",substr(names(RR),2,ifelse(substr(names(RR),ii,ii)==".",ii-1,nchar(names(RR))))))
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

  # return(list(CFO=RR[ind],PSE=MD[nnn],bwL=ifelse(!model%in%c("linear","xgb"),list(bwL),NA)[[1]]))
  return(list(CFO=RR[ind],PSE=MD[nnn]))
}

BWselect=function(RR){
  si=sapply(RR,class)
  RRR=unlist(RR[!si%in%c("sibandwidth","numeric")],recursive = F)
  ssi=sapply(RRR,class)
  if(all(ssi%in%c("sibandwidth","numeric"))){
    return(c(RRR,RR[si%in%c("sibandwidth","numeric")]))
  }else{
    return(BWselect(c(RRR,RR[si%in%c("sibandwidth","numeric")])))
  }
}

BWselect2=function(RR){
  si=sapply(RR,class)
  RRR=unlist(RR[si=="list"],recursive = F)
  ssi=sapply(RRR,class)
  if(all(ssi!="list")){
    return(c(RRR,RR[si!="list"]))
  }else{
    return(BWselect2(c(RRR,RR[si!="list"])))
  }
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

BErf=function(Z,Y,M,X,Xmean,z0,z1,bw){
  if(is.null(M)){
    mY=ll_regression_forest(cbind(Z,X), Y, num.trees = 3000)
    MM_z1=predict(mY,newdata=cbind(Z=rep(z1,length(Y)),Xmean)) #Y(za,M1,M2)
    MM_z0=predict(mY,newdata=cbind(Z=rep(z0,length(Y)),Xmean))
  }else{
    mY=ll_regression_forest(cbind(Z,M,X), Y, num.trees = 3000)
    MM_z1=predict(mY,newdata=cbind(Z=rep(z1,length(Y)),M,Xmean)) #Y(za,M1,M2)
    MM_z0=predict(mY,newdata=cbind(Z=rep(z0,length(Y)),M,Xmean))
  }
  if(is.null(M)){
    return(list(y1=MM_z1$predictions[1],y0=MM_z0$predictions[1]))
  }else if(ncol(M)>1){
    return(list(y1=BErf(Z,MM_z1$predictions,M[1:(ncol(M)-1)],X,Xmean,z0,z1),
                y0=BErf(Z,MM_z0$predictions,M[1:(ncol(M)-1)],X,Xmean,z0,z1)))
  }else{
    return(list(y1=BErf(Z,MM_z1$predictions,NULL,X,Xmean,z0,z1,bw=bw[1]),
                y0=BErf(Z,MM_z0$predictions,NULL,X,Xmean,z0,z1,bw=bw[2])))
  }
}

BExg=function(Z,Y,M,X,Xmean,z0,z1){
  if(is.null(M)){
    data=as.matrix(cbind(Z,X))
    data_train=xgb.DMatrix(data=data[1:(nrow(data)*4/5),], label=Y[1:(nrow(data)*4/5)])
    data_valid=xgb.DMatrix(data=data[(nrow(data)*4/5+1):nrow(data),], label=Y[(nrow(data)*4/5+1):nrow(data)])
    mY=xgb.train(
      params = list(eta=0.05,max_depth=4), #booster="gblinear"),#
      data = data_train,
      nrounds = 2000,
      watchlist = list(train=data_train,eval=data_valid),
      early_stopping_rounds = 10,   # stop if no improvement in 10 rounds
      #print_every_n = 500           # show progress every 10 rounds
      verbose = 0
    )
    # resY=Y-predict(mY,newdata=data)
    # MM_z1=predict(mY,newdata=as.matrix(cbind(Z=rep(z1,length(Y)),Xmean))) + resY#Y(za,M1,M2)
    # MM_z0=predict(mY,newdata=as.matrix(cbind(Z=rep(z0,length(Y)),Xmean))) + resY
    MM_z1=predict(mY,newdata=as.matrix(cbind(Z=rep(z1,length(Y)),Xmean))) #Y(za,M1,M2)
    MM_z0=predict(mY,newdata=as.matrix(cbind(Z=rep(z0,length(Y)),Xmean)))
  }else{
    data=as.matrix(cbind(Z,M,X))
    data_train=xgb.DMatrix(data=data[1:(nrow(data)*4/5),], label=Y[1:(nrow(data)*4/5)])
    data_valid=xgb.DMatrix(data=data[(nrow(data)*4/5+1):nrow(data),], label=Y[(nrow(data)*4/5+1):nrow(data)])
    mY=xgb.train(
      params = list(eta=0.05,max_depth=4), #booster="gblinear"),#
      data = data_train,
      nrounds = 2000,
      watchlist = list(train=data_train,eval=data_valid),
      early_stopping_rounds = 10,   # stop if no improvement in 10 rounds
      #print_every_n = 500           # show progress every 10 rounds
      verbose = 0
    )
    # resY=Y-predict(mY,newdata=data)
    # MM_z1=predict(mY,newdata=as.matrix(cbind(Z=rep(z1,length(Y)),M,Xmean))) + resY#Y(za,M1,M2)
    # MM_z0=predict(mY,newdata=as.matrix(cbind(Z=rep(z0,length(Y)),M,Xmean))) + resY
    MM_z1=predict(mY,newdata=as.matrix(cbind(Z=rep(z1,length(Y)),M,Xmean))) #Y(za,M1,M2)
    MM_z0=predict(mY,newdata=as.matrix(cbind(Z=rep(z0,length(Y)),M,Xmean)))
  }
  if(is.null(M)){
    return(list(y1=MM_z1[1],y0=MM_z0[1]))
  }else if(ncol(M)>1){
    return(list(y1=BExg(Z,MM_z1,M[1:(ncol(M)-1)],X,Xmean,z0,z1),
                y0=BExg(Z,MM_z0,M[1:(ncol(M)-1)],X,Xmean,z0,z1)))
  }else{
    return(list(y1=BExg(Z,MM_z1,NULL,X,Xmean,z0,z1),
                y0=BExg(Z,MM_z0,NULL,X,Xmean,z0,z1)))
  }
}

BEgm=function(Z,Y,M,X,Xmean,z0,z1,bw){
  if(is.null(M)){
    fu=as.formula(paste0("Y~Z+",paste0("s(",colnames(X),")",collapse = "+")))
    if(!is.null(bw)){
      mY=gam(fu, data=cbind(Z,X), sp=bw[[length(bw)]])
    }else{
      mY=gam(fu, data=cbind(Z,X))
    }
    MM_z1=predict(mY,newdata=cbind(Z=rep(z1,length(Y)),Xmean)) #Y(za,M1,M2)
    MM_z0=predict(mY,newdata=cbind(Z=rep(z0,length(Y)),Xmean))
  }else{
    fu=as.formula(paste0("Y~Z+",paste0("s(",colnames(M),")",collapse = "+"),"+",paste0("s(",colnames(X),")",collapse = "+")))
    if(!is.null(bw)){
      mY=gam(fu, data=cbind(Z,M,X), sp=bw[[length(bw)]])
    }else{
      mY=gam(fu, data=cbind(Z,M,X), method = "REML")
    }
    MM_z1=predict(mY,newdata=cbind(Z=rep(z1,length(Y)),M,Xmean)) #Y(za,M1,M2)
    MM_z0=predict(mY,newdata=cbind(Z=rep(z0,length(Y)),M,Xmean))
  }
  if(is.null(M)){
    return(list(y1=MM_z1[1],y0=MM_z0[1],bw=mY$sp))
  }else if(ncol(M)>1){
    return(list(y1=BEgm(Z,MM_z1,M[1:(ncol(M)-1)],X,Xmean,z0,z1,bw=bw[1:((length(bw)-1)/2)]),
                y0=BEgm(Z,MM_z0,M[1:(ncol(M)-1)],X,Xmean,z0,z1,bw=bw[((length(bw)-1)/2+1):(length(bw)-1)]),bw=list(mY$sp)))
  }else{
    return(list(y1=BEgm(Z,MM_z1,NULL,X,Xmean,z0,z1,bw=bw[1]),
                y0=BEgm(Z,MM_z0,NULL,X,Xmean,z0,z1,bw=bw[2]),bw=list(mY$sp)))
  }
}

BEss=function(Z,Y,M,X,Xmean,z0,z1){
  if(is.null(M)){
    mY=simsl(Y, Z, X)
    MM_z1=pred.simsl(mY, newA=rep(z1,length(Y)), newX=Xmean)
    MM_z0=pred.simsl(mY, newA=rep(z0,length(Y)), newX=Xmean)
  }else{
    mY=simsl(Y, Z, cbind(M,X))
    MM_z1=pred.simsl(mY, newA=rep(z1,length(Y)), newX=cbind(M,Xmean))
    MM_z0=pred.simsl(mY, newA=rep(z0,length(Y)), newX=cbind(M,Xmean))
  }
  if(is.null(M)){
    return(list(y1=MM_z1$pred.new[1],y0=MM_z0$pred.new[1]))
  }else if(ncol(M)>1){
    return(list(y1=BEss(Z,MM_z1$pred.new,M[1:(ncol(M)-1)],X,Xmean,z0,z1),y0=BEss(Z,MM_z0$pred.new,M[1:(ncol(M)-1)],X,Xmean,z0,z1)))
  }else{
    return(list(y1=BEss(Z,MM_z1$pred.new,NULL,X,Xmean,z0,z1),y0=BEss(Z,MM_z0$pred.new,NULL,X,Xmean,z0,z1)))
  }
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
    return(list(y1=MM_z1[1],y0=MM_z0[1]))
  }else if(ncol(M)>1){
    return(list(y1=BE(Z,MM_z1,M[1:(ncol(M)-1)],X,Xmean,z0,z1),y0=BE(Z,MM_z0,M[1:(ncol(M)-1)],X,Xmean,z0,z1)))
  }else{
    return(list(y1=BE(Z,MM_z1,NULL,X,Xmean,z0,z1),y0=BE(Z,MM_z0,NULL,X,Xmean,z0,z1)))
  }
}

BEiw=function(Z,Y,M,X,Xmean,z0,z1,bw){
  if(is.null(M)){
    fu=as.formula(paste0("Y~Z+",paste0(colnames(X),collapse = "+")))
    # print(fu)
    if(!is.null(bw)){
      mY=npindexbw(fu, data=cbind(Z,X), bws=bw[[length(bw)]], only.optimize.beta=TRUE)
    }else{
      mY=npindexbw(fu, data=cbind(Z,X))
      if(mY$bw>1){
        mY=npindexbw(fu, data=cbind(Z,X),bws=rep(length(Y)^(-1/5),length(all.vars(fu))),only.optimize.beta=T)
        print("bandwidth not converge, only optimize beta with (sample size)^(-1/5)")
      }
    }
    MM_z1=predict(npindex(mY),newdata=cbind(Z=rep(z1,length(Y)),Xmean)) #Y(za,M1,M2)
    MM_z0=predict(npindex(mY),newdata=cbind(Z=rep(z0,length(Y)),Xmean))
  }else{
    fu=as.formula(paste0("Y~Z+",paste0(colnames(M),collapse = "+"),"+",paste0(colnames(X),collapse = "+")))
    # print(fu)
    if(!is.null(bw)){
      mY=npindexbw(fu, data=cbind(Z,M,X), bws=bw[[length(bw)]], only.optimize.beta=TRUE)
    }else{
      mY=npindexbw(fu, data=cbind(Z,M,X))
      if(mY$bw>1){
        mY=npindexbw(fu, data=cbind(Z,M,X),bws=rep(length(Y)^(-1/5),length(all.vars(fu))),only.optimize.beta=T)
        print("bandwidth not converge, only optimize beta with (sample size)^(-1/5)")
      }
    }
    MM_z1=predict(npindex(mY),newdata=cbind(Z=rep(z1,length(Y)),M,Xmean)) #Y(za,M1,M2)
    MM_z0=predict(npindex(mY),newdata=cbind(Z=rep(z0,length(Y)),M,Xmean))
  }
  if(is.null(M)){
    return(list(y1=MM_z1[1],y0=MM_z0[1],bw=mY))
  }else if(ncol(M)>1){
    return(list(y1=BEiw(Z,MM_z1,M[1:(ncol(M)-1)],X,Xmean,z0,z1,bw=bw[1:((length(bw)-1)/2)]),
                y0=BEiw(Z,MM_z0,M[1:(ncol(M)-1)],X,Xmean,z0,z1,bw=bw[((length(bw)-1)/2+1):(length(bw)-1)]),bw=list(mY)))
  }else{
    return(list(y1=BEiw(Z,MM_z1,NULL,X,Xmean,z0,z1,bw=bw[1]),
                y0=BEiw(Z,MM_z0,NULL,X,Xmean,z0,z1,bw=bw[2]),bw=list(mY)))
  }
}

### problematic for three and above mediators (finding sequential estimation order)
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

stTB=function(BB,TV,app,bt="B",pr=F){
  EE=lapply(BB,function(bb){
    return(as.data.table(as.list(setNames(bb[[app]]$effect,bb[[app]]$path))))
  }) %>% rbindlist()
  if(bt=="B"){
    CC=lapply(BB,function(bb){
      return(as.data.table(as.list(setNames(bb[[app]][,4]<TV&TV<bb[[app]][,5],bb[[app]]$path))))
    }) %>% rbindlist()
  }else{
    CC=lapply(BB,function(bb){
      return(as.data.table(as.list(setNames(bb[[app]][,7]<TV&TV<bb[[app]][,8],bb[[app]]$path))))
    }) %>% rbindlist()
  }
  statALL=lapply(1:ncol(EE),ss,EE=EE,CI=CC,TV=TV) %>% rbindlist()
  if(pr==T){
    apply(statALL,1,function(l){
      cat(paste0(paste0(l,collapse = " & ")," \\\\"))
      cat("\n")
    })
  }
  return(statALL)
}

ss=function(c,EE,CI,TV){
  options(scipen = -2)
  PSEhat=EE[[c]]
  return(list(
    effect=names(EE)[c],
    true=signif(TV[c],3),
    MEAN=signif(mean(PSEhat),3),
    RB=signif(mean(PSEhat-TV[c])/TV[c],3),
    SD=signif(sd(PSEhat),3),
    RMSE=signif(sqrt(sum((PSEhat-TV[c])^2)/nrow(EE)),3),
    CP=mean(CI[[c]])
  ))
}

generate_data <- function(sample_size, para, forms, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Step 1: Base covariates
  Z  <- rnorm(sample_size)
  # Z  <- rbinom(sample_size,1,0.4)
  X1 <- rnorm(sample_size)
  X2 <- rnorm(sample_size)
  dt <- data.frame(Z, X1, X2)

  # Step 2: Sequential generation
  for (i in seq_along(forms)) {
    f <- forms[[i]]
    b <- para[[i]]

    lhs <- all.vars(f[[2]])      # e.g., "M2"
    rhs <- all.vars(f[[3]])      # e.g., "Z", "X1", "X2"

    # Ensure RHS variables exist
    missing_vars <- setdiff(rhs, names(dt))
    if (length(missing_vars) > 0) {
      stop(sprintf("Equation %d (%s): missing variables in data: %s",
                   i, lhs, paste(missing_vars, collapse = ", ")))
    }

    # Construct RHS design matrix manually (with intercept)
    rhs_formula <- reformulate(rhs)
    X <- model.matrix(rhs_formula, data = dt)

    # Check coefficient length
    if (length(b) != ncol(X)) {
      stop(sprintf("Equation %d (%s): expected %d coefficients, got %d",
                   i, lhs, ncol(X), length(b)))
    }

    # Simulate response
    y <- as.numeric(X %*% b) + rnorm(sample_size)

    # Add new variable to dataset
    dt[[lhs]] <- y
  }

  return(dt)
}

estimate_true <- function(para, forms, sample_size=1e+6, seed = NULL, std_data=T) {
  if (!is.null(seed)) set.seed(seed)

  # Step 1: Create base predictors
  Z  <- rnorm(sample_size)
  X1 <- rnorm(sample_size)
  X2 <- rnorm(sample_size)
  dt <- data.frame(Z, X1, X2)

  # Storage for results
  results <- vector("list", length(forms))

  # Step 2: Sequentially generate and fit models
  for (i in seq_along(forms)) {
    f <- forms[[i]]
    b <- para[[i]]

    lhs <- all.vars(f[[2]])  # e.g., "M2"
    rhs <- all.vars(f[[3]])  # e.g., "Z", "X1", "X2"

    # Ensure RHS variables exist
    missing_vars <- setdiff(rhs, names(dt))
    if (length(missing_vars) > 0) {
      stop(sprintf("Equation %d (%s): missing vars: %s",
                   i, lhs, paste(missing_vars, collapse = ", ")))
    }

    # Design matrix with intercept
    rhs_formula <- reformulate(rhs)
    X <- model.matrix(rhs_formula, data = dt)

    if (length(b) != ncol(X)) {
      stop(sprintf("Equation %d (%s): expected %d coefficients, got %d",
                   i, lhs, ncol(X), length(b)))
    }

    # Generate outcome
    y <- as.numeric(X %*% b) + rnorm(sample_size)
    dt[[lhs]] <- y

    # Fit lm model
    if(std_data){
      fit <- lm(f, data = as.data.frame(scale(dt)))
    }else{
      fit <- lm(f, data = dt)
    }
    # Store results
    results[[i]]= coef(fit)
  }

  return(results)
}
