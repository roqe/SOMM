#' Apply mediation analysis for non-rare binary outcome with two continuous mediators
#'
#' @param dt Input data.
#' @param cnfd Confounder values.
#' @param nb Number of bootstrapping. Default is 0 (no bootstrapping applied).
#' @param intv Number of intervention, 3 or 4. Default is 3.
#' @keywords Mediation analysis, Causal inference.
#' @import data.table lme4 foreach snow doSNOW doParallel survival
#' @export
#' @examples
#' para=c(rep(-0.5,9),1,1)
#' dat1=sim_mediation_data(0.5,1000,para) #binary exposure
#' apply(dat1,2,mean) #the proportion of Y should not be too skew (nearly 0 or 1)
#' res1_3=mediation_analysis(dat1)
#' tru1_3=calc_true_value(para)
#' dat2=sim_mediation_data(c(0,1),1000,para) #continous exposure
#' res2_4=mediation_analysis(dat2, intv=4, nb=500)
#' tru2_4=calc_true_value(para, intv=4)
#' summary(eIVF) #subset eIVF data for demo
#' res3_3=mediation_analysis(eIVF, cnfd=c(age=log(36),bmi=log(26)))
#' res3_4=mediation_analysis(eIVF, cnfd=c(age=log(36),bmi=log(26)), intv=4)
#' para_sg=c(rep(-0.5,3),rep(0,4),rep(-0.5,2),1,1) # for single mediator, set exposure-related parameters into 0
#' dat_sg=sim_mediation_data(0.5,1000,para_sg)
#' dat_sg$S=0 # set the second mediator into 0
#' res_sg_3=mediation_analysis(dat_sg)

mediation_analysis=function(dt,cnfd=c(),dt2=NULL,cnfd2=c(),dt3=NULL,cnfd3=c(),nb=0,intv=4,unit=1,reNAME=NULL,grpID=NULL,mc=5,autoR=T,stra=1,seed=217){
  if(!is.null(cnfd)){
    if(!all(names(cnfd)%in%names(dt)) | is.null(names(cnfd)))
      return("Colnames of covariates should a subset of colnames of dataset.")
  }
  print(paste0("Note: Estimating PSE on ",colnames(dt)[2],"(W) > ",colnames(dt)[3],"(Q) > ",colnames(dt)[4],"(S) > ",colnames(dt)[1],"(Y)"))
  colnames(dt)[1:4]=c("Y","W","Q","S")
  if(!is.null(reNAME)){
    print("Warn: The asymptotic CI/p-values using random effect models may be inaccurate.")
    if(nb==0){ nb=500; print("Note: Applying bootstrapping 500 times.") }
    colnames(dt)[colnames(dt)==reNAME]="id"
  }
  if(!is.null(grpID)){
    print("Warn: The asymptotic CI/p-values using conditional logistic models may be inaccurate.")
    if(nb==0){ nb=500; print("Note: Applying bootstrapping 500 times.") }
    colnames(dt)[colnames(dt)==grpID]="grp"
  }
  if(!is.null(dt2)){
    print(paste0("Note: Using another population for estimating ",colnames(dt2)[1],"(W) > ",colnames(dt2)[2],"(Q) > ",colnames(dt2)[3],"(S)"))
    colnames(dt2)[1:3]=c("W","Q","S")
    colnames(dt2)[colnames(dt2)==grpID]="grp"
  }
  if(!is.null(dt3)){
    print(paste0("Note: Using another population for estimating ",colnames(dt3)[1],"(W) > ",colnames(dt3)[2],"(Q)"))
    colnames(dt3)[1:2]=c("W","Q")
    colnames(dt3)[colnames(dt3)==grpID]="grp"
  }
  if(unit=="IQR"){
    QQ=quantile(dt$W,na.rm=T);
    x0=QQ[[2]]; x1=QQ[[4]]
  }else{
    x0=0; x1=1
  }
  GT=get_theta(dt,dt2,dt3,reNAME,grpID)
  BF=ifelse(all(dt$Y%in%c(0,1)),T,F)
  fid=unique(c(names(dt)[which(sapply(dt,is.factor))],names(dt)[which(sapply(dt,is.character))]))
  eid1=names(cnfd)%in%c(fid,names(GT$theta_hat$bc))
  eid2=names(cnfd2)%in%c(fid,names(GT$theta_hat$ac))
  eid3=names(cnfd3)%in%c(fid,names(GT$theta_hat$dc))
  if(sum(!eid1)!=0){ print(paste("model Y remove covariate:",paste(names(cnfd)[!eid1],collapse = ", "))) }
  if(length(cnfd2)!=0&sum(!eid2)!=0){ print(paste("model S remove covariate:",paste(names(cnfd2)[!eid2],collapse = ", "))) }
  if(length(cnfd3)!=0&sum(!eid3)!=0){ print(paste("model Q remove covariate:",paste(names(cnfd3)[!eid3],collapse = ", "))) }
  cnfd=cnfd[eid1]; cnfd2=cnfd2[eid2]; cnfd3=cnfd3[eid3]
  if(BF){
    o11=pnorm(sum(GT$total[c("(Intercept)","W",names(cnfd))]*c(1,x1,cnfd)))
    o10=pnorm(sum(GT$total[c("(Intercept)","W",names(cnfd))]*c(1,x0,cnfd)))
  }else{
    o11=sum(GT$total[c("(Intercept)","W",names(cnfd))]*c(1,x1,cnfd))/GT$theta_hat$st
    o10=sum(GT$total[c("(Intercept)","W",names(cnfd))]*c(1,x0,cnfd))/GT$theta_hat$st
  }
  nnaa=c("lower(a)","upper(a)","pv(a)")
  nnbb=c("lower(b)","upper(b)","pv(b)")
  nnba=c("lower(n)","upper(n)","pv(n)")
  V.matrix=create_vmatrix(GT)
  if(sum(dt$S,na.rm = T)==0){
    GT$theta_hat$bc[is.na(GT$theta_hat$bc)]=0
    V.matrix[is.na(V.matrix)]=0
  }

  if(nb>0){
    cl = snow::makeCluster(mc)
    doSNOW::registerDoSNOW(cl)
    pb = txtProgressBar(max = nb, style = 3)
    progress = function(n) setTxtProgressBar(pb, n)
    opts = list(progress = progress)
    var.boot=data.table::rbindlist(foreach::foreach(boot.count=1:nb, .options.snow = opts, .packages = "data.table") %dopar% {
      create_var_boot(dt,cnfd,dt2,cnfd2,dt3,cnfd3,nb,intv,reNAME,grpID,x0,x1,BF,stra,seed=seed*boot.count) })
    snow::stopCluster(cl)

    set.seed(seed)
    bsRD1=bss( 1,var.boot,F,F);bsRD2=bss( 2,var.boot,F,F);bsRD3=bss( 3,var.boot,F,F);bsRDT=bss(13,var.boot,F,F)
    pmRD1=bss(16,var.boot,F,T);pmRD2=bss(17,var.boot,F,T);pmRD3=bss(18,var.boot,F,T)
    if(BF){
      bsRR1=bss( 4,var.boot,T,T);bsRR2=bss( 5,var.boot,T,T);bsRR3=bss( 6,var.boot,T,T);bsRRT=bss(14,var.boot,T,T)
      bsOR1=bss( 7,var.boot,T,T);bsOR2=bss( 8,var.boot,T,T);bsOR3=bss( 9,var.boot,T,T);bsORT=bss(15,var.boot,T,T)
      pmRR1=bss(19,var.boot,F,T);pmRR2=bss(20,var.boot,F,T);pmRR3=bss(21,var.boot,F,T)
      pmOR1=bss(22,var.boot,F,T);pmOR2=bss(23,var.boot,F,T);pmOR3=bss(24,var.boot,F,T)
    }
  }

  if(intv==3){
    PP=PSE_three(GT,x0,x1,cnfd,cnfd2,cnfd3,V.matrix,BF)
    psel=c("W>Y","W>S>Y","W>QY","total")
    if(BF){
      pse_values=data.table::data.table(stat=c(rep("RD",4),rep("RR",4),rep("OR",4)),path=rep(psel,3),rbind(PP$RD,PP$RR,PP$OR))
    }else{
      pse_values=data.table::data.table(stat=c(rep("mean",4)),path=psel,PP$RD)
    }
    colnames(pse_values)[3:ncol(pse_values)]=c("effect",nnaa,"pm_effect")
    pme_values=pse_values[,c(1,2,7)]
    if(nb>0){
      if(BF){
        pse_values=cbind(pse_values[,1:6],rbind(bsRD1,bsRD2,bsRD3,bsRDT,bsRR1,bsRR2,bsRR3,bsRRT,bsOR1,bsOR2,bsOR3,bsORT))
      }else{
        pse_values=cbind(pse_values[,1:6],rbind(bsRD1,bsRD2,bsRD3,bsRDT))
      }
      colnames(pse_values)[7:ncol(pse_values)]=c(nnbb,nnba)
      if(abs(sum(sign(pse_values$effect[1:3])))==sum(pse_values$effect[1:3]!=0)){
        if(BF){
          pme_values=cbind(pme_values,rbind(pmRD1,pmRD2,pmRD3,rep(NA,3),pmRR1,pmRR2,pmRR3,rep(NA,3),pmOR1,pmOR2,pmOR3,rep(NA,3)))
        }else{
          pme_values=cbind(pme_values,rbind(pmRD1,pmRD2,pmRD3,rep(NA,3)))
        }
        colnames(pme_values)[4:ncol(pme_values)]=c(nnbb,nnba)
      }
    }else{
      pse_values=pse_values[,1:6]
    }
  }else if(intv==4){
    PP=PSE_four(GT,x0,x1,cnfd,cnfd2,cnfd3,V.matrix,BF)
    psel=c("W>Y","W>S>Y","W>Q>Y","W>Q>S>Y","total")
    if(BF){
      pse_values=data.table::data.table(stat=c(rep("RD",5),rep("RR",5),rep("OR",5)),path=rep(psel,3),rbind(PP$RD,PP$RR,PP$OR))
    }else{
      pse_values=data.table::data.table(stat=c(rep("mean",5)),path=psel,PP$RD)
    }
    colnames(pse_values)[3:ncol(pse_values)]=c("effect",nnaa,"pm_effect")
    pme_values=pse_values[,c(1,2,7)]
    if(nb>0){
      bsRD4=bss(10,var.boot,F,F);bsRR4=bss(11,var.boot,T,F);bsOR4=bss(12,var.boot,T,F)
      pmRD4=bss(25,var.boot,F,T);pmRR4=bss(26,var.boot,F,T);pmOR4=bss(27,var.boot,F,T)
      if(BF){
        pse_values=cbind(pse_values[,1:6],rbind(bsRD1,bsRD2,bsRD3,bsRD4,bsRDT,bsRR1,bsRR2,bsRR3,bsRR4,bsRRT,bsOR1,bsOR2,bsOR3,bsOR4,bsORT))
      }else{
        pse_values=cbind(pse_values[,1:6],rbind(bsRD1,bsRD2,bsRD3,bsRD4,bsRDT))
      }
      colnames(pse_values)[7:ncol(pse_values)]=c(nnbb,nnba)
      if(abs(sum(sign(pse_values$effect[1:4])))==sum(pse_values$effect[1:4]!=0)){
        if(BF){
          pme_values=cbind(pme_values,rbind(pmRD1,pmRD2,pmRD3,pmRD4,rep(NA,3),pmRR1,pmRR2,pmRR3,pmRR4,rep(NA,3),pmOR1,pmOR2,pmOR3,pmOR4,rep(NA,3)))
        }else{
          pme_values=cbind(pme_values,rbind(pmRD1,pmRD2,pmRD3,pmRD4,rep(NA,3)))
        }
        colnames(pme_values)[4:ncol(pme_values)]=c(nnbb,nnba)
      }
    }else{
      pse_values=pse_values[,1:6]
    }
  }
  if(autoR){
    pse_values=pse_values[!is.na(pse_values$`pv(a)`),]
    pme_values=pme_values[pme_values$pm_effect!=0,]
    pme_values$`lower(b)`[pme_values$`lower(b)`<0]=0
    pme_values$`upper(b)`[pme_values$`upper(b)`>1]=1
    pme_values$`lower(n)`[pme_values$`lower(n)`<0]=0
    pme_values$`upper(n)`[pme_values$`upper(n)`>1]=1
  }
  if(BF){
    return(list(TOTAL=data.table(crude_RD=o11-o10,crude_RR=o11/o10,crude_OR=(o11/(1-o11))/(o10/(1-o10))),
                OMEGA=PP$omega,PSE=pse_values,PM=pme_values))
  }else{
    return(list(TOTAL=data.table(crude_mean_diff=o11-o10),OMEGA=PP$omega,PSE=pse_values,PM=pme_values))
  }
}
