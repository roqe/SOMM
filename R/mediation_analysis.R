#' Apply mediation analysis for non-rare binary outcome with two continuous mediators
#'
#' @param dt Input data.
#' @param confounders Confounder values.
#' @param nb Number of bootstrapping. Default is 0 (no bootstrapping applied).
#' @param intv Number of intervention, 3 or 4. Default is 3.
#' @keywords Mediation analysis, Causal inference.
#' @import data.table lme4 foreach snow doSNOW doParallel
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
#' res3_3=mediation_analysis(eIVF, confounders=c(log(36),log(26)))
#' res3_4=mediation_analysis(eIVF, confounders=c(log(36),log(26)), intv=4)
#' para_sg=c(rep(-0.5,3),rep(0,4),rep(-0.5,2),1,1) # for single mediator, set exposure-related parameters into 0
#' dat_sg=sim_mediation_data(0.5,1000,para_sg)
#' dat_sg$S=0 # set the second mediator into 0
#' res_sg_3=mediation_analysis(dat_sg)

mediation_analysis=function(dt,confounders=c(),nb=0,intv=3,unit=1,reNAME=NULL,mc=5,autoR=T,stra=1,seed=217){
  colnames(dt)[1:4]=c("Y","W","Q","S")
  colnames(dt)[colnames(dt)==reNAME]="id"
  if(unit=="IQR"){
    QQ=quantile(dt$W,na.rm=T);
    x0=QQ[[2]]; x1=QQ[[4]]
  }else{
    x0=0; x1=1
  }
  GT=get_theta(dt,reNAME)
  o11=pnorm(sum(GT$total*c(1,x1,confounders)))
  o10=pnorm(sum(GT$total*c(1,x0,confounders)))
  nnaa=c("lower(a)","upper(a)","pv(a)")
  nnbb=c("lower(b)","upper(b)","pv(b)")
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
      set.seed(seed*boot.count)
      create_var_boot(nb, dt, confounders=confounders, intv=intv, reNAME=reNAME, x0, x1, stra) })
    stopCluster(cl)

    set.seed(seed)
    bsRD1=bss(1,var.boot,F);bsRD2=bss(2,var.boot,F);bsRD3=bss(3,var.boot,F);bsRDT=bss(13,var.boot,F)
    bsRR1=bss(4,var.boot,T);bsRR2=bss(5,var.boot,T);bsRR3=bss(6,var.boot,T);bsRRT=bss(14,var.boot,T)
    bsOR1=bss(7,var.boot,T);bsOR2=bss(8,var.boot,T);bsOR3=bss(9,var.boot,T);bsORT=bss(15,var.boot,T)
    pmRD1=bss(16,var.boot,F,T);pmRD2=bss(17,var.boot,F,T);pmRD3=bss(18,var.boot,F,T)
    pmRR1=bss(19,var.boot,F,T);pmRR2=bss(20,var.boot,F,T);pmRR3=bss(21,var.boot,F,T)
    pmOR1=bss(22,var.boot,F,T);pmOR2=bss(23,var.boot,F,T);pmOR3=bss(24,var.boot,F,T)
  }

  if(intv==3){
    PP=PSE_three(GT,x0,x1,confounders,V.matrix)
    psel=c("W>Y","W>S>Y","W>QY","total")
    pse_values=data.table::data.table(stat=c(rep("RD",4),rep("RR",4),rep("OR",4)),path=rep(psel,3),rbind(PP$RD,PP$RR,PP$OR))
    colnames(pse_values)[3:ncol(pse_values)]=c("effect",nnaa,"pm_effect")
    pme_values=pse_values[,c(1,2,7)]
    if(nb>0){
      pse_values=cbind(pse_values[,1:6],
        rbind(btbd(pse_values$effect[1],bsRD1,F),btbd(pse_values$effect[2],bsRD2,F),
              btbd(pse_values$effect[3],bsRD3,F),btbd(pse_values$effect[4],bsRDT,F),
              btbd(pse_values$effect[5],bsRR1,T),btbd(pse_values$effect[6],bsRR2,T),
              btbd(pse_values$effect[7],bsRR3,T),btbd(pse_values$effect[8],bsRRT,T),
              btbd(pse_values$effect[9],bsOR1,T),btbd(pse_values$effect[10],bsOR2,T),
              btbd(pse_values$effect[11],bsOR3,T),btbd(pse_values$effect[12],bsORT,T)))
      colnames(pse_values)[7:ncol(pse_values)]=nnbb
      pme_values=cbind(pme_values,
        rbind(btbd(pme_values$pm_effect[1],pmRD1,F),btbd(pme_values$pm_effect[2],pmRD2,F),
              btbd(pme_values$pm_effect[3],pmRD3,F),rep(NA,3),
              btbd(pme_values$pm_effect[5],pmRR1,F),btbd(pme_values$pm_effect[6],pmRR2,F),
              btbd(pme_values$pm_effect[7],pmRR3,F),rep(NA,3),
              btbd(pme_values$pm_effect[9],pmOR1,F),btbd(pme_values$pm_effect[10],pmOR2,F),
              btbd(pme_values$pm_effect[11],pmOR3,F),rep(NA,3)))
      colnames(pme_values)[4:ncol(pme_values)]=nnbb
    }else{
      pse_values=pse_values[,1:6]
    }
  }else if(intv==4){
    PP=PSE_four(GT,x0,x1,confounders,V.matrix)
    psel=c("W>Y","W>S>Y","W>Q>Y","W>Q>S>Y","total")
    pse_values=data.table::data.table(stat=c(rep("RD",5),rep("RR",5),rep("OR",5)),path=rep(psel,3),rbind(PP$RD,PP$RR,PP$OR))
    colnames(pse_values)[3:ncol(pse_values)]=c("effect",nnaa,"pm_effect")
    pme_values=pse_values[,c(1,2,7)]
    if(nb>0){
      bsRD4=bss(10,var.boot,F);bsRR4=bss(11,var.boot,T);bsOR4=bss(12,var.boot,T)
      pmRD4=bss(25,var.boot,F,T);pmRR4=bss(26,var.boot,F,T);pmOR4=bss(27,var.boot,F,T)
      pse_values=cbind(pse_values[,1:6],
                       rbind(btbd(pse_values$effect[1],bsRD1,F),btbd(pse_values$effect[2],bsRD2,F),
                             btbd(pse_values$effect[3],bsRD3,F),btbd(pse_values$effect[4],bsRD4,F),btbd(pse_values$effect[5],bsRDT,F),
                             btbd(pse_values$effect[6],bsRR1,T),btbd(pse_values$effect[7],bsRR2,T),
                             btbd(pse_values$effect[8],bsRR3,T),btbd(pse_values$effect[9],bsRR4,T),btbd(pse_values$effect[10],bsRRT,T),
                             btbd(pse_values$effect[11],bsOR1,T),btbd(pse_values$effect[12],bsOR2,T),
                             btbd(pse_values$effect[13],bsOR3,T),btbd(pse_values$effect[14],bsOR4,T),btbd(pse_values$effect[15],bsORT,T)))
      colnames(pse_values)[7:ncol(pse_values)]=nnbb
      pme_values=cbind(pme_values,
                       rbind(btbd(pme_values$pm_effect[1],pmRD1,F),btbd(pme_values$pm_effect[2],pmRD2,F),
                             btbd(pme_values$pm_effect[3],pmRD3,F),btbd(pme_values$pm_effect[4],pmRD4,F),rep(NA,3),
                             btbd(pme_values$pm_effect[6],pmRR1,F),btbd(pme_values$pm_effect[7],pmRR2,F),
                             btbd(pme_values$pm_effect[8],pmRR3,F),btbd(pme_values$pm_effect[9],pmRR4,F),rep(NA,3),
                             btbd(pme_values$pm_effect[11],pmOR1,F),btbd(pme_values$pm_effect[12],pmOR2,F),
                             btbd(pme_values$pm_effect[13],pmOR3,F),btbd(pme_values$pm_effect[14],pmOR4,F),rep(NA,3)))
      colnames(pme_values)[4:ncol(pme_values)]=nnbb
    }else{
      pse_values=pse_values[,1:6]
    }
  }
  if(autoR){
    pse_values=pse_values[!is.na(pse_values$`pv(a)`),]
    pme_values=pme_values[pme_values$pm_effect!=0,]
    pme_values$`lower(b)`[pme_values$`lower(b)`<0]=0
    pme_values$`upper(b)`[pme_values$`upper(b)`>1]=1
  }
  return(list(TOTAL=data.table(RD=o11-o10,RR=o11/o10,OR=(o11/(1-o11))/(o10/(1-o10))),DAG=PP$omega,PSE=pse_values,PM=pme_values))
}
