#' Apply mediation analysis for non-rare binary outcome with two continuous mediators
#'
#' @param dt Input data.
#' @param confounders Confounder values.
#' @param nb Number of bootstrapping. Default is 0 (no bootstrapping applied).
#' @param intv Number of intervention, 3 or 4. Default is 3.
#' @keywords Mediation analysis, Causal inference.
#' @import data.table lme4 foreach snow doSNOW
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

mediation_analysis=function(dt,confounders=c(),nb=0,intv=3,unit=1,reNAME=NULL,mc=5,autoR=T){
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
  bdnp=c("lower(a)","upper(a)","pv(a)")
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
    var.boot=data.table::rbindlist(foreach::foreach(boot.count=1:nb, .options.snow = opts) %dopar% {
      create_var_boot(nb, dt, confounders=confounders, intv=4, reNAME=reNAME, x0, x1) })

    bsRD1=bss(1,var.boot,F);bsRD2=bss(2,var.boot,F);bsRD3=bss(3,var.boot,F)
    bsRR1=bss(4,var.boot,T);bsRR2=bss(5,var.boot,T);bsRR3=bss(6,var.boot,T)
    bsOR1=bss(7,var.boot,T);bsOR2=bss(8,var.boot,T);bsOR3=bss(9,var.boot,T)
    bdnp=c(bdnp,"lower(b)","upper(b)","pv(b)")
  }

  if(intv==3){
    PP=PSE_three(GT,x0,x1,confounders,V.matrix)
    psel=c("W>Y","W>S>Y","W>QY")
    pse_values=data.table::data.table(stat=c(rep("RD",3),rep("RR",3),rep("OR",3)),
                                      path=rep(psel,3),rbind(PP$RD,PP$RR,PP$OR))
    if(nb>0){
      pse_values=cbind(pse_values,
        rbind(unlist(bsRD1),unlist(bsRD2),unlist(bsRD3),
              unlist(bsRR1),unlist(bsRR2),unlist(bsRR3),
              unlist(bsOR1),unlist(bsOR2),unlist(bsOR3)))
    }
  }else if(intv==4){
    PP=PSE_four(GT,x0,x1,confounders,V.matrix)
    psel=c("W>Y","W>S>Y","W>Q>Y","W>Q>S>Y")
    pse_values=data.table::data.table(stat=c(rep("RD",4),rep("RR",4),rep("OR",4)),
                                      path=rep(psel,3),rbind(PP$RD,PP$RR,PP$OR))
    if(nb>0){
      bsRD4=bss(10,var.boot,F); bsRR4=bss(11,var.boot,T); bsOR4=bss(12,var.boot,T)
      pse_values=cbind(pse_values,
        rbind(unlist(bsRD1),unlist(bsRD2),unlist(bsRD3),unlist(bsRD4),
              unlist(bsRR1),unlist(bsRR2),unlist(bsRR3),unlist(bsRR4),
              unlist(bsOR1),unlist(bsOR2),unlist(bsOR3),unlist(bsOR4)))
    }
  }
  colnames(pse_values)[3:ncol(pse_values)]=c("effect",bdnp)
  if(autoR){ pse_values=pse_values[complete.cases(pse_values)] }
  return(list(DAG=PP$omega,TOTAL=data.table(RD=o11-o10,RR=o11/o10,OR=(o11/(1-o11))/(o10/(1-o10))),PSE=pse_values))
}
