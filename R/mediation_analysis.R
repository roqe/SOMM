#' Apply mediation analysis for non-rare binary outcome with two continuous mediators
#'
#' @param dt Input data.
#' @param confounders Confounder values.
#' @param nb Number of bootstrapping. Default is 0 (no bootstrapping applied).
#' @param intv Number of intervention, 3 or 4. Default is 3.
#' @keywords Mediation analysis, Causal inference.
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

mediation_analysis=function(dt,confounders=c(),nb=0,intv=3){
  colnames(dt)[1:4]=c("Y","W","Q","S")
  if(all(dt$Y%in%c(0,1))){
    y.reg=glm(Y~., family = binomial(link="probit"), data=dt)
    total.effect=summary(glm(Y~.-Q-S, family = binomial (link="probit"), data=dt))$coefficients[,1]
  }else{
    y.reg=glm(Y~., family = gaussian, data=dt)
    total.effect=summary(glm(Y~.-Q-S, data=dt))$coefficients[,1]
  }
  s.reg=lm(S~.-Y, data=dt)
  q.reg=lm(Q~.-S-Y, data=dt)
  beta.hat=y.reg$coefficients
  alpha.hat=s.reg$coefficients
  delta.hat=q.reg$coefficients
  ss.hat=sqrt(mean((s.reg$residuals)^2))
  sq.hat=sqrt(mean((q.reg$residuals)^2))
  V.matrix=create_vmatrix(y.reg,s.reg,q.reg)
  if(sum(dt$S,na.rm = T)==0){
    beta.hat[is.na(beta.hat)]=0
    V.matrix[is.na(V.matrix)]=0
  }
  o11=pnorm(sum(total.effect*c(1,1,confounders)))
  o10=pnorm(sum(total.effect*c(1,0,confounders)))
  total.rd=o11-o10 #when S=1 compared to S=0
  total.rr=o11/o10
  total.or=(o11/(1-o11))/(o10/(1-o10))
  bdnp=c("lower(a)","upper(a)","pv(a)")

  if(nb>0){
    var.boot=create_var_boot(nb, dt, confounders=confounders, intv=intv)
    total=dim(var.boot)[1]
    bsRD1=bss(1,F);bsRD2=bss(2,F);bsRD3=bss(3,F)
    bsRR1=bss(4,T);bsRR2=bss(5,T);bsRR3=bss(6,T)
    bsOR1=bss(7,T);bsOR2=bss(8,T);bsOR3=bss(9,T)
    bdnp=c("lower(a)","upper(a)","pv(a)","lower(b)","upper(b)","pv(b)")
  }

  theta_hat=list(beta.hat, alpha.hat, delta.hat, sq.hat, ss.hat)
  if(intv==3){
    p000=omega(theta_hat, c(0,0,0), confounders)
    p100=omega(theta_hat, c(1,0,0), confounders) #first part of difference
    p110=omega(theta_hat, c(1,1,0), confounders) #first part of difference
    p111=omega(theta_hat, c(1,1,1), confounders) #first part of difference
    omega_values=c(p000[1],p100[1],p110[1],p111[1],total.rd,total.rr,total.or)
    names(omega_values)=c("p000","p100","p110","p111","total RD","total RR","total OR")
    RD1=rd(p100,p000,V.matrix)
    RD2=rd(p110,p100,V.matrix)
    RD3=rd(p111,p110,V.matrix)
    RR1=rr(p100,p000,V.matrix)
    RR2=rr(p110,p100,V.matrix)
    RR3=rr(p111,p110,V.matrix)
    OR1=rr(p100/(1-p100),p000/(1-p000),V.matrix)
    OR2=rr(p110/(1-p110),p100/(1-p100),V.matrix)
    OR3=rr(p111/(1-p111),p110/(1-p110),V.matrix)
    if(nb>0){
      pse_values=c(RD1,unlist(bsRD1),RD2,unlist(bsRD2),RD3,unlist(bsRD3),
                   RR1,unlist(bsRR1),RR2,unlist(bsRR2),RR3,unlist(bsRR3),
                   OR1,unlist(bsOR1),OR2,unlist(bsOR2),OR3,unlist(bsOR3))
    }else{
      pse_values=c(RD1,RD2,RD3,RR1,RR2,RR3,OR1,OR2,OR3)
    }
    names(pse_values)=c("RD W>Y",bdnp,"RD W>S>Y",bdnp,"RD W>QY",bdnp,
                        "RR W>Y",bdnp,"RR W>S>Y",bdnp,"RR W>QY",bdnp,
                        "OR W>Y",bdnp,"OR W>S>Y",bdnp,"OR W>QY",bdnp)
  }else if(intv==4){
    p0000=omega(theta_hat, c(0,0,0,0), confounders)
    p1000=omega(theta_hat, c(1,0,0,0), confounders) #first part of difference
    p1100=omega(theta_hat, c(1,1,0,0), confounders) #first part of difference
    p1110=omega(theta_hat, c(1,1,1,0), confounders) #first part of difference
    p1111=omega(theta_hat, c(1,1,1,1), confounders) #first part of difference
    omega_values=c(p0000[1],p1000[1],p1100[1],p1110[1],p1111[1],total.rd,total.rr,total.or)
    names(omega_values)=c("p0000","p1000","p1100","p1110","p1111","total RD","total RR","total OR")
    RD1=rd(p1000,p0000,V.matrix)
    RD2=rd(p1100,p1000,V.matrix)
    RD3=rd(p1110,p1100,V.matrix)
    RD4=rd(p1111,p1110,V.matrix)
    RR1=rr(p1000,p0000,V.matrix)
    RR2=rr(p1100,p1000,V.matrix)
    RR3=rr(p1110,p1100,V.matrix)
    RR4=rr(p1111,p1110,V.matrix)
    OR1=rr(p1000/(1-p1000),p0000/(1-p0000),V.matrix)
    OR2=rr(p1100/(1-p1100),p1000/(1-p1000),V.matrix)
    OR3=rr(p1110/(1-p1110),p1100/(1-p1100),V.matrix)
    OR4=rr(p1111/(1-p1111),p1110/(1-p1110),V.matrix)
    if(nb>0){
      bsRD4=bss(10,F); bsRR4=bss(11,T); bsOR4=bss(12,T)
      pse_values=c(RD1,unlist(bsRD1),RD2,unlist(bsRD2),RD3,unlist(bsRD3),RD4,unlist(bsRD4),
                   RR1,unlist(bsRR1),RR2,unlist(bsRR2),RR3,unlist(bsRR3),RR4,unlist(bsRR4),
                   OR1,unlist(bsOR1),OR2,unlist(bsOR2),OR3,unlist(bsOR3),OR4,unlist(bsOR4))
    }else{
      pse_values=c(RD1,RD2,RD3,RD4,RR1,RR2,RR3,RR4,OR1,OR2,OR3,OR4)
    }
    names(pse_values)=c("RD W>Y",bdnp,"RD W>S>Y",bdnp,"RD W>Q>Y",bdnp,"RD W>Q>S>Y",bdnp,
                        "RR W>Y",bdnp,"RR W>S>Y",bdnp,"RR W>Q>Y",bdnp,"RR W>Q>S>Y",bdnp,
                        "OR W>Y",bdnp,"OR W>S>Y",bdnp,"OR W>Q>Y",bdnp,"OR W>Q>S>Y",bdnp)
  }
  return(c(omega_values,pse_values))
}
