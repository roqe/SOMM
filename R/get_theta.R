#' @import data.table lme4
#' @export
get_theta=function(dt,reNAME){
  covnames=colnames(dt)[!colnames(dt)%in%c("Y","W","Q","S","id")]
  if(is.null(reNAME)){
    if(all(dt$Y%in%c(0,1))){
      y.reg=glm(Y~., family = binomial(link="probit"), data=dt)
      total.effect=summary(glm(Y~.-Q-S, family = binomial(link="probit"), data=dt))$coef[,1]
    }else{
      y.reg=glm(Y~., family = gaussian, data=dt)
      total.effect=summary(glm(Y~.-Q-S, data=dt))$coef[,1]
    }
    s.reg=lm(S~.-Y, data=dt)
    q.reg=lm(Q~.-S-Y, data=dt)
  }else{
    if(all(dt$Y%in%c(0,1))){
      y.reg=lme4::glmer(as.formula(paste0("Y~W+Q+S+(1|id)+",paste0(covnames,collapse = "+"))), family=binomial(link="probit"), data=dt)
      total.effect=summary(lme4::glmer(as.formula(paste0("Y~W+(1|id)+",paste0(covnames,collapse = "+"))), family = binomial(link="probit"), data=dt))$coef[,1]
    }else{
      y.reg=lme4::lmer(as.formula(paste0("Y~W+Q+S+(1|id)+",paste0(covnames,collapse = "+"))), data=dt)
      total.effect=summary(lme4::lmer(as.formula(paste0("Y~W+(1|id)+",paste0(covnames,collapse = "+"))), data=dt))$coef[,1]
    }
    s.reg=lme4::lmer(as.formula(paste0("S~W+Q+(1|id)+",paste0(covnames,collapse = "+"))), data=dt)
    q.reg=lme4::lmer(as.formula(paste0("Q~W+(1|id)+",paste0(covnames,collapse = "+"))), data=dt)
  }

  yc=data.table::data.table(t(summary(y.reg)$coef[,1]))
  covnames=names(yc)[which(!names(yc)%in%c("(Intercept)","Y","W","Q","S","id"))]
  beta.hat=cbind(data.table::data.table(yc[1,1:3],S=ifelse("S"%in%names(yc),yc[["S"]],0),yc[1,..covnames]))
  if("S"%in%names(yc)){
    alpha.hat=summary(s.reg)$coef[,1]
  }else{
    yc[1,]=0; alpha.hat=yc[1,]
  }
  delta.hat=summary(q.reg)$coef[,1]
  ss.hat=ifelse("S"%in%names(yc),sqrt(mean((summary(s.reg)$residuals)^2)),0)
  sq.hat=sqrt(mean((summary(q.reg)$residuals)^2))
  theta_hat=list(bc=unlist(beta.hat), ac=unlist(alpha.hat), dc=unlist(delta.hat), sq=sq.hat, ss=ss.hat)
  return(list(theta_hat=theta_hat,reg=list(y.reg=y.reg,s.reg=s.reg,q.reg=q.reg),total=total.effect))
}
