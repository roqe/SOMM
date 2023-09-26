#' @import data.table lme4 survival
#' @export
get_theta=function(dt,dt2,dt3,reNAME,grpID){
  if(!is.null(reNAME)){
    print("Note: Apply random-effct models.")
    covnames=colnames(dt)[!colnames(dt)%in%c("Y","W","Q","S","id","grp")]
    if(length(covnames)==0){ nnn=NULL }else{ nnn=paste0("+",paste0(covnames,collapse = "+")) }
    if(all(dt$Y%in%c(0,1))){
      y.reg=lme4::glmer(as.formula(paste0("Y~W+Q+S+(1|id)",nnn)), family=binomial(link="probit"), data=dt)
      total.effect=summary(lme4::glmer(as.formula(paste0("Y~W+(1|id)",nnn)), family = binomial(link="probit"), data=dt))$coef[,1]
    }else{
      y.reg=lme4::lmer(as.formula(paste0("Y~W+Q+S+(1|id)",nnn)), data=dt)
      total.effect=summary(lme4::lmer(as.formula(paste0("Y~W+(1|id)",nnn)), data=dt))$coef[,1]
    }
    s.reg=lme4::lmer(as.formula(paste0("S~W+Q+(1|id)",nnn)), data=dt)
    q.reg=lme4::lmer(as.formula(paste0("Q~W+(1|id)",nnn)), data=dt)
  }else if(!is.null(grpID)){
    if(all(dt$Y%in%c(0,1))){
      print("Note: Apply conditional logistic models.")
      total.effect=NA
      y.reg=survival::clogit(Y~.+strata(grp), data=dt)
      if(!is.null(dt2) & all(dt2$S%in%c(0,1))){
        s.reg=survival::clogit(S~.+strata(grp), data=dt2)
      }else{ s.reg=NA }
      if(!is.null(dt3) & all(dt3$Q%in%c(0,1))){
        q.reg=survival::clogit(Q~.+strata(grp), data=dt3)
      }else{ q.reg=NA }
    }else{
      print("Warn: Only supporting conditional logistic regression.")
    }
  }else{
    if(all(dt$Y%in%c(0,1))){
      print("Note: Apply GLM models, binomial outcome, probit link.")
      y.reg=glm(Y~., family = binomial(link="probit"), data=dt)
      total.effect=summary(glm(Y~.-Q-S, family = binomial(link="probit"), data=dt))$coef[,1]
    }else{
      print("Note: Apply GLM models, gaussian outcome.")
      y.reg=glm(Y~., family = gaussian, data=dt)
      total.effect=summary(glm(Y~.-Q-S, data=dt))$coef[,1]
    }
    s.reg=lm(S~.-Y, data=dt)
    q.reg=lm(Q~.-S-Y, data=dt)
  }

  yc=data.table::data.table(t(summary(y.reg)$coef[,1]))
  covnames=names(yc)[which(!names(yc)%in%c("(Intercept)","Y","W","Q","S","id","grp"))]
  beta.hat=data.table::data.table(I=ifelse("(Intercept)"%in%names(yc),yc[["(Intercept)"]],0),yc[1,c("W","Q")],
                                  S=ifelse("S"%in%names(yc),yc[["S"]],0),yc[1,..covnames])
  data.table::setnafill(beta.hat,fill = 0)
  if(ifelse("S"%in%names(yc),!is.na(yc[["S"]]),F)){
    alpha.hat=summary(s.reg)$coef[,1]
  }else{
    alpha.hat=beta.hat[,!"S"]
    alpha.hat[1,]=0
  }
  dc=data.table::data.table(t(summary(q.reg)$coef[,1]))
  delta.hat=data.table::data.table(I=ifelse("(Intercept)"%in%names(dc),dc[["(Intercept)"]],0),dc[1,c("W")],dc[1,..covnames])
  ss.hat=ifelse("S"%in%names(yc),ifelse(!is.na(yc[["S"]]),sqrt(mean(resid(s.reg)^2)),0),0)
  sq.hat=sqrt(mean(resid(q.reg)^2))
  theta_hat=list(bc=unlist(beta.hat), ac=unlist(alpha.hat), dc=unlist(delta.hat), sq=sq.hat, ss=ss.hat)
  return(list(theta_hat=theta_hat,reg=list(y.reg=y.reg,s.reg=s.reg,q.reg=q.reg),total=total.effect))
}
