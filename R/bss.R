#' @import data.table
bss=function(ind,var.boot,log_scale,pos=F,exact=F){
  dd=var.boot[[ind]]
  if(all(is.na(dd))){ return(list(pp=NA,btsd=NA)) }
  sddd=sd(dd,na.rm = T)
  if(pos&sddd!=0){ dd=dd[dd>=0] }
  if(log_scale){ dd=log(dd); sddd=sd(dd,na.rm = T) }
  if(exact|sddd==0){
    cc=length(dd)
  }else{
    cc=1e+6; dd=rnorm(cc, mean=mean(dd,na.rm = T), sd=sddd)
  }
  vv=sum(dd>0)/cc
  pp=min(vv,(1-vv))*2
  return(list(pp=pp,btsd=sddd))
}

btbd=function(trmn,bss,log_scale){
  if(log_scale){
    return(c(exp(log(trmn)-1.96*bss$btsd),exp(log(trmn)+1.96*bss$btsd),bss$pp))
  }else{
    return(c(trmn-1.96*bss$btsd,trmn+1.96*bss$btsd,bss$pp))
  }
}
