#' @import data.table
bss=function(ind,var.boot,log_scale,pos=F,exact=F){
  dd=var.boot[[ind]]
  if(pos&(sd(dd)!=0)){ dd=dd[dd>0] }
  if(log_scale){ dd=log(dd) }
  if(exact|sd(dd)==0){
    cc=length(dd)
  }else{
    cc=1e+6; dd=rnorm(cc, mean=mean(dd), sd=sd(dd))
  }
  vv=sum(dd>0)/cc
  pp=min(vv,(1-vv))*2
  return(list(pp=pp,btsd=sd(dd)))
}

btbd=function(trmn,bss,log_scale){
  if(log_scale){
    return(c(exp(log(trmn)-1.96*bss$btsd),exp(log(trmn)+1.96*bss$btsd),bss$pp))
  }else{
    return(c(trmn-1.96*bss$btsd,trmn+1.96*bss$btsd,bss$pp))
  }
}
