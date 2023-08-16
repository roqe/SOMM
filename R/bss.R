#' @import data.table
bss=function(ind,var.boot,log_scale=T,pos=F,exact=F){
  dd=var.boot[[ind]]
  if(pos&(sd(dd)!=0)){ dd=dd[dd>0] }
  cc=length(dd)
  bb=c(quantile(dd, probs=c(0.025, 0.975), na.rm=T))
  if(log_scale){
    if((!exact)&(sd(dd)!=0)){
      cc=1e+6; dd=rnorm(cc, mean=mean(log(dd)), sd=sd(log(dd))); vv=sum(dd>0)/cc
    }else{ vv=sum(log(dd)>0)/cc }
  }else{
    if((!exact)&(sd(dd)!=0)){
      cc=1e+6; dd=rnorm(cc, mean=mean(dd), sd=sd(dd))
    }
    vv=sum(dd>0)/cc
  }
  pp=min(vv,(1-vv))*2
  return(list(bb=bb,pp=pp))
}
