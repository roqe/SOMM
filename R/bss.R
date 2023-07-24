#' @import data.table
bss=function(ind,var.boot,log_scale=T,exact=F){
  cc=nrow(var.boot)
  dd=unlist(var.boot[,..ind])
  if(!exact){ cc=1e+6; dd=rnorm(cc, mean=mean(dd), sd=sd(dd)) }

  bb=c(quantile(dd, probs=c(0.025, 0.975), na.rm=T))
  if(log_scale){
    vv=sum(log(dd)>0)/cc
  }else{
    vv=sum(dd>0)/cc
  }
  pp=min(vv,(1-vv))*2
  return(list(bb=bb,pp=pp))
}
