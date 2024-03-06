#' @import data.table
bss=function(ind,var.boot,log_scale,pos){
  dd=var.boot[[ind]]
  if(all(is.na(dd))){ return(list(pp=NA,btsd=NA)) }else{ dd=dd[!is.na(dd)] }
  if(pos){ dd=dd[dd>=0] }
  cc=length(dd)
  if(log_scale){ vv=sum(dd>1)/cc }else{ vv=sum(dd>0)/cc }
  pp=min(vv,(1-vv))*2
  exact=c(quantile(dd,c(.05,.95)),pp)

  if(log_scale){ dd=log(dd) }
  cc=1e+6; dd=rnorm(cc, mean=mean(dd,na.rm = T), sd=sd(dd,na.rm = T))
  vv=sum(dd>0)/cc;
  pp=min(vv,(1-vv))*2
  if(log_scale){ appxn=c(exp(quantile(dd,c(.05,.95))),pp)  }else{ appxn=c(quantile(dd,c(.05,.95)),pp) }

  return(c(exact,appxn))
}
