bss=function(ind,log_scale=T){
  bb=c(quantile(var.boot[,ind], probs=c(0.025, 0.975)))
  if(log_scale){
    vv=sum(log(var.boot[,ind])>0)/total
  }else{
    vv=sum(var.boot[,ind]>0)/total
  }
  pp=min(vv,(1-vv))*2
  return(list(bb=bb,pp=pp))
}
