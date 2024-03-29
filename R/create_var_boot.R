strbootQ=function(dt,stra,byGroup=NA){
  if(!is.na(byGroup)){
    dt.raw=dt
    dt=aggregate(as.formula(paste0(".~",byGroup)), data=dt, mean)
  }

  q.range=as.data.frame(sapply(c("Y","Q","W"), function(x) {
    if (length(unique(dt[,x]))==stra) { sort(unique(dt[,x]))
    }else{ unname(quantile(dt[,x], probs=seq(0, 1, 1/(stra-1)))) }
  }))

  dt.BBB=c()
  for (i in 1:stra){
    if (i==1) {
      dt.Y=dt[(dt$Y<=q.range[i,1] & dt$Y>=range(dt.raw$Y)[1]),]
    }else{
      dt.Y=dt[(dt$Y<=q.range[i,1] & dt$Y>q.range[(i-1),1]),]
    }
    dt.BB=c()
    for (j in 1:stra){
      if (j==1) {
        dt.Q=dt[(dt.Y$Q<=q.range[j,2] & dt.Y$Q>=range(dt.raw$Q)[1]),]
      }else{
        dt.Q=dt[(dt.Y$Q<=q.range[j,2] & dt.Y$Q>q.range[(j-1),2]),]
      }
      dt.B=c()
      for (k in 1:stra){
        if (k==1) {
          w_range=which(dt.Q$W<=q.range[k,3] & dt.Q$W>=range(dt.raw$W)[1])
        }else{
          w_range=which(dt.Q$W<=q.range[k,3] & dt.Q$W>q.range[(k-1),3])
        }
        dt.b=as.data.frame(dt.Q[w_range,][sample(1:nrow(dt.Q[w_range,]), replace=T),])
        dt.B=rbind(dt.B, dt.b)
      }
      dt.BB=rbind(dt.BB, dt.B)
    }
    dt.BBB=rbind(dt.BBB, dt.BB)
  }

  if (!is.null(byGroup)) {
    dt.BBB=dt.BB[!is.na(dt.BBB[[byGroup]]),]
    a=table(dt.BBB[[byGroup]])
    dt.boot=c()
    for (AA in min(a):max(a)){
      b=dt.raw[which(dt.raw[[byGroup]] %in% names(which(a>=AA))),]
      dt.boot=rbind(dt.boot,b)
    }
    dt.boot=dt.boot[!is.na(dt.boot$Y),]
  } else {
    dt.boot=dt.BBB[!is.na(dt.BBB$Y),]
  }

  return(dt.boot)
}

# strbootQ=function(dt,stra,byGroup=NA){
#   if(!is.na(byGroup)){
#     dt.raw=dt
#     dt=aggregate(as.formula(paste0(".~",byGroup)), data=dt, mean)
#   }
#   uQ=diff(range(dt$Q))/(stra-1)
#   q_range=c(seq(min(dt$Q),max(dt$Q),by=uQ),max(dt$Q)+0.01)
#   dt.BB=data.table::rbindlist(lapply(1:stra,function(i){
#     dt.Q=dt[dt$Q<q_range[i+1] & dt$Q>=q_range[i],]
#     uW=diff(range(dt.Q$W))/(stra-1)
#     w_range=c(seq(min(dt.Q$W),max(dt.Q$W),by=uW),max(dt.Q$W)+0.01)
#     dt.B=data.table::rbindlist(lapply(1:stra,function(j){
#       ww=(dt.Q$W<w_range[j+1] & dt.Q$W>=w_range[j])
#       return(dt.Q[ww,][sample(1:sum(ww), replace=T),])
#     }))
#   }))
#   if(!is.na(byGroup)){
#     dt.BB=dt.BB[!is.na(dt.BB[[byGroup]]),]
#     a=table(dt.BB[[byGroup]])
#     dt.boot=data.table::rbindlist(lapply(min(a):max(a),function(AA){
#       return(dt.raw[which(dt.raw$id %in% names(which(a>=AA))),])
#     }))
#   }else{
#     dt.boot=dt.BB[!is.na(dt.BB[,ncol(dt.BB)]),]
#   }
#   return(dt.boot)
# }

create_var_boot=function(dt,cnfd=c(),dt2=NULL,cnfd2=c(),dt3=NULL,cnfd3=c(),nb,intv=4,reNAME=NULL,grpID=NULL,x0,x1,BF,stra=1,seed){
  set.seed(seed)
  if(stra>1){
    byGroup=ifelse(is.null(reNAME),ifelse(is.null(grpID),NA,"grp"),"id")
    dt1.boot=strbootQ(dt,stra,byGroup)
    if(!is.null(dt2)){ dt2.boot=strbootQ(dt2,stra,byGroup) }else{ dt2.boot=NULL }
    if(!is.null(dt3)){ dt3.boot=strbootQ(dt3,stra,byGroup) }else{ dt3.boot=NULL }
  }else{
    dt1.boot=as.data.frame(dt[sample(1:nrow(dt), replace=TRUE),])
    if(!is.null(dt2)){ dt2.boot=as.data.frame(dt2[sample(1:nrow(dt2), replace=TRUE),]) }else{ dt2.boot=NULL }
    if(!is.null(dt3)){ dt3.boot=as.data.frame(dt3[sample(1:nrow(dt3), replace=TRUE),]) }else{ dt3.boot=NULL }
  }
  GT=try(get_theta(dt1.boot,dt2.boot,dt3.boot,reNAME,grpID),silent = T)
  if(class(GT)=="try-error"){ return(NULL) }
  if(intv==3){
    PP=PSE_three(GT,x0,x1,cnfd,cnfd2,cnfd3,V.matrix=NULL,BF)
    if(BF){
      if(sum(PP$omega<0)>0){ return(NULL) }else{
        return(data.table(t(c(PP$RD$V1[1:3],PP$RR$V1[1:3],PP$OR$V1[1:3],NA,NA,NA,
                              PP$RD$V1[4],PP$RR$V1[4],PP$OR$V1[4],PP$RD$PM[1:3],PP$RR$PM[1:3],PP$OR$PM[1:3],NA,NA,NA))))
      }
    }else{
      return(data.table(t(c(PP$RD$V1[1:3],rep(NA,9),PP$RD$V1[4],NA,NA,PP$RD$PM[1:3],rep(NA,9)))))
    }
  }else if(intv==4){
    PP=PSE_four(GT,x0,x1,cnfd,cnfd2,cnfd3,V.matrix=NULL,BF)
    if(BF){
      if(sum(PP$omega<0)>0){ return(NULL) }else{
        return(data.table(t(c(PP$RD$V1[1:3],PP$RR$V1[1:3],PP$OR$V1[1:3],PP$RD$V1[4],PP$RR$V1[4],PP$OR$V1[4],
                              PP$RD$V1[5],PP$RR$V1[5],PP$OR$V1[5],PP$RD$PM[1:3],PP$RR$PM[1:3],PP$OR$PM[1:3],PP$RD$PM[4],PP$RR$PM[4],PP$OR$PM[4]))))
      }
    }else{
      return(data.table(t(c(PP$RD$V1[1:3],rep(NA,6),PP$RD$V1[4],NA,NA,PP$RD$V1[5],NA,NA,PP$RD$PM[1:3],rep(NA,6),PP$RD$PM[4],NA,NA))))
    }
  }
 }
