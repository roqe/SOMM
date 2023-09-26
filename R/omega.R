omega=function(SC,w,cnfd,cnfd2,cnfd3,BF){
  b0=SC$bc[["I"]]; bw=SC$bc[["W"]]; bq=SC$bc[["Q"]]; bs=SC$bc[["S"]]
  a0=SC$ac[["I"]]; aw=SC$ac[["W"]]; aq=SC$ac[["Q"]]
  d0=SC$dc[["I"]]; dw=SC$dc[["W"]]
  sq=SC$sq
  ss=SC$ss
  bl=length(SC$bc)
  if(bl>4){
    bX=SC$bc[5:bl]
    aX=SC$ac[4:(bl-1)]
    dX=SC$dc[3:(bl-2)]
  }else{ bX=aX=dX=0 }
  if(length(cnfd2)==0){ cnfd2=cnfd }
  if(length(cnfd3)==0){ cnfd3=cnfd }
  mu.q=d0+dw*w[3]+sum(dX*cnfd3)
  if(length(w)==3){
    mu.s=a0+aw*w[2]+aq*mu.q+sum(aX*cnfd2)
    mu=b0+bw*w[1]+bq*mu.q+bs*mu.s+sum(bX*cnfd)
    sd=sqrt(((aq*bs+bq)^2)*(sq^2)+(ss^2)*(bs^2)+1)
    y=mu/sd
    dB0=exp(-y^2/2)/(sqrt(2*pi)*sd)
    dbw=dB0*w[1]
    dbq=dB0*(mu.q-y*(bq+bs*aq)*sq^2/sd)
    dbs=dB0*(mu.s-y*(bs*ss^2+(bs*aq+bq)*aq*sq^2)/sd)
    da0=dB0*bs
    daw=da0*w[2]
    daq=da0*(mu.q-y*(bq+bs*aq)*sq^2/sd)
    dd0=dB0*(bq+bs*aq)
    ddw=dd0*w[3]
  }else if(length(w)==4){
    mu.q.star=d0+dw*w[4]+sum(dX*cnfd3)
    mu.s=a0+aw*w[2]+aq*mu.q.star+sum(aX*cnfd2)
    mu=b0+bw*w[1]+bq*mu.q+bs*mu.s+sum(bX*cnfd)
    sd=sqrt((aq^2*bs^2+bq^2)*(sq^2)+(ss^2)*(bs^2)+1)
    y=mu/sd
    dB0=exp(-y^2/2)/(sqrt(2*pi)*sd) #derivative beta 0
    dbw=dB0*w[1] #derivative beta s
    dbq=dB0*(mu.q-y*bq*sq^2/sd) #derivative beta m
    dbs=dB0*(mu.s-y*(bs*ss^2+bs*aq^2*sq^2)/sd) #derivative beta g
    da0=dB0*bs#derivative for alpha 0
    daw=da0*w[2]#derivative for alpha s
    daq=da0*(mu.q.star-y*sq^2*aq*bs^2/sd) #derivative for alpha m
    dd0=dB0*(bq+bs*aq)#derivative for delta 0
    ddw=dB0*(bq*w[3]+bs*aq*w[4])#derivative for delta s
  }
  dBX=dB0*cnfd
  daX=da0*cnfd2
  ddX=dd0*cnfd3
  derivatives=c(dB0, dbw, dbq, dbs, dBX, da0, daw, daq, daX, dd0, ddw, ddX)
  return(c(ifelse(BF,pnorm(y),y), derivatives))
}
