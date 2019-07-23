omega=function(SC,w,confounders=c()){
  b0=SC[1]
  bw=SC[2]
  bq=SC[3]
  bs=SC[4]
  a0=SC[5]
  aw=SC[6]
  aq=SC[7]
  d0=SC[8]
  dw=SC[9]
  sq=SC[10]
  ss=SC[11]
  if(length(SC)>11){
    bX=SC[12]
    aX=SC[13]
    dX=SC[14]
  }else{ bX=aX=dX=0 }
  mu.q=d0+dw*w[3]+sum(dX*confounders)
  if(length(w)==3){
    mu.s=a0+aw*w[2]+aq*mu.q+sum(aX*confounders)
    mu=b0+bw*w[1]+bq*mu.q+bs*mu.s+sum(bX*confounders)
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
    mu.q.star=d0+dw*w[4]+sum(dX*confounders)
    mu.s=a0+aw*w[2]+aq*mu.q.star+sum(aX*confounders)
    mu=b0+bw*w[1]+bq*mu.q+bs*mu.s+sum(bX*confounders)
    sd=sqrt((aq^2*bs^2+bq^2)*(sq^2)+(ss^2)*(bs^2)+1)
    y=mu/sd
    dB0=exp(-y^2/2)/(sqrt(2*pi)*sd) #derivative beta 0
    dbw=dB0*w[1] #derivative beta s
    dbq=dB0*(mu.q-y*bq*sq^2/sd) #derivative beta m
    dbs=dB0*(mu.s-y*(bs*ss^2+bs*aq^2*sq^2)/sd) #derivative beta g
    da0=dB0*bs#derivative for alyha 0
    daw=da0*w[2]#derivative for alyha s
    daq=da0*(mu.q.star-y*sq^2*aq*bs^2/sd) #derivative for alyha m
    dd0=dB0*(bq+bs*aq)#derivative for delta 0
    ddw=dB0*(bq*w[3]+bs*aq*w[4])#derivative for delta s
  }
  dBX=dB0*confounders
  daX=da0*confounders
  ddX=dd0*confounders
  derivatives=c(dB0, dbw, dbq, dbs, dBX, da0, daw, daq, daX, dd0, ddw, ddX)
  return(c(pnorm(y), derivatives))
}
