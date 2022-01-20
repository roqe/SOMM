create_vmatrix=function(y.reg, s.reg, q.reg){
  beta.hat=length(y.reg$coefficients)
  alpha.hat=length(s.reg$coefficients)
  delta.hat=length(q.reg$coefficients)
  dimensions=sum(beta.hat+alpha.hat+delta.hat)
  V.matrix=matrix(0, nr=dimensions, nc=dimensions)
  V.matrix[c(1:beta.hat),c(1:beta.hat)]=vcov(y.reg)
  V.matrix[c((1+beta.hat):(beta.hat+alpha.hat)),c((1+beta.hat):(beta.hat+alpha.hat))]=vcov(s.reg)
  V.matrix[c((1+beta.hat+alpha.hat):(beta.hat+alpha.hat+delta.hat)),c((1+beta.hat+alpha.hat):(beta.hat+alpha.hat+delta.hat))]=vcov(q.reg)
  return(V.matrix)
}
