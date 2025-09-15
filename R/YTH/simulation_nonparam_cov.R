library(np)
set.seed(12345)

#setwd("C:/Users/ythuang/Desktop/SI_mediation/simulation_cov/")

## Only linear model works, i.e., type=1

n <- 500
type <- 1
n.sim <- 1000
a10<-0.0; a11<-0.5
a20<-0.0; a21<-0.5
a30<-0.0; a31<-1.0
b1<-1.0; b2<-1.0; b3<-1.0
x<-rnorm(n)

trueY<-function(za=1, zb=0, type){
	if (type==1){
		eta<-za+(b1*a10+b2*a20+b3*a30)+(a11*b1+a21*b2+a31*b3)*zb
	}
	if (type==2){
		eta<-za+b1*(a10+a11*zb)+b2*(a20+a21*zb)+b3*(a30+a31*zb)+0.5*(b1^2+b2^2+b3^2)
		eta<-exp(eta)
	} 
	if (type==3){
		eta<-b1^2+b2^2+b3^2 + (za+b1*(a10+a11*zb)+b2*(a20+a21*zb)+b3*(a30+a31*zb))^2
	}
	return(eta)
}

True<-NULL
for (za in seq(-0.0, 1.0, 0.1)){
	true<-NULL
	for (zb in seq(-0.0, 1.0, 0.1)){
		true<-c(true, trueY(za, zb, type=type))
	}
	True<-cbind(True, true)
}

param<-NULL
EST<-list()
for (k in 1:n.sim){

	z <- x+runif(n, -0.5, 1.5)
	m1 <- a10+x+a11*z+rnorm(n)
	m2 <- a20+x+a21*z+rnorm(n)
	m3 <- a30+x+a31*z+rnorm(n)
	
	if (type==1) y <- x + 1*z + b1*m1 + b2*m2 + b3*m3 + rnorm(n)
	if (type==2) y <- x + exp(1*z + b1*m1 + b2*m2 + b3*m3) + rnorm(n)
	if (type==3) y <- x + (1*z + b1*m1 + b2*m2 + b3*m3)^2 + rnorm(n)

	bwy <- npindexbw(formula=y~z+m1+m2+m3+x)
	
	npmed<-function(za=0, zb=0){
	
		dat1 <- data.frame(z=za, m1=m1, m2=m2, m3=m3, x=0)
		ghat<-predict(npindex(bwy, bw=bwy), newdata=dat1)
		bwm <- npindexbw(formula=ghat~z+x)
		dat2 <- data.frame(z=rep(zb, bwm$nobs), x=0)
		EY.zazb<-predict(npindex(bwm, bw=bwm), newdata=dat2)[1]
	
		return(EY.zazb)
	}
	
	Estimate<-NULL
	for (za in seq(-0.0, 1.0, 0.1)){
		est<-NULL
		for (zb in seq(-0.0, 1.0, 0.1)){
			est<-c(est, npmed(za, zb))
		}
		Estimate<-cbind(Estimate, est)
	}
	EST[[k]]<-Estimate
	
	param<-rbind(param, c(bwy$beta, npmed(1, 1), npmed(1, 0), trueY(zb=1, type=type), trueY(zb=0, type=type)))
	print(apply(param, 2, mean))
	
}

EST.sum<-0
for (l in 1:length(EST)) EST.sum<-EST.sum+EST[[l]]
EST.mean<-EST.sum/length(EST)

col.pal<-colorRampPalette(c("blue", "red"))
colors<-col.pal(100)
zz<-EST.mean
z.facet.center <- (zz[-1, -1] + zz[-1, -ncol(zz)] + zz[-nrow(zz), -1] + zz[-nrow(zz), -ncol(zz)])/4
z.facet.range<-cut(z.facet.center, 100)

par(mfrow=c(2,1))
persp(x=seq(-0.0, 1.0, length.out=nrow(zz)), 
	y=seq(-0.0, 1.0, length.out=ncol(zz)),
	z=zz, zlim=range(True), xlab="zb", ylab="za", zlab="EYzazb",
	col=colors[z.facet.range], box=TRUE, border="grey", 
	theta=-30, phi=10, r=4
)

persp(x=seq(-0.0, 1.0, length.out=nrow(True)), 
	y=seq(-0.0, 1.0, length.out=ncol(True)),
	z=True, zlim=range(True), xlab="zb", ylab="za", zlab="EYzazb",
	col=colors[z.facet.range], box=TRUE, border="grey", 
	theta=-30, phi=10, r=4
)

save(list = ls(all.names = TRUE), file = paste("sim_cov_type", type, "_n", n, "_m0.RData", sep=""), envir = .GlobalEnv)


