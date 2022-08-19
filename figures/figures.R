
#' @title Draws histogram.
#'
#' @description
#'  Utility function to draw histograms with more options than \code{hist} allows.
#'  
#' @param height Height of histogram bars.
#' @param breaks Locations of boundaries of histogram bins (must be 1 longer than \code{height}).
#' @param lineonly If TRUE, uses \code{\link{lines}} to draw lines on current plot; else uses 
#' \code{\link{plot}} to draw lines on new plot.
#' @param outline If TRUE, draws only the outline (profile) of the histogram; else draws each 
#' complete bar.
#' @param fill If TRUE, uses polygon() to fill barsl in this case valid arguments to polygon() 
#' are passed via argument(s) "...". If fill==FALSE, valid arguments to plot() or lines() are 
#' passed via argument(s) "..."
#' @param ylim Range of y-axis.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param ... See aargument \code{fill}.
#' 
#' @export histline
histline = function(height,breaks,lineonly=FALSE,outline=FALSE,fill=FALSE,ylim=range(height),
                    xlab="x",ylab="y",...)
{
  n=length(height)
  if(length(breaks)!=(n+1)) stop("breaks must be 1 longer than height")
  if(outline) {
    y=c(0,rep(height,times=rep(2,n)),0)
    x=rep(breaks,times=rep(2,(n+1)))
  }   else {
    y=rep(0,4*n)
    x=rep(0,4*n+2)
    for(i in 1:n) {
      y[((i-1)*4+1):(i*4)]=c(0,rep(height[i],2),0)
      x[((i-1)*4+1):(i*4)]=c(rep(breaks[i],2),rep(breaks[i+1],2))
    }
    x=x[1:(4*n)]
  }
  if(lineonly) {
    if(!fill) lines(x,y,...)
    else polygon(x,y,...)
  } else {
    if(!fill) plot(x,y,type="l",ylim=ylim,xlab=xlab,ylab=ylab,...)
    else {
      plot(x,y,type="n",ylim=ylim,xlab=xlab,ylab=ylab)
      polygon(x,y,...)
    }
  }
}


# geometric mean
geomean = function(x) prod(x)^(1/length(x))

#================================= End of Functions ===================================


# Integral figure
# ===============
cube = function(x) x^3 + 2
xmin = -1
xmax = 1.5
x = seq(xmin, xmax, length=500)
plot(x,cube(x)+3,type="l",ylab=expression(f(x)))

pdf(file="./revision/figures/integralfigure.pdf",h=4,w=14)
par(mfrow=c(1,3))
ylim = c(0,cube(xmax))

nbar = 20
dx = (xmax-xmin)/nbar
xs = seq(xmin,xmax,dx)
xmids = xs[-1] - dx/2
barht = cube(xmids)
histline(barht,xs,lineonly=FALSE,outline=FALSE,fill=TRUE,col="gray",ylab=expression(f(x)),ylim=ylim)
lines(x,cube(x))
title("(a)")
approxint1 = sum(barht)*dx

nbar = 50
dx = (xmax-xmin)/nbar
xs = seq(xmin,xmax,dx)
xmids = xs[-1] - dx/2
barht = cube(xmids)
histline(barht,xs,lineonly=FALSE,outline=TRUE,fill=TRUE,col="gray",ylab=expression(f(x)),ylim=ylim)
lines(x,cube(x))
title("(b)")
approxint2 = sum(barht)*dx

nbar = 1000
dx = (xmax-xmin)/nbar
xs = seq(xmin,xmax,dx)
xmids = xs[-1] - dx/2
barht = cube(xmids)
histline(barht,xs,lineonly=FALSE,outline=TRUE,fill=TRUE,col="gray",ylab=expression(f(x)),ylim=ylim)
lines(x,cube(x))
title("(c)")
approxint3 = sum(barht)*dx

dev.off()

approxint = c(approxint1, approxint2, approxint3)
true = (xmax^4/4+2*1.5) - (xmin^4/4+2*xmin)
100*(approxint-true)/true


# Sampling distributions
# ======================
set.seed(1)
y = rnorm(10000000,10,5)
hist(y)
hist(exp(y/10),nclass=30)
smth = density(y,bw=0.75)
indices = rep(1:2000000,times=rep(5,2000000))
ymeans = aggregate(y,by=list(indices=indices),FUN="mean")$x
sdymeans = aggregate(y,by=list(indices=indices),FUN="sd")$x
expymeans = aggregate(exp(y/7),by=list(indices=indices),FUN="mean")$x
geoymeans = aggregate(exp(y/7),by=list(indices=indices),FUN="geomean")$x
meansmth = density(ymeans,bw=0.75)
expmeansmth = density(expymeans,bw=0.1)
geomeansmth = density(geoymeans,bw=0.1)
sdsmth = density(sdymeans,bw=0.1)

pdf(file="./revision/figures/samplingdbns.pdf",h=4,w=12)
par(mfrow=c(1,3))
plot(meansmth$x,meansmth$y,type="l",xlab=expression(hat(mu)),ylab="",main="(a)")
plot(sdsmth$x,sdsmth$y,type="l",xlab=expression(hat(sigma)),ylab="",main="(b)")
#plot(expmeansmth$x,expmeansmth$y,type="l",xlab=expression(bar(exp(y/7))),ylab="",main="(c)")
plot(geomeansmth$x,geomeansmth$y,type="l",xlab=expression((prod(y[i]))^(1/n)),ylab="",main="(c)")
dev.off()

# Prior and posterior figures
# ===========================
N = 100
y = 10
theta = seq(0,1,length=300)
dtheta = unique(diff(theta))[1]
prior = dbeta(theta,3,6)*dtheta
lik = dbinom(y,N,theta)
mle = y/N
lik.mle = dbinom(y,N,mle)
posterior = lik*prior/sum(lik*prior)

pdf(file="./revision/figures/thetalikelihoodmax.pdf",h=4,w=6)
plot(theta,lik,type="l",xlab=expression(theta),ylab=expression(L(theta)),col="blue")
segments(mle,0,mle,lik.mle,col="red",lty=2)
dev.off()


pdf(file="./revision/figures/thetaloglikelihoodmax.pdf",h=4,w=6)
plot(theta,log(lik),type="l",xlab=expression(theta),ylab=expression(l(theta)),col="blue")
segments(mle,-500,mle,log(mle,lik),col="red",lty=2)
dev.off()


pdf(file="./revision/figures/thetalikelihood.pdf",h=4,w=6)
plot(theta,lik,type="l",xlab=expression(theta),ylab=expression(L(theta)),col="blue")
dev.off()


pdf(file="./revision/figures/thetaprior.pdf",h=4,w=6)
plot(theta,prior,type="l",xlab=expression(theta),ylab=expression(P(theta)),col="darkgrey")
dev.off()

ylim=range(c(0,prior,lik/sum(lik),posterior))
pdf(file="./revision/figures/thetaposterior.pdf",h=4,w=6)
plot(theta,prior,type="l",xlab=expression(theta),ylab="",ylim=ylim,col="darkgrey",lty=2)
lines(theta,posterior,col="black")
lines(theta,lik/sum(lik),col="blue",lty=2)
legend("topright",legend=c("prior","likelihood","posterior"),lty=c(2,2,1),col=c("darkgray","blue","black"),bty="n")
dev.off()




