### coloc class: simple class to hold results of coloc.test()

validColoc <- function(object) {
  if(length(object@result) != 4 ||
     !all(c("eta.hat","statistic","df","ppp") %in% names(object@result))) {
    return("result slot in coloc objects should be a named numeric vector with 4 elements: eta.hat, statistic, df, ppp")
  }
}
setClass("coloc",
         representation(result="numeric"),
         validity=validColoc)

show.coloc <- function(object) {
  print(object@result)
}
setMethod("show","coloc",show.coloc)

eta <- function(object) {
  object@result["eta.hat"]
}
statistic <- function(object) {
  object@result["statistic"]
}
df <- function(object) {
  object@result["df"]
}
ppp.value <- function(object) {
  object@result["ppp"]
}
p.value <- function(object) {
  pchisq(object@result["statistic"],df=object@result["df"],lower=FALSE)
}
setGeneric("eta",function(object) standardGeneric("eta"))
setMethod("eta","coloc",function(object) object@result["eta.hat"])
setGeneric("theta",function(object) standardGeneric("theta"))
setMethod("theta","coloc",function(object) {
          theta <- atan(object@result["eta.hat"])
          names(theta) <- "theta.hat"
          theta
          })
setGeneric("chisquare",function(object) standardGeneric("chisquare"))
setMethod("chisquare","coloc",function(object) object@result["statistic"])
setGeneric("df",function(object) standardGeneric("df"))
setMethod("df","coloc",function(object) object@result["df"])
setGeneric("ppp.value",function(object) standardGeneric("ppp.value"))
setMethod("ppp.value","coloc",function(object) object@result["ppp"])
setGeneric("p.value",function(object) standardGeneric("p.value"))
setMethod("p.value","coloc",function(object)
          pchisq(object@result["statistic"],df=object@result["df"],lower=FALSE))



coloc.test <- function(X,Y,k=1,plot.coeff=TRUE,plots.extra=NULL,
                       vars.drop=NULL,...) {
  ## X and Y are glm objects, fitted to the same snps, with different outcome variables
  ## return values are
  ## return(c(eta.hat=eta.hat,chisq=X2,ppp=ppp$value))
  ## where
  ## eta.hat is the estimated slope
  ## chisq is the test statistic (degrees of freedom <= number of snps)
  ## ppp is the posterior predictive p value

  vars.drop <- c(vars.drop,"(Intercept)")
  snps <- setdiff(intersect(names(coefficients(X)),names(coefficients(Y))),
                  vars.drop)
  snps.dropX <- setdiff(names(coefficients(X)),c(snps,vars.drop))
  snps.dropY <- setdiff(names(coefficients(Y)),c(snps,vars.drop))
  if(length(snps.dropX))
    cat("Warning:",length(snps.dropX),"variables dropped from regression X:\n\t",snps.dropX,"\n")
  if(length(snps.dropY))
    cat("Warning:",length(snps.dropY),"variables dropped from regression Y:\n\t",snps.dropY,"\n")
  
  if(length(snps)<=1) { # 1 common coef => equal already
    cat("Only 1 factor,",snps," in common.  Skipping\n")
    return(c(NA,NA,NA))
  }
  b1 <- coefficients(X)[snps]
  b2 <- coefficients(Y)[snps]
  V1 <- vcov(X)[snps,snps]
  V2 <- vcov(Y)[snps,snps]

  ## posterior for theta
  S1 <- solve(V1)
  S2 <- solve(V2)
  p <- length(b1)
  const <- ( sqrt(2*pi)^p * det(V1) * det(V2) )^(-1)
  M <- function(theta) { solve(S1 + S2 * tan(theta)^2) }
  mu <- function(theta) { t( (t(b1) %*% S1 + tan(theta) * t(b2) %*% S2) %*% M(theta) ) }

  ## beta prior
  pr <- function(theta) { tt <- tan(theta);
                          k*(1+tt^2) / (2*pi*(1 + k^2 * tt^2)) }
  L <- function(theta) {
    const * pr(theta) * det(M(theta)) * exp( -0.5 * (t(b1) %*% S1 %*% b1 + t(b2) %*% S2 %*% b2 -
                                        t(mu(theta)) %*% solve(M(theta)) %*% mu(theta)) )
  }
  LV <- Vectorize(L,"theta")
  LV.int <- integrate(LV,lower=0,upper=pi)
  post <- function(theta) { LV(theta) / LV.int$value }

  ## chisq values
  d <- function(theta) { b1 - tan(theta)^(-1) * b2 }
  S <- function(theta) { V1 + tan(theta)^(-2) * V2 }
  T <- function(theta) { t(d(theta)) %*% solve(S(theta)) %*% d(theta) }
  TV <- Vectorize(T,"theta")
  pv <- function(theta) { pchisq(T(theta),df=p,lower=FALSE) }
  pval <- Vectorize(pv,"theta")

  ## final p value
  toint <- function(theta) { pval(theta) * post(theta) }
  ppp <- integrate(toint,lower=0,upper=pi)
  x <- seq(0,pi,length=1001)
  t <- TV(x)
  theta.hat <- x[which.min(TV(x))]
  eta.hat <- tan(theta.hat)
  X2 <- TV(theta.hat)

  ## plots
  if(plot.coeff) {
    bplot(b1,b2,diag(V1),diag(V2),eta=eta.hat,
          main="Coefficients",
          sub=paste("ppp =",format.pval(ppp$value)),xlab="beta_1",ylab="beta_2")
  }
  
  x <- seq(-pi,pi,length=101)
  
  if(!is.null(plots.extra)) {
      plot.data <- list(theta=x,
                        eta=tan(x),
                        chisq=TV(x),
                        post.theta=post(x),
                        lhood=LV(x))
    if(!is.list(plots.extra) || length(plots.extra)!=2 ||
       !("x" %in% names(plots.extra)) || !("y" %in% names(plots.extra)) ||
#       length(plots.extra$x)!=length(plots.extra$y) ||
       !all(plots.extra$x %in% names(plot.data)) ||
       !all(plots.extra$y %in% names(plot.data))) {
      warning("plots.extra must be of the format list(x=..., y=...) where x and y are equal length character vectors with elements from theta, eta, lhood, chisq, theta.post.  Skipping plots.extra.")
    } else {
      if(length(plots.extra$y)>length(plots.extra$x)) {
        plots.extra$x <- rep(plots.extra$x,length=length(plots.extra$y))
      }
      if(length(plots.extra$x)>length(plots.extra$y)) {
        plots.extra$y <- rep(plots.extra$y,length=length(plots.extra$x))
      }
        
      for(i in 1:length(plots.extra$x)) {        
      plot(plot.data[[ plots.extra$x[i] ]],
           plot.data[[ plots.extra$y[i] ]],
           type="l",axes=FALSE,
           xlab=plots.extra$x[i],ylab=plots.extra$y[i],
           main=paste(plots.extra$y[i],"vs",plots.extra$x[i]),...)
      if(plots.extra$x[i]=="theta") {
        axis(1,at=seq(0,pi,length=5),labels=c("0","pi/4","pi/2","3pi/4","pi"))
      } else {
        axis(1)
      }
      axis(2); box()
      if(plots.extra$x[i]=="theta") 
        abline(v=theta.hat,col="blue",lty=3)
      if(plots.extra$x[i]=="eta") 
        abline(v=eta.hat,col="blue",lty=3)
      
    }
    }
    }

  ## return
  new("coloc",
      result=c(eta.hat=eta.hat,statistic=X2,df=length(snps)-1,ppp=ppp$value))  
}
