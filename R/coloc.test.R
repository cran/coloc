### coloc class: simple class to hold results of coloc.test()
validColoc <- function(object) {
  if(length(object@result) != 4 ||
     !all(c("eta.hat","chisquare","n","ppp") %in% names(object@result))) {
    return("result slot in coloc objects should be a named numeric vector with 4 elements: eta.hat, chisquare, n, ppp, p.boot")
  }
  ## if(length(object@u) != nrow(object@V))
  ##   return("dim of V should equal length of u")
}
setClass("coloc",
         representation(result="numeric", credible.interval="list"),
         validity=validColoc)

show.coloc <- function(object) {
  res <- c(object@result,  p.value=p.value(object), ci(object))
  print(res)
}
setMethod("show","coloc",show.coloc)

## eta <- function(object) {
##   object@result["eta.hat"]
## }
## chisquare <- function(object) {
##   object@result["chisquare"]
## }
## df <- function(object) {
##   object@result["df"]
## }
## df <- function(object) {
##   object@result["n"]
## }
## ppp.value <- function(object) {
##   object@result["ppp"]
## }
## p.value <- function(object) {
##   pchisq(object@result["chisquare"],df=object@result["df"],lower.tail=FALSE)
## }
setGeneric("eta",function(object) standardGeneric("eta"))
setMethod("eta","coloc",function(object) object@result["eta.hat"])
setGeneric("theta",function(object) standardGeneric("theta"))
setMethod("theta","coloc",function(object) {
          theta <- atan(object@result["eta.hat"])
          names(theta) <- "theta.hat"
          theta
          })
setGeneric("chisquare",function(object) standardGeneric("chisquare"))
setMethod("chisquare","coloc",function(object) object@result["chisquare"])
setGeneric("ci",function(object) standardGeneric("ci"))
setMethod("ci","coloc",function(object) object@credible.interval[c("mode","lower","upper","level.observed.value","interior")])
setGeneric("df",function(object) standardGeneric("df"))
setMethod("df","coloc",function(object) object@result["n"]-1)
setGeneric("n",function(object) standardGeneric("n"))
setMethod("n","coloc",function(object) object@result["n"])
setGeneric("ppp.value",function(object) standardGeneric("ppp.value"))
setMethod("ppp.value","coloc",function(object) object@result["ppp"])
setGeneric("p.value",function(object) standardGeneric("p.value"))
setMethod("p.value","coloc",function(object)
          pchisq(object@result["chisquare"],df=object@result["n"]-1,lower.tail=FALSE))


### function to do colocalisation testing

coloc.test <- function(X,Y,k=1,plot.coeff=TRUE,plots.extra=NULL,
                       vars.drop=NULL, bayes.ppp=TRUE, bayes.ci=TRUE,n.ci=1001,level.ci=0.95) {
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
  S1 <- solve(V1)
  S2 <- solve(V2)
  theta.min <- 0
  theta.max <- pi

  ## -2LL = Fieller's chisq
  d <- function(theta,b1,b2) { sin(theta) * b1 - cos(theta) * b2 }
  Vstar <- function(theta) { sin(theta)^2 * V1 + cos(theta)^2 * V2 }
  chisq <- function(theta,b1,b2) { t(d(theta,b1,b2)) %*% solve(Vstar(theta)) %*% d(theta,b1,b2) }
  chisqV <- Vectorize(chisq, "theta")

  findmin <- function(b1,b2) {
    ## there are at most two minima, and never both on the same side of pi/2
    o.left <- optimize(chisq,interval=c(0,pi/2),b1=b1,b2=b2)
    o.right <- optimize(chisq,interval=c(pi/2,pi),b1=b1,b2=b2)
    if(o.left$objective < o.right$objective) {
      return(o.left)
    } else {
      return(o.right)
    }
  }
  fm <- findmin(b1,b2)
  theta.hat <- fm$minimum; eta.hat=tan(theta.hat)
  X2 <- fm$objective[1,1]

  if(bayes.ppp || bayes.ci) {
    ## cauchy prior for theta
    prior <- function(theta) { tt <- tan(theta);
                               k*(1+tt^2) / (2*pi*(1 + k^2 * tt^2)) }
    ## posterior for theta
    p <- length(b1)
    const <- ( sqrt(2*pi)^p * det(V1) * det(V2) )^(-1)
    M <- function(theta) { solve(cos(theta)^2 * S1 + sin(theta)^2 * S2) }
    mu <- function(theta) { t( (cos(theta) * t(b1) %*% S1 +
                                sin(theta) * t(b2) %*% S2) %*% M(theta) ) }
  L <- function(theta) {
    const * prior(theta) * det(M(theta))^(-0.5) *
      exp( -0.5 * (t(b1) %*% S1 %*% b1 + t(b2) %*% S2 %*% b2 -
                   t(mu(theta)) %*% solve(M(theta)) %*% mu(theta)) )
  }
  LV <- Vectorize(L,"theta")
  LV.int <- integrate(LV,lower=0,upper=pi)
  post <- function(theta) { LV(theta) / LV.int$value }
  }

  if(bayes.ppp) {
  ##  posterior predictive p value
  pv <- function(theta) { pchisq(chisq(theta,b1,b2),df=p,lower.tail=FALSE) }
  pval <- Vectorize(pv,"theta")
  toint <- function(theta) { pval(theta) * post(theta) }
  ppp <- integrate(toint,lower=theta.min,upper=theta.max)
} else {
  ppp <- NA
}

  if(bayes.ci) {
  ## credible interval
  
  ## minimum, so that we can reset theta.min, theta.max and centre on the mode
  minimum <- optimize(post, interval=c(theta.min, theta.max))
  t.min <- minimum$minimum
  t.max <- minimum$minimum + pi
  ## mode, within that range
  mode <- optimize(post, interval=c(t.min, t.max), maximum=TRUE)

  pvec <- seq(t.min, t.max, length=n.ci)
  pval <- post(pvec)
  pval <- pval/sum(pval)

  ## search
  centre <- which.max(pval)
  v <- pval[centre-1]
  l <- which.min(abs(pval[1:(centre+1)] - v))
  u <- which.min(abs(pval[(centre+1):n.ci] - v)) + centre
  while(sum(pval[l:u])<level.ci) { # go down in twenties
    cat(".")
    v <- max(pval[l-20], pval[u+20])
    l <- which.min(abs(pval[1:(centre+1)] - v))
    u <- which.min(abs(pval[(centre+1):n.ci] - v)) + centre
  }
  while(sum(pval[l:u])>level.ci) { # go up in ones
    cat(".")
    v <- min(pval[l+1], pval[u-1])
    l <- which.min(abs(pval[1:(centre+1)] - v))
    u <- which.min(abs(pval[(centre+1):n.ci] - v)) + centre
  }
  ci <- c(mode=tan(mode$maximum), lower=tan(pvec[l]), upper=tan(pvec[u]), level=level.ci, level.observed=integrate(post, lower=pvec[l], upper=pvec[u]), interior=pvec[l] %/% pi == pvec[u] %/% pi)
} else {
  ci <- NA
}
   
################################################################################
  ## plots
  if(plot.coeff) {
    coeff.plot(b1,b2,diag(V1),diag(V2),eta=eta.hat,
          main="Coefficients",
#         sub=paste("ppp =",format.pval(ppp$value,digits=2),"p =",format.pval(pchisq(X2,df=length(snps)-1,lower.tail=FALSE),digits=2)),
          xlab=expression(b[1]),ylab=expression(b[2]))
  }

  x <- seq(theta.min,theta.max,length=1001)

  if(!is.null(plots.extra)) {
    plot.data <- list(theta=x,
                      eta=tan(x),
                      chisq=chisqV(x,b1,b2),
                      post.theta=post(x),
                      lhood=chisqV(x,b1,b2))
    extra.plot(plot.data, plots.extra, theta.hat=theta.hat, eta.hat=eta.hat)   
  }

  ## return
  return(new("coloc",
             result=c(eta.hat=eta.hat,chisquare=X2,n=length(snps),
               ppp=ppp$value),
             credible.interval=ci))
             ## u=d(theta.hat, b1, b2),
             ## V=Vstar(theta.hat),
             ## post=post(seq(theta.min,theta.max,0.01))))
}

