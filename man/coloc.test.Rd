\name{coloc.test}
\alias{coloc.test}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to do colocalisation tests of two traits
}
\description{
Performs the colocalisation tests described in Plagnol et al (2009) and Wallace et al, in preparation.
}
\usage{
coloc.test(X, Y, k=1,
plot.coeff=TRUE,
plots.extra=NULL,
vars.drop=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{X, Y}{ Each is either an lm or glm object.  The intersection
  of \code{names(coefficients(X))} and \code{names(coefficients(Y))}
  is used to identify SNPs in common which will be tested for
  colocalisation.  Any Intercept term is dropped, but other covariates
  should have distinct names or be listed in \code{vars.drop} to avoid
  them being included in the colocalisation test.}
\item{k}{Theta has a Cauchy(0,k) prior.  The default, k=1, is equivalent
  to a uniform (uninformative) prior.  We have found
  varying k to have little effect on the results.}
\item{plot.coeff}{
  \code{TRUE} if you want to generate a plot showing the coefficients
  from the two regressions together with confidence regions.}
\item{plots.extra}{
  list with 2 named elements, x and y, equal length character vectors
  containing the names of the quantities to be plotted on the x and y
  axes.

  \code{x} is generally a sequence of \code{theta} and \code{eta}, with
  \code{y} selected from 
  \code{theta.post}, the posterior density of theta, 
  \code{chisq}, the chi-square values
  of the test, and \code{lhood}, the likelihood function.}
\item{vars.drop}{Character vector naming additional variables in either
  regression which are not  SNPs and should not be used in the
  colocalisation test.  They should  appear in
  \code{c(names(coefficients(X)),names(coefficients(Y)))}}
\item{...}{Additional arguments passed to plotting functions}
}
\details{
  Please see Wallace et al, in preparation, for more details of the
  model and test.
}
\value{
  a numeric vector with 3 named elements:
  \item{eta.hat}{The estimated slope.}
  \item{statistic}{The test statistic}
  \item{n}{The number of snps used in the test.  If eta were known, this
    would be the degrees of freedom of the test. Because
    eta has been replaced by its ML estimate, Plagnol et al suggest we expect the
    degrees of freedom to be n-1, but this requires the likelihood to be
    well behaved which is not always the case.  We prefer to consider
    the posterior predictive p value.}
  \item{ppp}{The posterior predictive p value}
}
\references{
  Plagnol et al (2009).
  Statistical independence of the colocalized association signals for
  type 1 diabetes and RPS26 gene expression on chromosome
  12q13. Biostatistics
  10:327-34. \url{http://www.ncbi.nlm.nih.gov/pubmed/19039033}

  Wallace et al (in preparation).
%% ~put references to the literature/web site here ~
}
\author{
  Chris Wallace.
}
\note{
  For Plagnol et al's original test only, please see an original package by Vincent
  Plagnol: 
  \code{QTLMatch v0.8} http://www-gene.cimr.cam.ac.uk/vplagnol/software.shtml
}

\examples{
  ## simulate covariate matrix (X) and continuous response vector (Y)
  ## for two populations/triats Y1 and Y2 depend equally on f1 and f2
  ## within each population, although their distributions differ between
  ## populations.  They are compatible with a null hypothesis that they
  ## share a common causal variant
  X1 <- matrix(rbinom(1000,1,0.4),ncol=2)
  Y1 <- rnorm(500,apply(X1,1,sum),2)
  X2 <- matrix(rbinom(1000,1,0.6),ncol=2)
  Y2 <- rnorm(500,2*apply(X2,1,sum),5)
  
  boxplot(list(Y1,Y2),names=c("Y1","Y2"))
  
  ## fit and store linear model objects
  colnames(X1) <- colnames(X2) <- c("f1","f2")
  summary(lm1 <- lm(Y1~f1+f2,data=as.data.frame(X1)))
  summary(lm2 <- lm(Y2~f1+f2,data=as.data.frame(X2)))
  
  ## test whether the traits are compatible with colocalisation
  ### ppp should be large (>0.05, for example), indicating that they are.
  par(mfrow=c(2,2))
  coloc.test(lm1,lm2,plot.coeff=TRUE,
             plots.extra=list(x=c("eta","theta"),
                              y=c("lhood","theta.post")))

}
%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
