\name{coloc-class}
\Rdversion{1.1}
\docType{class}
\alias{coloc-class}

\title{Class \code{"coloc"}}
\description{
  %%  ~~ A concise (1-5 lines) description of what the class is. ~~
  Class designed to hold objects returned by function
  \code{\link{coloc.test}} which performs a test of the null
  hypothesis that two genetic traits colocalise - that they share a
  common causal variant.
}
\section{Objects from the Class}{
Objects can be created by calls to the function \code{\link{coloc.test}()}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{result}:}{named numeric vector containing maximum
      likelihood estimate and associated chisquare and posterior
      predictive p values}
    \item{\code{credible.interval}:}{list of named elements describing
      the credible interval for \code{eta}}
  }
}
\section{Methods}{
  \describe{
    \item{theta}{\code{signature(object = "coloc")}: returns theta.hat, the maximum likelihood value of theta. }
    \item{eta}{\code{signature(object = "coloc")}: returns eta.hat, the maximum likelihood value of eta. }
    \item{chisquare}{\code{signature(object = "coloc")}: returns the value of the chisquare statistic derived by
  Wallace et al (in preparation). }
    \item{df}{\code{signature(object = "coloc")}: returns the associated degrees of freedom. }
    \item{p.value}{\code{signature(object = "coloc")}: returns the associated p value. }
    \item{ppp.value}{\code{signature(object = "coloc")}: returns the
      posterior predicted p value, or \code{NA} if not calculated. }
    \item{ci}{\code{signature(object = "coloc")}: returns the Bayesian
      credible interval for eta, or \code{NA} if not calculated. }
    \item{show}{\code{signature(object = "coloc")}: print a coloc object. }
  }
}
\references{
  Wallace et al (2012).
  Statistical colocalisation of monocyte gene expression and genetic risk variants for type 1 diabetes.
  Hum Mol Genet 21:2815-2824.
  \url{http://europepmc.org/abstract/MED/22403184}

  Plagnol et al (2009).
  Statistical independence of the colocalized association signals for
  type 1 diabetes and RPS26 gene expression on chromosome
  12q13. Biostatistics
  10:327-34. \url{http://www.ncbi.nlm.nih.gov/pubmed/19039033}
}
\author{
  Chris Wallace.
}
\seealso{
  \code{\link{coloc.test}}
}
\examples{
showClass("coloc")
}
\keyword{classes}
