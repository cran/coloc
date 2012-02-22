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
    \item{\code{result}:}{Named numeric vector containing maximum
      likelihood estimate and associated chisquare and posterior
      predictive p values}.
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
    \item{ppp.value}{\code{signature(object = "coloc")}: returns the posterior predicted p value. }
    \item{show}{\code{signature(object = "coloc")}: print a coloc object. }
  }
}
\references{
  Plagnol et al (2009).
  Statistical independence of the colocalized association signals for
  type 1 diabetes and RPS26 gene expression on chromosome
  12q13. Biostatistics
  10:327-34. \url{http://www.ncbi.nlm.nih.gov/pubmed/19039033}

  Wallace et al (in preparation).
%%  ~~put references to the literature/web site here~~
}
\author{
  Chris Wallace.
%%  ~~who you are~~
}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{coloc.test}}
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("coloc")
}
\keyword{classes}
