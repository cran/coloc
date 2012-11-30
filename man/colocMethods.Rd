\name{Methods for coloc}
\alias{eta}
\alias{theta}
\alias{df}
\alias{ppp.value}
\alias{p.value}
\alias{chisquare}
\alias{df,coloc-method}
\alias{ci}
\alias{ci,coloc-method}
\alias{n,coloc-method}
\alias{eta,coloc-method}
\alias{ppp.value,coloc-method}
\alias{p.value,coloc-method}
\alias{show,coloc-method}
\alias{chisquare,coloc-method}
\alias{theta,coloc-method}

\title{
Functions to extract information from a \code{coloc} object
}
\usage{
\S4method{eta}{coloc}(object)
}
\description{
  Extract information from a \code{coloc} object.
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{ Object returned by \code{coloc.test()} function.}
}
\details{
  \code{eta()} returns eta.hat, the maximum likelihood value of eta.
  
  \code{theta()} returns theta.hat, the maximum likelihood value of eta.
  
  \code{lhood()} returns -2 times the log-likelihood ratio..
  
  \code{chisquare()} returns the value of the chisquare statistic derived by
  Wallace et al (in preparation).
  
  \code{df()} returns the associated degrees of freedom.
  
  \code{p.value()} returns the associated p value.
  
  \code{ppp.value()} returns the posterior predicted p value, or \code{NA} if
  not calculated.

  \code{ci()} returns the credible interval, or \code{NA} if not calculated.  
}

\author{
  Chris Wallace.
}

\seealso{\code{\link{coloc.test}}}
