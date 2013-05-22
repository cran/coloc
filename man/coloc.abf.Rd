\name{coloc.abf}
\alias{coloc.abf}
\title{Fully Bayesian colocalisation analysis using p values}
\usage{
  coloc.abf(pvalues.dataset1, pvalues.dataset2, MAF,
    N.dataset1, N.dataset2, type.dataset1 = "quant",
    type.dataset2 = "quant", p1 = 1e-04, p2 = 1e-04,
    p12 = 1e-05, s.dataset1 = 0.5, s.dataset2 = 0.5)
}
\arguments{
  \item{pvalues.dataset1}{single variant P-values in
  dataset 1}

  \item{pvalues.dataset2}{single variant P-values in
  dataset 2}

  \item{MAF}{minor allele frequency of the variants}

  \item{N.dataset1}{number of individuals in dataset 1}

  \item{N.dataset2}{number of individuals in dataset 2}

  \item{type.dataset1}{the type of data in dataset 1 -
  either "quant" or "cc" to denote quantitative or
  case-control}

  \item{type.dataset2}{the type of data in dataset 2}

  \item{p1}{prior probability a SNP is associated with
  trait 1}

  \item{p2}{prior probability a SNP is associated with
  trait 2}

  \item{p12}{prior probability a SNP is associated with
  both traits}

  \item{s.dataset1}{the proportion of samples in dataset 1
  that are cases (only relevant for case control samples)}

  \item{s.dataset2}{the proportion of samples in dataset 2
  that are cases}
}
\value{
  a list of two \code{data.frame}s: \itemize{ \item results
  is a vector giving the number of SNPs analysed, and the
  posterior probabilities of H0 (no causal variant), H1
  (causal variant for trait 1 only), H2 (causal variant for
  trait 2 only), H3 (two distinct causal variants) and H4
  (one common causal variant) \item merged.df is an
  annotated version of the input \code{data.frame} }
}
\description{
  Bayesian colocalisation analysis using summary p values
}
\details{
  This function makes a data frame obtained by merging
  p-values for two traits over the same set of SNPs and
  calculates posterior probabilities of different causal
  variant configurations under the assumption of a single
  causal variant for each trait.

  It uses an approximation to calculate Bayes factors for
  association at each SNP that depends on the SNP's MAF and
  ignores any uncertainty in imputation, and is therefore
  better suited to well imputed SNPs.
  \code{\link{coloc.abf.imputed}} is preferred when
  regression coefficients and their standard errors are
  available.
}
\author{
  Claudia Giambartolomei, Chris Wallace
}

