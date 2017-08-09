\name{test.CausalANOVA}
\alias{test.CausalANOVA}
\title{Estimating the AMEs and AMIEs after Regularization with the CausalANOVA.}
\description{
  \code{test.CausalANOVA} estimates the AMEs and AMIEs with confidence intervals after regularization
  with \code{CausalANOVA}function.
}
\usage{
test.CausalANOVA(fit, newdata,collapse.level=TRUE, diff=FALSE, pair.id=NULL,cluster=NULL)
}
\arguments{
  \item{fit}{The output from \code{CausalANOAV} function.}
  \item{newdata}{A data frame to use for re-estimating the AMEs and AMIEs with confidence intervals.}
  \item{collapse.level}{A logical indicating whether to collapse insignificant levels within factors as suggested by the \code{CausalANOVA} output users provide.}
  \item{diff}{A logical indicating whether the outcome is the choice between a pair. If \code{diff=TRUE}, \code{pair.id} should specify a pair of comparison. Default is \code{FALSE}.}
  \item{pair.id}{(optional).Unique identifiers for each pair of comparison. This option is used when \code{diff=TRUE}.}
  \item{cluster}{Unique identifies with which cluster standard errors are computed.}
}

\details{See Details in \code{CausalANOVA}.}
\value{
  \item{fit}{The output of class \code{CausalANOVA}.}
}
\references{
Egami, Naoki and Kosuke Imai. 2016+. ``Causal Interaction in Factorial Experiments:
Application to Conjoint Analysis.'' Working paper. \url{http://imai.princeton.edu/research/files/int.pdf}

Lim, M. and Hastie, T. 2015. Learning interactions via hierarchical group-lasso regularization. Journal of Computational and Graphical Statistics 24, 3, 627--654.

Post, J. B. and Bondell, H. D. 2013. ``Factor selection and structural identification
in the interaction anova model.'' Biometrics 69, 1, 70--79.
}

\author{Naoki Egami and Kosuke Imai.}

\seealso{\link{CausalANOVA}.}