\name{cv.CausalANOVA}
\alias{cv.CausalANOVA}
\alias{plot.cv.CausalANOVA}
\title{Cross validation for the CausalANOVA.}
\description{
  \code{cv.CausalANOVA} implements cross-validation for \code{CausalANOVA} to select the \code{collapse.cost} parameter.
  \code{CausalANOVA} runs this function internally when defaults when \code{collapse.type=cv.min} or \code{collapse.type=cv.1Std}.
}
\usage{
cv.CausalANOVA(formula, int2.formula=NULL, int3.formula=NULL,
	       data, nway=1,pair.id=NULL, diff=FALSE,
	       cv.collapse.cost=c(0.1,0.3,0.7), nfolds=5,
	       screen=TRUE, screen.type="fixed", screen.num.int=3,
	       family="binomial", cluster=NULL,
	       maxIter=50, eps=1e-5, seed=1234,
	       fac.level=NULL,ord.fac=NULL, verbose=TRUE)
}
\arguments{
  \item{formula}{a formula that specifies outcome and treatment variables.}
  \item{int2.formula}{(optional). A formula that specifies two-way interactions.}
  \item{int3.formula}{(optional). A formula that specifies three-way interactions.}
  \item{data}{an optional data frame, list or environment (or object
              coercible by 'as.data.frame' to a data frame) containing the
              variables in the model.  If not found in 'data', the
              variables are taken from 'environment(formula)', typically
              the environment from which 'CausalANOVA' is called.}
 \item{nway}{With \code{nway=1}, the function estimates the Average Marginal Effects (AMEs) only.
  With \code{nway=2},
  the function estimates the AMEs and the two-way Average Marginal Interaction Effects (AMIEs).
  With \code{nway=3},
  the function estiamtes the AMEs, the two-way and three-way AMIEs. Default is 1.}
  \item{diff}{A logical indicating whether the outcome is the choice between a pair.
  If \code{diff=TRUE}, \code{pair.id} should specify a pair of comparison. Default is \code{FALSE}.}
  \item{pair.id}{(optional).Unique identifiers for each pair of comparison. This option is used when
  \code{diff=TRUE}.}
 \item{cv.collapse.cost}{A vector containing candidates for a cost parameter ranging from 0 to 1.
  1 corresponds to no regularization and the smaller value corresponds to the stronger regularization.
  Default is \code{c(0.1,0.3,0.7)}.}
  \item{nfolds}{number of folds - default is 5. Although nfolds can be as large as the sample
  size (leave-one-out CV), it is not recommended for large datasets.}
  \item{screen}{A logical indicating whether select significant factor interactions with \code{glinternet}.
  When users specify interactions using \code{int2.formula} or \code{int3.formula}, this option is ignored.
  \code{screen} should be used only when users want data-driven selection of factor-interactions.
  With \code{screen.type},
  users can specify how to screen factor interactions. We recommend to use this option when the number of
  factors is large,
  e.g., more than 6. Default is \code{FALSE}.}
  \item{screen.type}{Type for screening factor interactions. (1) \code{"fixed"} select the fixed number (specified by
  \code{screen.num.int}) of factor interactions. (2) \code{"cv.min"} selects factor-interactions with the
  tuning parameter
  giving the minimum cross-validation error. (3) \code{"cv.1Std"} selects factor-interactions with the
  tuning parameter
  giving a cross-validation error that is within 1 standard deviation of the minimum cv error.}
  \item{screen.num.int}{(optional).The number of factor interactions to select. This option is used when and
  \code{screen=TRUE} and \code{screen.type="fixed"}. Default is 3.}
  \item{family}{A family of outcome varialbes. \code{"gaussian"} when continuous outcomes \code{"binomial"} when binary outcomes.
  Default is \code{"binomial"}.}
  \item{cluster}{Unique identifies with which cluster standard errors are computed.}
  \item{maxIter}{The number of maximum iteration for \code{glinternet}.}
  \item{eps}{A tolerance parameter in the internal optimization algorithm.}
  \item{seed}{an argument for \code{set.seed()}.}
    \item{fac.level}{optional. A vector containing the number of levels in each factor. The order of \code{fac.level} should match to the order of columns in the data. For example, when the first and second columns of the design matrix is "Education" and "Race", the first and second element of \code{fac.level} should be the number of levels in "Education" and "Race", respectively.}
  \item{ord.fac}{optional. logical vectors indicating whether each factor has ordered (\code{TRUE}) or unordered (\code{FALSE}) levels. When levels are ordered, the function uses the order given by function \code{levels()}. If levels are ordered, the function places penalties on the differences between adjacent levels. If levels are unordered, the function places penalties on the differences based on every pairwise comparison.}
  \item{verbose}{whether it prints the value of a cost parameter used.}
}

\details{See Details in \code{CausalANOVA}.}
\value{
  \item{cv.error}{The mean cross-validated error - a vector of length \code{length(cv.t)}.}
  \item{cv.min}{A value of \code{t} that gives minimum \code{cv.missclass}.}
  \item{cv.1sd}{The largest value of \code{t} such that error is within 1 standard error of the minimum.}
  \item{cv.each.mat}{A matrix containing cross-validation errors for each fold and cost parameter.}
  \item{cv.cost}{The \code{cv.collapse.cost} used in the function.}
}
\references{
Post, J. B. and Bondell, H. D. 2013. ``Factor selection and structural identification
in the interaction anova model.'' Biometrics 69, 1, 70--79.

Egami, Naoki and Kosuke Imai. 2016+. ``Causal Interaction in Factorial Experiments:
Application to Conjoint Analysis.'' Working paper. \url{http://imai.princeton.edu/research/files/int.pdf}
}

\author{Naoki Egami and Kosuke Imai.}

\seealso{\link{CausalANOVA}.}