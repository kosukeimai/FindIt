\name{CausalANOVA}
\alias{CausalANOVA}
\alias{summary.CausalANOVA}
\alias{plot.CausalANOVA}
\title{Estimating the AMEs and AMIEs with the CausalANOVA.}
\description{
  \code{CausalANOVA} estimates coefficients of the specified ANOVA with regularization. By taking differences in coefficients, the function recovers the AMEs and AMIEs. 
}
\usage{	
CausalANOVA(formula,data,cost,pair.id=NULL,nway=2,diff=TRUE,
	    select.prob=FALSE,boot=100,block.id=NULL,seed=1234,
	    eps=1e-5,fac.level=NULL,ord.fac=NULL,verbose=TRUE)
}
\arguments{
  \item{formula}{a formula that specifies outcome and treatment variables.}	
  \item{data}{an optional data frame, list or environment (or object
              coercible by 'as.data.frame' to a data frame) containing the
              variables in the model.  If not found in 'data', the
              variables are taken from 'environment(formula)', typically
              the environment from which 'CausalANOVA' is called.}
  \item{cost}{a cost parameter ranging from 0 to 1. 1 corresponds to no regularization. It is recommended to choose the value using cross validation with \code{cv.GashANOVA}.}
  \item{pair.id}{unique identifiers for each pair of comparison. This option is used when \code{dif=TRUE}.}
  \item{nway}{"2" when the two way causal interactions are of interest and "3" when the three-way and two-way causal interactions are of interest. Default is 2.}
  \item{diff}{a logical indicating whether the outcome is the choice between a pair. If \code{diff=TRUE}, \code{pair.id} should specify a pair of comparison.}
  \item{select.prob}{a logical indicating whether selection probabilities are computed. We recommend to use this option when the model is finalized.}
  \item{boot}{the number of bootstrap replicates.}
  \item{block.id}{unique identifies for blocks that the function uses for block bootstrap.}
  \item{seed}{seed for bootstrap.}
  \item{eps}{a tolerance parameter in the internal optimization algorithm.}
  \item{fac.level}{optional. A vector containing the number of levels in each factor. The order of \code{fac.level} should match to the order of columns in the data. For example, when the first and second columns of the design matrix is "Education" and "Race", the first and second element of \code{fac.level} should be the number of levels in "Education" and "Race", respectively.}
  \item{ord.fac}{optional. logical vectors indicating whether each factor has ordered (\code{TRUE}) or unordered (\code{FALSE}) levels. When levels are ordered, the function uses the order given by function \code{levels()}. If levels are ordered, the function places penalties on the differences between adjacent levels. If levels are unordered, the function places penalties on the differences based on every pairwise comparison.}
  \item{verbose}{whether it prints the value of a cost parameter used.}
}

\details{Suggested workflow.
\enumerate{
\item{Specify the order of levels within each factor using \code{levels()}.}{   Since the function places penalties on the differences between adjacent levels when levels are ordered, it is crucial to specify the order of levels within each factor carefully.}
\item{Run \code{cv.CausalANOVA}.}{   Select the cost parameter minimizing the cross-validation error. Or choose largest value of \code{cost} such that error is within 1 standard error of the minimum. \code{plot.cv.CausalANOVA} can be used to investigate how cross-validation errors vary depending on cost parameters.}
\item{Run \code{CausalANOVA}.}{   Run the main model with the chosen cost parameter and see summary by \code{summary.CausalANOVA}. If researchers want to compute selection probabilities, set \code{select.prob=TRUE}. Given it is computationally intensive, we recommend to compute selection probabilities when the model is finalized. The selection probability for the range of the AME (AMIE) is one minus the proportion of bootstrap replicates in which all coefficients for the corresponding factor (factor interaction) are estimated to be zero. The selection probability of the AME (AMIE) is the proportion of bootstrap replicates in which the sign of the effect is the same as the point estimate.}
\item{Investigate two-way interactions.}{   Run \code{plot.CausalANOVA} and visualize the AMIEs by choosing two factors of interest. Run \code{AMIE} to examine decomposition of the average combination effect into the AMIE and AMEs.}
}
}
\value{
  \item{intercept}{an intercept of the estimated ANOVA model.If \code{diff=TRUE}, this should be close to 0.5.}
  \item{coefs}{a named vector of coefficients of the estimated ANOVA model.}
  \item{formula}{the \code{formula} used in the function.}
  \item{cost}{the cost parameter used in the function.}
  \item{...}{arguments passed to the function or arguments only for the internal use.}	
}
\references{
Post, J. B. and Bondell, H. D. 2013. ``Factor selection and structural identification
in the interaction anova model.'' Biometrics 69, 1, 70--79.

Egami, Naoki and Kosuke Imai. 2016+. ``Causal Interaction in Factorial Experiments: Application to Conjoint Analysis.'' Working paper. \url{http://imai.princeton.edu/research/files/int.pdf}
}

\author{Naoki Egami and Kosuke Imai.}

\seealso{\link{cv.CausalANOVA},\link{AMIE}}

\examples{
data(Carlson)
## Specify the order of each factor
Carlson$newRecordF<- factor(Carlson$newRecordF,ordered=TRUE,
                         levels=c("YesLC", "YesDis","YesMP",
                             "noLC","noDis","noMP","noBusi"))
Carlson$promise <- factor(Carlson$promise,ordered=TRUE,levels=c("jobs","clinic","education"))
Carlson$coeth_voting <- factor(Carlson$coeth_voting,ordered=FALSE,levels=c("0","1"))
Carlson$relevantdegree <- factor(Carlson$relevantdegree,ordered=FALSE,levels=c("0","1"))

## Run cv.CausalANOVA
\dontrun{
cv.fit <- cv.CausalANOVA(won ~ newRecordF + promise + coeth_voting + relevantdegree,
                         data=Carlson,
                         pair.id=Carlson$contestresp,diff=TRUE, nway=2)

cv.fit
plot(cv.fit)
}

fit <- CausalANOVA(won ~ newRecordF + promise + coeth_voting + relevantdegree,
                    data=Carlson,
                    pair.id=Carlson$contestresp,diff=TRUE, nway=2,cost=0.15)
## Or when we need selection probabilities.
\dontrun{
fit <- CausalANOVA(won ~ newRecordF + promise + coeth_voting + relevantdegree,
                    data=Carlson,
                    pair.id=Carlson$contestresp,diff=TRUE,nway=2,cost=0.15,
                    select.prob=TRUE,boot=500,block.id=Carlson$respcodeS)
}
summary(fit)

\dontrun{
## plot 
plot(fit,fac.name=c("newRecordF","coeth_voting"))
}

## compute AMIEs
amie1 <- AMIE(fit,fac.name=c("promise","newRecordF"),
              level.name=c("jobs","noLC"),
              base.name=c("jobs","YesLC"))

amie2 <- AMIE(fit,fac.name=c("newRecordF","coeth_voting"),
              level.name=c("noBus","1"),
              base.name=c("noMP","0"))

}