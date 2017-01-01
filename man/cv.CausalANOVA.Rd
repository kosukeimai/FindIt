\name{cv.CausalANOVA}
\alias{cv.CausalANOVA}
\alias{plot.cv.CausalANOVA}
\title{Cross validation for the CausalANOVA.}
\description{
  \code{cv.CausalANOVA} implements cross-validation for \code{CausalANOVA} to select the cost parameter. 
}
\usage{
cv.CausalANOVA(formula,data,cv.cost=c(0.1,0.3,0.5,0.7,1.0), type="bin",
               pair.id=NULL,nway=2,diff=TRUE,eps=1e-5,nfolds=10,seed=1234)
}
\arguments{
  \item{formula}{a formula that specifies outcome and treatment variables.}	
  \item{data}{an optional data frame, list or environment (or object
              coercible by 'as.data.frame' to a data frame) containing the
              variables in the model.  If not found in 'data', the
              variables are taken from 'environment(formula)', typically
              the environment from which 'CausalANOVA' is called.}
  \item{cv.cost}{a vector containing candidates for a cost parameter ranging from 0 to 1. 1 corresponds to no regularization and the smaller value corresponds to the stronger regularization. Default is \code{c(0.1,0.3,0.5,0.7,1.0)}.}
  \item{type}{When the outcome is binary, set \code{type} to \code{"bin"}. Cross-validation error is based on misclassification. WHen the outcome is continuous, set \code{type} to \code{"cont"}. Cross-validation error is based on the mean squared error.}
\item{pair.id}{unique identifiers for each pair of comparison. This option is used when \code{dif=TRUE}.}
  \item{nway}{"2" when the two way causal interactions are of interest and "3" when the three-way and two-way causal interactions are of interest. Default is 2.}
  \item{diff}{a logical indicating whether the outcome is the choice between a pair. If \code{diff=TRUE}, \code{pair.id} should specify a pair of comparison.}	
  \item{eps}{a tolerance parameter in the internal optimization algorithm.}
  \item{nfolds}{number of folds - default is 10. Although nfolds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets.}
  \item{seed}{an argument for \code{set.seed()}.}
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
  \item{cv.error}{The mean cross-validated error - a vector of length \code{length(cv.t)}.}
  \item{cv.min}{value of \code{t} that gives minimum \code{cv.missclass}.}
  \item{cv.1sd}{largest value of \code{t} such that error is within 1 standard error of the minimum.}
  \item{cv.each.mat}{a matrix containing cross-validation errors for each fold and cost parameter.}
  \item{cv.cost}{the \code{cv.cost} used in the function.}
}
\references{
Post, J. B. and Bondell, H. D. 2013. ``Factor selection and structural identification
in the interaction anova model.'' Biometrics 69, 1, 70--79.

Egami, Naoki and Kosuke Imai. 2016+. ``Causal Interaction in Factorial Experiments: Application to Conjoint Analysis.'' Working paper. \url{http://imai.princeton.edu/research/files/int.pdf}
}

\author{Naoki Egami and Kosuke Imai.}

\seealso{\link{CausalANOVA},\link{AMIE}.}

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