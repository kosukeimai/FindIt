\name{AMIE}
\alias{AMIE}
\title{Decomposing the Combination Effect into the AMEs and the AMIE.}
\description{
  \code{AMIE} shows the decomposition of the average combination effect into the AMEs and the AMIE. The function can decompose the conditional effect by setting the level of one factor equal to its baseline.
}
\usage{	
	AMIE(object,fac.name,level.name,base.name=NULL,verbose=TRUE,...)
}
\arguments{
  \item{object}{an object from \code{CausalANOVA}.}	
  \item{fac.name}{a character 2 dimensional vector indicating two factors of interest.}
  \item{level.name}{a character 2 dimensional vector indicating levels of each factor. When one of \code{level.name} is equal to \code{base.name}, the function decomposes the conditional effect instead of the combination effect.}
  \item{base.name}{a character 2 dimensional vector indicating baseline levels of each factor. If not specified, the function uses the last level defined by \code{levels()}.}	
  \item{verbose}{whether it prints the results.}
  \item{...}{arguments passed to the function or arguments only for the internal use.}	
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
  \item{AMIE.main}{an estimated AMIE.}
  \item{baseline}{baseline levels in each factor.}
  \item{decompose}{the decomposition of the combination effect into the AMEs and the AMIE. When one of \code{level.name} is equal to \code{base.name}, the function decomposes the conditional effect.}
}
\references{
Post, J. B. and Bondell, H. D. 2013. ``Factor selection and structural identification
in the interaction anova model.'' Biometrics 69, 1, 70--79.

Egami, Naoki and Kosuke Imai. 2016+. ``Causal Interaction in Factorial Experiments: Application to Conjoint Analysis.'' Working paper. \url{http://imai.princeton.edu/research/files/int.pdf}
}

\author{Naoki Egami and Kosuke Imai.}

\seealso{\link{cv.CausalANOVA},\link{CausalANOVA}}

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