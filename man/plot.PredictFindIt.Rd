\name{plot.PredictFindIt}
\alias{plot.PredictFindIt}
\title{
Plot estimated treatment effects or predicted outcomes for each treatment combination.
}
\description{
Plot estimated treatment effects when \code{treat.type="single"} and predicted outcomes for each treatment combination when \code{treat.type="multiple"}.
}
\usage{ 
\method{plot}{PredictFindIt}(x,main,xlab,ylab, ...)
}
\arguments{
  \item{x}{
  output from \code{predict.FindIt}.
}
  \item{main}{
  the argument specifying the main title of the plot.
}
  \item{xlab}{
  the argument specifying the name of x axis.
}
  \item{ylab}{
  the argument specifying the name of y axis.
}
  \item{\dots}{
  further arguments passed to or from other methods.
}
}
\details{
	Plot estimated treatment effects when \code{treat.type="single"} and predicted outcomes for each treatment combination when \code{treat.type="multiple"}.
}
\value{
   \item{plot}{Plot estimated treatment effects when \code{treat.type="single"} and predicted outcomes for each treatment combination when \code{treat.type="multiple"}.}
}

\author{Naoki Egami, Marc Ratkovic and Kosuke Imai.}

%\seealso{\link{FindIt}}

\examples{
## See the help page for FindIt() for an example.
	
}
