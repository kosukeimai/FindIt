#' @export
plot.cv.CausalANOVA <- function(x,main=NULL,xlab=NULL,ylab=NULL,...){
    object <- x
    rm(x)
    cv.cost <- object$cv.cost
    cv.error <- object$cv.error
    cv.min.id <- which.min(cv.error)
    cv.sd.each <- apply(object$cv.each.mat,2,sd)
    cv.sd.upper <- cv.error + cv.sd.each
    cv.sd.lower <- cv.error - cv.sd.each
    y.max <- max(cv.sd.upper)
    y.min <- min(cv.sd.lower)

    plot(cv.cost,cv.error,pch=19,ylim=c(y.min,y.max),main=main,xlab="cv.t",
         ylab="Cross-validation error")
    arrows(cv.cost,cv.sd.upper,cv.cost,cv.sd.lower,length=0.05,angle=90,code=3)
    abline(h=cv.sd.upper[cv.min.id],lty=2)
}

