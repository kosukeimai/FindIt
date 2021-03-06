#' Estimating the Conditional Effects with the CausalANOVA.
#' 
#' \code{ConditionalEffect} estimates a variety of conditional effects using
#' the ouput from \code{CausalANOVA}.
#' 
#' See Details in \code{CausalANOVA}.
#' 
#' @param object The output from \code{CausalANOAV} function.
#' @param treat.fac The name of factor acting as the main treatment variable.
#' @param cond.fac The name of factor acting as the conditioning (moderating)
#' variable.
#' @param base.ind An indicator for the baseline of the treatment factor.
#' Default is 1.
#' @param round Digits to round estimates. Default is 3.
#' @param inference (optional). This argument is mainly for internal use. It
#' indicates whether \code{CausalANOVA} has done inference or not.
#' @param verbose Whether it prints the progress.
#' @return \item{CondtionalEffects}{The summary of estimated conditional
#' effects.} \item{...}{Arguments for the internal use.}
#' @author Naoki Egami and Kosuke Imai.
#' @seealso \link{CausalANOVA}.
#' @references Egami, Naoki and Kosuke Imai. 2019. Causal Interaction in
#' Factorial Experiments: Application to Conjoint Analysis, Journal of the American Statistical Association.
#' \url{http://imai.fas.harvard.edu/research/files/int.pdf}
#' 
#' Lim, M. and Hastie, T. 2015. Learning interactions via hierarchical
#' group-lasso regularization. Journal of Computational and Graphical
#' Statistics 24, 3, 627--654.
#' 
#' Post, J. B. and Bondell, H. D. 2013. ``Factor selection and structural
#' identification in the interaction anova model.'' Biometrics 69, 1, 70--79.
#' @examples
#' 
#' data(Carlson)
#' ## Specify the order of each factor
#' Carlson$newRecordF<- factor(Carlson$newRecordF,ordered=TRUE,
#'                             levels=c("YesLC", "YesDis","YesMP",
#'                                      "noLC","noDis","noMP","noBusi"))
#' Carlson$promise <- factor(Carlson$promise,ordered=TRUE,levels=c("jobs","clinic","education"))
#' Carlson$coeth_voting <- factor(Carlson$coeth_voting,ordered=FALSE,levels=c("0","1"))
#' Carlson$relevantdegree <- factor(Carlson$relevantdegree,ordered=FALSE,levels=c("0","1"))
#' 
#' ## ####################################### 
#' ## Without Screening and Collapsing
#' ## ####################################### 
#' #################### AMEs and two-way AMIEs ####################
#' fit2 <- CausalANOVA(formula=won ~ newRecordF + promise + coeth_voting + relevantdegree,
#'                     int2.formula = ~ newRecordF:coeth_voting,
#'                     data=Carlson, pair.id=Carlson$contestresp,diff=TRUE,
#'                     cluster=Carlson$respcodeS, nway=2)
#' summary(fit2)
#' plot(fit2, type="ConditionalEffect", fac.name=c("newRecordF","coeth_voting"))
#' ConditionalEffect(fit2, treat.fac="newRecordF", cond.fac="coeth_voting")
#' @export
ConditionalEffect <- function(object,treat.fac=NULL, cond.fac=NULL,
                              base.ind=1, round=3,
                              inference=NULL, verbose=TRUE){
    

    if(is.null(inference)==TRUE){
        inference <- object$inference
    }
    ## House Keeping
    ## only for Two-ways
    ## check mod
    ## check base.ind
    ## If selected interactions are not in indTwo.
    
    formula <- object$formula
    AME <- object$AME
    AMIE2 <- object$AMIE2
    fac.level <- object$fac.level
    ord.fac <- object$ord.fac
    Gorder <- object$Gorder
    indTwo <- object$indTwo
    indThree <- object$indThree
    data.main <- object$data[,-1]
    level.name <- lapply(data.main, levels)
    terms.f <- terms(formula,data=object$data)
    order.f <- attr(terms.f, "order")
    var.name <- attr(terms.f, "term.labels")[order.f==1]

    coefs.u <- unlist(object$coefs)
    
    if(inference==TRUE){
        vcov.u <- unlist(object$vcov)
    }

    ## LevelIndex
    ## Index
    levelIndex <- CreatelevelIndex(fac.level=fac.level, ord.fac=ord.fac, Gorder=Gorder,
                                   indTwo=indTwo, indThree=indThree)
    use.ind <-  (levelIndex$plus==1)*(levelIndex$dif==0)
    Index.use <- levelIndex[use.ind==1,]
    Index.use$start <- c(1,cumsum(Index.use$length)[-nrow(Index.use)]+1)
    Index.use$end <-  cumsum(Index.use$length)
    
    
    ## Find Index for Factors 
    treat.ind <- which(treat.fac==var.name)
    cond.ind <- which(cond.fac==var.name)
    treat.level <- levels(data.main[,treat.ind])
    cond.level <- levels(data.main[, cond.ind])
    rm(data.main)

    if(is.null(base.ind)==TRUE){
        base.ind <- length(treat.level)
    }
    
    Fac.ind <- c(cond.ind, treat.ind)
    Norotate <- all(order(Fac.ind) == c(1,2))
    Fac.ind <- Fac.ind[order(Fac.ind)]
    var.ind.mat <- matrix(seq(1:(fac.level[Fac.ind[1]]*fac.level[Fac.ind[2]])),
                          nrow=fac.level[Fac.ind[1]], ncol=fac.level[Fac.ind[2]])
    if(Norotate==FALSE){ var.ind.mat <- t(var.ind.mat) }


    ## Find the index for Main Effect
    AME.var.ind <- Index.use$start[treat.ind]:Index.use$end[treat.ind]
    AME.var.ind.mat <- as.matrix(cbind(AME.var.ind, AME.var.ind[base.ind]))
    ## rownames(AME.var.ind.mat) <- level.name[[treat.ind]]
    AME.var.ind.mat.f <- rep(list(AME.var.ind.mat), fac.level[cond.ind])
    ## names(AME.var.ind.mat.f) <- paste(var.name[cond.ind], "=", level.name[[cond.ind]],sep="")
    
    ## Find the index for Interaction 
    INT.ind <- which(apply(c(Fac.ind[1], Fac.ind[2]) == indTwo, 2, all)==TRUE)
    if(length(INT.ind)==0) stop("Specified Interactions are not in the model.")
    INT.ind.u <- sum(order.f==1) + INT.ind ## z
    AMIE.var.ind <- Index.use$start[INT.ind.u]:Index.use$end[INT.ind.u]
    
    ind.var <- list()
    for(i in 1:nrow(var.ind.mat)){
        ind.var0 <- matrix(NA, ncol=2, nrow=ncol(var.ind.mat))
        for(j in 1:ncol(var.ind.mat)){
            ind.var0[j, 1:2] <- c(var.ind.mat[i,j], var.ind.mat[i,base.ind])
        }
        ind.var[[i]] <- ind.var0
    }
    ## ## always row=Conditional, col=Treatment
    AMIE.var.ind.mat.f <- list()
    for(i in 1:length(ind.var)){
        AMIE.var.ind.mat.f[[i]] <- cbind(AMIE.var.ind[ind.var[[i]][,1]],AMIE.var.ind[ind.var[[i]][,2]])
    }

    ## Compute Conditional Effect
    if(inference==TRUE){
        CE.l <- list()
        for(i in 1:length(ind.var)){
            CE.tab <- matrix(NA, ncol=4, nrow=nrow(AMIE.var.ind.mat.f[[i]]))
            for(ce in 1:nrow(AMIE.var.ind.mat.f[[i]])){
                point <- (coefs.u[AME.var.ind.mat.f[[i]][ce,1]] - coefs.u[AME.var.ind.mat.f[[i]][ce,2]]) +
                    (coefs.u[AMIE.var.ind.mat.f[[i]][ce,1]] - coefs.u[AMIE.var.ind.mat.f[[i]][ce,2]]) 
                vari <- VarCondEffect(AME.var.ind.mat.f[[i]][ce,1], AME.var.ind.mat.f[[i]][ce,2],
                                      AMIE.var.ind.mat.f[[i]][ce,1], AMIE.var.ind.mat.f[[i]][ce,2], vcov.full=vcov.u)
                if(vari>0) std <- sqrt(vari) else std <- 0 
                CE.tab[ce,1:4] <- c(point, std, point - 1.96*std, point + 1.96*std)
            }
            rownames(CE.tab) <- level.name[[treat.ind]]
            colnames(CE.tab) <- c("ConditionalEffect", "sd", "2.5% CI", "97.5% CI")
            CE.tab <- round(CE.tab, digits=round)
            CE.l[[i]] <- CE.tab
        }
        names(CE.l) <- paste(var.name[cond.ind], "=", level.name[[cond.ind]],sep="")
    }else{
        CE.l <- list()
        for(i in 1:length(ind.var)){
            CE.tab <- c()
            for(ce in 1:nrow(AMIE.var.ind.mat.f[[i]])){
                point <- (coefs.u[AME.var.ind.mat.f[[i]][ce,1]] - coefs.u[AME.var.ind.mat.f[[i]][ce,2]]) +
                    (coefs.u[AMIE.var.ind.mat.f[[i]][ce,1]] - coefs.u[AMIE.var.ind.mat.f[[i]][ce,2]])                 
                CE.tab[ce] <- point
            }
            names(CE.tab) <- level.name[[treat.ind]]
            CE.tab <- round(CE.tab, digits=round)
            CE.l[[i]] <- CE.tab
        }
        names(CE.l) <- paste(var.name[cond.ind], "=", level.name[[cond.ind]],sep="")
    }
    
    
    ## ## base.ind
    ## ## if(missing(base.ind)==TRUE) base.ind <- 1

    ## AME.u0 <- AME[treat.ind][[1]] - AME[treat.ind][[1]][base.ind]
    ## AME.u <- matrix(rep(AME.u0, fac.level[cond.ind]), nrow=fac.level[cond.ind], byrow=TRUE)
    
    ## AMIE.mat <- matrix(unlist(AMIE2[INT.ind]), nrow=fac.level[Fac.ind[1]], ncol=fac.level[Fac.ind[2]])
    ## if(Norotate==FALSE){
    ##     AMIE.mat <- t(AMIE.mat)
    ## }else{
    ##         
    ## AMIE.u  <- AMIE.mat - AMIE.mat[,base.ind]
    
    ## CE <- AME.u + AMIE.u
    
    if(verbose==TRUE){
        cat(paste("\nTreatment Factor is ", var.name[treat.ind], " and ",
                  "Conditioning Factor is ", var.name[cond.ind],"\n",sep=""))
        print(CE.l)
    }
    
    output <- list("ConditionalEffects"=CE.l,
                   "treat.fac"=treat.fac,
                   "cond.fac"=cond.fac,
                   "treat.level"=treat.level,
                   "cond.level"=cond.level)
}
