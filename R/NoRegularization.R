NoRegularization <- function(y, X.no, base.name, fac.level, ord.fac, Gorder, indTwo, indThree,
                             cluster){

    lm.fit.no <- lm(y ~ X.no)
    intercept <- coef(lm.fit.no)[1]
    coef.no <- coef(lm.fit.no)[-1]
    levelIndex <- CreatelevelIndex(fac.level=(fac.level-1),ord.fac=ord.fac,Gorder=Gorder,
                                   indTwo=indTwo, indThree=indThree)
    use.ind <-  (levelIndex$plus==1)*(levelIndex$dif==0)
    Index.use <- levelIndex[use.ind==1,]
    Fac.index <- levelIndex[,regexpr("Fac",colnames(levelIndex))>0]
    Fac.Ind.use <- Fac.index[use.ind==1,]
    Fac.level.useC <- (rep(1,nrow(Fac.Ind.use)) %x% t((fac.level-1))) * Fac.Ind.use
    Fac.level.use <- (rep(1,nrow(Fac.Ind.use)) %x% t((fac.level))) * Fac.Ind.use
    Index.use$start <- c(1,cumsum(Index.use$length)[-nrow(Index.use)]+1)
    Index.use$end <-  cumsum(Index.use$length)

    if(Gorder==1){
        Index.use$seq.order <- seq(1:sum(Index.use$order==1))
    }else if(Gorder==2){
        Index.use$seq.order <- c(seq(1:sum(Index.use$order==1)),
                                 seq(1:sum(Index.use$order==2)))
    }else if(Gorder==3){
        Index.use$seq.order <- c(seq(1:sum(Index.use$order==1)),
                                 seq(1:sum(Index.use$order==2)),
                                 seq(1:sum(Index.use$order==3)))
    }
    
    coefs.f <- CreateCoef(coef.use=coef.no, 
                          Index.use=Index.use,
                          Fac.level.use=Fac.level.useC,
                          Gorder=Gorder)$full
    
    if(is.null(cluster)==FALSE){
        vcov.c <- cluster_se_glm(model=lm.fit.no, cluster=cluster)
        vcov.u <- vcov.c[-1, -1]
    }else{
        vcov.u <- vcov(lm.fit.no)[-1,-1]
    }
    
    
    vcov.f <- FullVCOV(vcov=vcov.u, fac.level=fac.level, ord.fac=ord.fac,
                       indTwo=indTwo, indThree=indThree, Gorder=Gorder)
    
    ## Name
    levelIndex2 <- CreatelevelIndex(fac.level=fac.level,ord.fac=ord.fac,Gorder=Gorder,
                                    indTwo=indTwo, indThree=indThree)
    use.ind2 <-  (levelIndex2$plus==1)*(levelIndex2$dif==0)
    Index.use2 <- levelIndex2[use.ind2==1,]
    Index.use2$start <- c(1,cumsum(Index.use2$length)[-nrow(Index.use2)]+1)
    Index.use2$end <-  cumsum(Index.use2$length)
    
    coefs.noreg <- list()
    for(i in 1:nrow(Index.use2)){
        coefs.noreg[[i]] <- c(unlist(coefs.f[[i]]))
        names(coefs.noreg[[i]]) <- base.name[Index.use2$start[i]:Index.use2$end[i]]
    }

    output <- list("intercept"=intercept, "coefs"=coefs.noreg, "vcov"=vcov.f)
    
}
