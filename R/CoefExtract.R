CoefExtract <- function(BaseCoef, base.name, fac.level, ord.fac, Gorder, indTwo=NULL, indThree=NULL){
    ## input: BaseCoef
    ## I need to remove intercept first 
    
    n.fac <- length(fac.level)   
    levelIndex <- CreatelevelIndex(fac.level=fac.level,ord.fac=ord.fac,Gorder=Gorder,
                                   indTwo=indTwo, indThree=indThree)
    levelIndex$end <- cumsum(levelIndex$length)
    use.ind <- levelIndex$dif==0
    coefIndex <- levelIndex[use.ind==1,]
    plus.ind <- coefIndex$plus==1
    coefIndexPlus <- coefIndex[plus.ind==1,]
    coefIndexMinus <- coefIndex[plus.ind==0,]
    nameIndex <- coefIndexPlus
    nameIndex$start <- c(1, cumsum(nameIndex$length)[-nrow(nameIndex)]+1)
    nameIndex$end   <- cumsum(nameIndex$length)
    
    coefPlus <- c()
    Plus.list <- list()
    for(i in 1:nrow(coefIndexPlus)){
        ind <- coefIndexPlus$start[i]:coefIndexPlus$end[i]        
        Plus.list[[i]] <- BaseCoef[ind]       
        coefPlus <- c(coefPlus,BaseCoef[ind])
    }
    coefMinus <- c()
    Minus.list <- list()
    for(i in 1:nrow(coefIndexMinus)){
        ind <- coefIndexMinus$start[i]:coefIndexMinus$end[i]
        Minus.list[[i]] <- BaseCoef[ind]
        coefMinus <- c(coefMinus,BaseCoef[ind])
    }
    
    Coef.list <- list()
    for(i in 1:nrow(coefIndexPlus)){
        name.ind <- nameIndex$start[i]:nameIndex$end[i]
        coef <- Plus.list[[i]] - Minus.list[[i]]
        names(coef) <- base.name[name.ind]
        Coef.list[[i]] <- coef
    }    
    return(Coef.list)
}
