##################################################
## Naoki Egami
## 2016/08/22
##################################################

## Even for three-way interaction models, I can exactly the same funciton
## for Zonefunction, Ztwofunction


## I need to make this general for the general formula 
Zcombinefunction <- function(X,fac.level,ord.fac,Gorder){

    
    Zonefunction <-function(Xone,nlevel,ord=FALSE){
        if(ord==FALSE){
            qp <- nlevel*(nlevel-1)/2
        }else if(ord==TRUE){
            qp <- (nlevel-1)
        }
        zero.mat <- matrix(0,nrow=nrow(Xone),ncol=(2*qp))
        Zone <- cbind(Xone, -Xone, zero.mat)
        return(Zone)
    }


    Ztwofunction <- function(Xtwo,alevel,blevel,ordAB=c(FALSE,FALSE)){
        if(ordAB[1]==FALSE){
            r1 <- (blevel)*(alevel)*(alevel-1)/2
        }else if(ordAB[1]==TRUE){
            r1 <- (blevel)*(alevel-1)
        }
        if(ordAB[2]==FALSE){
            r2 <- (alevel)*(blevel)*(blevel-1)/2
        }else if(ordAB[2]==TRUE){
            r2 <- (alevel)*(blevel-1)
        }    
        zero.mat1 <- matrix(0,nrow=nrow(Xtwo),ncol=(2*r1))
        zero.mat2 <- matrix(0,nrow=nrow(Xtwo),ncol=(2*r2))
        Ztwo <- cbind(Xtwo, -Xtwo, zero.mat1, zero.mat2)
        return(Ztwo)
    }

    Zthreefunction <- function(Xthree,alevel,blevel,clevel,ordABC=c(FALSE,FALSE,FALSE)){
        if(ordABC[1]==FALSE){
            r1 <- (blevel*clevel)*(alevel)*(alevel-1)/2
        }else if(ordABC[1]==TRUE){
            r1 <- (blevel*clevel)*(alevel-1)
        }
        if(ordABC[2]==FALSE){
            r2 <- (alevel*clevel)*(blevel)*(blevel-1)/2
        }else if(ordABC[2]==TRUE){
            r2 <- (alevel*clevel)*(blevel-1)
        }
        if(ordABC[3]==FALSE){
            r3 <- (alevel*blevel)*(clevel)*(clevel-1)/2
        }else if(ordABC[3]==TRUE){
            r3 <- (alevel*blevel)*(clevel-1)
        }
        zero.mat1 <- matrix(0,nrow=nrow(Xthree),ncol=(2*r1))
        zero.mat2 <- matrix(0,nrow=nrow(Xthree),ncol=(2*r2))
        zero.mat3 <- matrix(0,nrow=nrow(Xthree),ncol=(2*r3))
        Zthree <- cbind(Xthree, -Xthree, zero.mat1, zero.mat2, zero.mat3)
        return(Zthree)
    }
    ## Input should not have a intercept
    ## Add intercept later
    ## This does not change with ordered categorical variables.

    ## fac.level <- c(4,7,4)
    ## ord <- c(TRUE,TRUE,TRUE)

    levelIndex <- CreatelevelIndex(fac.level=fac.level,ord.fac=ord.fac,Gorder=Gorder)
    use.ind <- (levelIndex$plus==1) * (levelIndex$dif==0)
    Index.use <- levelIndex[use.ind==1,]
    Index.use$start <- c(1,cumsum(Index.use$length)[-nrow(Index.use)]+1)
    Index.use$end <-  cumsum(Index.use$length)

    Fac.index <- Index.use[,regexpr("Fac",colnames(Index.use))>0]

    Z.list <- list()

    for(i in 1:nrow(Index.use)){
        if(Index.use$order[i]==1){
            Xone <- X[,Index.use$start[i]:Index.use$end[i]]
            level.main <- fac.level[which(Fac.index[i,]==1)]
            ord.main <- ord.fac[which(Fac.index[i,]==1)]
            Z.list[[i]] <- Zonefunction(Xone=Xone,nlevel=level.main,ord=ord.main)
        }else if(Index.use$order[i]==2){
            Xtwo <- X[,Index.use$start[i]:Index.use$end[i]]
            level.main <- fac.level[which(Fac.index[i,]==1)]
            ord.main <- ord.fac[which(Fac.index[i,]==1)]
            Z.list[[i]] <- Ztwofunction(Xtwo=Xtwo,alevel=level.main[1],
                                        blevel=level.main[2],ordAB=ord.main)
        }else if(Index.use$order[i]==3){
            Xthree <- X[,Index.use$start[i]:Index.use$end[i]]
            level.main <- fac.level[which(Fac.index[i,]==1)]
            ord.main <- ord.fac[which(Fac.index[i,]==1)]
            Z.list[[i]] <- Zthreefunction(Xthree=Xthree,
                                          alevel=level.main[1],blevel=level.main[2],
                                          clevel=level.main[3],
                                          ordABC=ord.main)
        }
    }

    Zcombine.temp <- do.call(cbind,Z.list)
    Zcombine <- cbind(1, Zcombine.temp)
    
    return(Zcombine)
}

