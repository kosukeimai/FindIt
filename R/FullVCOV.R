########################################
## Naoki Egami
## 2017/08/02
########################################

## Index to Index Function 


## output: numbers and matrix with 1 
IndexOne <- function(index, n.level){
    index.u <- list()
    base.u <- rep(0, n.level)
    if(index <=n.level - 1){
        index.u[[1]] <- index   ## just an index in the trancated scale 
        base.u[index] <- 1      ## matrix form 
        index.u[[2]] <- base.u  ## matrix 
        index.u[[3]] <- 1  ## plus
    }else{
        index.u[[1]] <- seq(1:(n.level-1))
        base.u[-index] <- -1
        index.u[[2]] <- base.u  ## matrix 
        index.u[[3]] <- -1  ## plus
    }
    return(index.u)
}

IndexTwo <- function(index1, index2, n.level1, n.level2){
    index.u <- list()
    base.u <- matrix(0, nrow=n.level1, ncol=n.level2)
    if(index1 <= n.level1 - 1 & index2 <= n.level2 - 1){
        index.u[[1]] <- index1 + (index2 - 1)*(n.level1 - 1)   ## just an index in the trancated scale 
        base.u[index1, index2] <- 1      ## matrix form 
        index.u[[2]] <- base.u  ## matrix 
        index.u[[3]] <- 1  ## plus
    }else if(index1 == n.level1 & index2 <= n.level2 -1){
        index.u[[1]] <- seq(from=(1 + (index2 - 1)*(n.level1 - 1)), to=(index2*(n.level1 - 1)))
        base.u[1:(n.level1-1), index2] <- -1
        index.u[[2]] <- base.u  ## matrix 
        index.u[[3]] <- -1  ## plus
    }else if(index1 <= (n.level1 -1) & index2 == n.level2){
        index.u[[1]] <- seq(from=index1, to=(index1+(n.level1 - 1)*(n.level2-2)), by=(n.level1-1))
        base.u[index1, 1:(n.level2-1)] <- -1
        index.u[[2]] <- base.u  ## matrix 
        index.u[[3]] <- -1  ## plus
    }else if(index1 == n.level1 & index2 == n.level2){
        index.u[[1]] <- seq(from=1, to=(n.level1 - 1)*(n.level2-1))
        base.u[1:(n.level1-1), 1:(n.level2-1)] <- 1
        index.u[[2]] <- base.u  ## matrix 
        index.u[[3]] <- 1  ## plus
    }
    return(index.u)
}

IndexThree <- function(index1, index2, index3, n.level1, n.level2, n.level3){
    index.u <- list()
    base.u <- array(0, dim=c(n.level1,n.level2, n.level3))
    if(index3 <= n.level3 - 1){
        if(index1 <= n.level1 - 1 & index2 <= n.level2 - 1){
            index.u[[1]] <- index1 + (index2 - 1)*(n.level1 - 1) +
                (n.level1-1)*(n.level2-1)*(index3-1) ## just an index in the trancated scale 
            base.u[index1, index2, index3] <- 1      ## matrix form 
            index.u[[2]] <- base.u  ## matrix 
            index.u[[3]] <- 1  ## plus
        }else if(index1 == n.level1 & index2 <= n.level2 -1){
            index.u[[1]] <- seq(from=(1 + (index2 - 1)*(n.level1 - 1)), to=(index2*(n.level1 - 1))) +
                (n.level1-1)*(n.level2-1)*(index3-1)
            base.u[1:(n.level1-1), index2, index3] <- -1
            index.u[[2]] <- base.u  ## matrix 
            index.u[[3]] <- -1  ## plus
        }else if(index1 <= (n.level1 -1) & index2 == n.level2){
            index.u[[1]] <- seq(from=index1, to=(index1+(n.level1 - 1)*(n.level2-2)), by=(n.level1-1)) +
                (n.level1-1)*(n.level2-1)*(index3-1)
            base.u[index1, 1:(n.level2-1), index3] <- -1
            index.u[[2]] <- base.u  ## matrix 
            index.u[[3]] <- -1  ## plus
        }else if(index1 == n.level1 & index2 == n.level2){
            index.u[[1]] <- seq(from=1, to=(n.level1 - 1)*(n.level2-1)) +
                (n.level1-1)*(n.level2-1)*(index3-1)
            base.u[1:(n.level1-1), 1:(n.level2-1), index3] <- 1
            index.u[[2]] <- base.u  ## matrix 
            index.u[[3]] <- 1  ## plus
        }
    }else if(index3==n.level3){
        if(index1 <= n.level1 - 1 & index2 <= n.level2 - 1){
            index.u[[1]] <- seq(from=(index1 + (index2 - 1)*(n.level1 - 1)),
                                to=((index1 + (index2 - 1)*(n.level1 - 1)) +
                                    (n.level1-1)*(n.level2-1)*(n.level3-2)),
                                by=(n.level1-1)*(n.level2-1))  ## just an index in the trancated scale 
            base.u[index1, index2, seq(1:(n.level3-1))] <- -1      ## matrix form 
            index.u[[2]] <- base.u  ## matrix 
            index.u[[3]] <- -1  ## plus
        }else if(index1 == n.level1 & index2 <= n.level2 -1){
            b.i <- seq(from=(1 + (index2 - 1)*(n.level1 - 1)), to=(index2*(n.level1 - 1)))
            b.i.l <- c()
            for(z in 1:(n.level3-1)){
                b.i.l <- c(b.i.l, (n.level1-1)*(n.level2-1)*(z-1) + b.i)
            }
            index.u[[1]] <- b.i.l
            base.u[1:(n.level1-1), index2, seq(1:(n.level3-1))] <- 1
            index.u[[2]] <- base.u  ## matrix 
            index.u[[3]] <- 1  ## plus
        }else if(index1 <= (n.level1 -1) & index2 == n.level2){
            b.i <- seq(from=index1, to=(index1+(n.level1 - 1)*(n.level2-2)), by=(n.level1-1))
            b.i.l <- c()
            for(z in 1:(n.level3-1)){
                b.i.l <- c(b.i.l, (n.level1-1)*(n.level2-1)*(z-1) + b.i)
            }
            index.u[[1]] <- b.i.l
            base.u[index1, 1:(n.level2-1), seq(1:(n.level3-1))] <- 1
            index.u[[2]] <- base.u  ## matrix 
            index.u[[3]] <- -1  ## plus
        }else if(index1 == n.level1 & index2 == n.level2){
            index.u[[1]] <- seq(from=1, to=(n.level1 - 1)*(n.level2-1)*(n.level3-1))
            base.u[1:(n.level1-1), 1:(n.level2-1), seq(1:(n.level3-1))] <- -1
            index.u[[2]] <- base.u  ## matrix 
            index.u[[3]] <- -1  ## plus
        }
    }
    return(index.u)
}

IndConvert <- function(index.level, index.fac, term.index,
                       fac.level, Index.use){
    nway <- length(index.level)
    fac.use <- fac.level[index.fac]
    if(term.index==1) base.n <- 0 else base.n <- Index.use$end[(term.index-1)]
    if(nway==1) base <- IndexOne(index=index.level, n.level=fac.use)
    if(nway==2) base <- IndexTwo(index1=index.level[1], index2=index.level[2],
                                 n.level1=fac.use[1], n.level2=fac.use[2])
    if(nway==3) base <- IndexThree(index1=index.level[1], index2=index.level[2], index3=index.level[3],
                                   n.level1=fac.use[1], n.level2=fac.use[2], n.level3=fac.use[3])
    base.ind <- base[[1]]
    sign <- base[[3]]
    ## final ind
    index.u <- base.n + base.ind
    out <- list("index.u"=index.u, "sign"=sign)
    return(out)
}

ind2Var <- function(out, vcov){
    ## index.u from IndConvert
    index.u <- out$index.u
    if(length(index.u)==1){
        Var <- vcov[index.u, index.u]
    }else{
        ## Linear Combination
        Var <- sum(vcov[index.u, index.u])
    }
    return(Var)
}

ind2Cov <- function(out1, out2, vcov){
    ## index.u1 and index.u2 from IndConvert
    index.u1 <- out1$index.u
    index.u2 <- out2$index.u
    sign1 <- matrix(rep(out1$sign, length(index.u1)),nrow=1)
    sign2 <- matrix(rep(out2$sign, length(index.u2)),nrow=1)        
    base.COV <- vcov[index.u1, index.u2]
    Cov <- sign1 %*% base.COV %*% t(sign2)
    return(Cov)
}
## Use these to construct FULL Variance Matrix
## I want to fill in based on rows.

VarEffect <- function(ind1, ind2, vcov.full){    
    effect.var <- vcov.full[ind1,ind1] + vcov.full[ind2,ind2] - 2*vcov.full[ind1,ind2]
    return(effect.var)
}

VarCondEffect <- function(ame.ind1, ame.ind2,
                          amie.ind1, amie.ind2,
                          vcov.full){
    if(ame.ind1== ame.ind2 & amie.ind1==amie.ind2){
        effect.var <- 0
    }else{
        effect.var <- vcov.full[ame.ind1,ame.ind1] + vcov.full[ame.ind2,ame.ind2] +
            vcov.full[amie.ind1,amie.ind1] + vcov.full[amie.ind2,amie.ind2] -
            2*vcov.full[ame.ind1,ame.ind2] + 2*vcov.full[ame.ind1,amie.ind1] -
            2*vcov.full[ame.ind1,amie.ind2] - 2*vcov.full[ame.ind2,amie.ind1] +
            2*vcov.full[ame.ind2,amie.ind2] - 2*vcov.full[amie.ind1,amie.ind2]
    }
    return(effect.var)
}

FullVCOV <- function(vcov, fac.level, ord.fac, indTwo, indThree, Gorder=Gorder){
    ## Now I need to conect Effects to Index.
    levelIndexS <- CreatelevelIndex(fac.level=(fac.level-1),ord.fac=ord.fac, 
                                    Gorder=Gorder, indTwo=indTwo, indThree=indThree)
    use.indS <-  (levelIndexS$plus==1)*(levelIndexS$dif==0)
    Index.useS <- levelIndexS[use.indS==1,]
    levelIndex <- CreatelevelIndex(fac.level=fac.level,ord.fac=ord.fac, 
                                   Gorder=Gorder, indTwo=indTwo, indThree=indThree)
    Index.useS$end <-  cumsum(Index.useS$length)
    use.ind <-  (levelIndex$plus==1)*(levelIndex$dif==0)
    Index.use <- levelIndex[use.ind==1,]    
    Fac.index <- levelIndex[,regexpr("Fac",colnames(levelIndex))>0]
    Fac.Ind.use <- Fac.index[use.ind==1,]

    out <- list()
    sign <- c()
    for(z in 1:nrow(Index.use)){
        if(Index.use$order[z]==1){
            out0 <- list()
            sign0 <- c()
            index.fac <- which(Fac.Ind.use[z,]==1)
            for(i in 1:Index.use$length[z]){
                out0[[i]] <- IndConvert(index.level=i, index.fac=index.fac, term.index=z,
                                        fac.level=fac.level, Index.use=Index.useS)
            }
            out <- c(out, out0)
        }
        if(Index.use$order[z]==2){
            out0 <- list()
            sign0 <- c()
            index.fac <- which(Fac.Ind.use[z,]==1)
            exp2 <- as.matrix(expand.grid(seq(1:fac.level[index.fac[1]]), seq(1:fac.level[index.fac[2]])))
            for(i in 1:Index.use$length[z]){
                out0[[i]] <- IndConvert(index.level=c(exp2[i,1:2]), index.fac=c(index.fac), term.index=z,
                                        fac.level=fac.level, Index.use=Index.useS)                
            }
            out <- c(out, out0)
        }             
        if(Index.use$order[z]==3){
            out0 <- list()
            index.fac <- which(Fac.Ind.use[z,]==1)
            exp3 <- as.matrix(expand.grid(seq(1:fac.level[index.fac[1]]),
                                          seq(1:fac.level[index.fac[2]]),
                                          seq(1:fac.level[index.fac[3]])))
            for(i in 1:Index.use$length[z]){
                out0[[i]] <- IndConvert(index.level=c(exp3[i,1:3]), index.fac=c(index.fac), term.index=z,
                                        fac.level=fac.level, Index.use=Index.useS)

            }
            out <- c(out, out0)               
        }
    }
    
    vcov.f <- matrix(0, ncol=sum(Index.use$length), nrow=sum(Index.use$length))
    vcov.f[1,1] <- ind2Var(out=out[[1]], vcov=vcov)
    for(i in 2:length(out)){           
        main.out <- out[[i]]
        for(j in 1:(i-1)){
            vcov.f[i,j] <- ind2Cov(out1=main.out, out2=out[[j]], vcov=vcov)
        }
        vcov.f[i,i] <- ind2Var(out=main.out, vcov=vcov)
    }
    vcov.f2 <- vcov.f
    diag(vcov.f2) <- 0
    vcov.out <- vcov.f + t(vcov.f2)

    return(vcov.out)
}


