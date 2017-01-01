##################################################
## Naoki Egami
## 2016/08/22
##################################################


Lcombinefunction <- function(fac.level,ord.fac,facCons,Gorder=3){

    ## Make the difference in three-way interactions
    D3function <- function(alevel,blevel,clevel,type="A",ordABC=c(FALSE,FALSE,FALSE)){
        if(type=="A"){
            D3 <- D3Afunction(alevel=alevel,blevel=blevel,clevel=clevel,ord=ordABC[1])
        }else if(type=="B"){
            D3 <- D3Bfunction(alevel=alevel,blevel=blevel,clevel=clevel,ord=ordABC[2])
        }else if(type=="C"){
            D3 <- D3Cfunction(alevel=alevel,blevel=blevel,clevel=clevel,ord=ordABC[3])
        }
        return(D3)
    }

    D3Afunction <- function(alevel,blevel,clevel,ord=FALSE){
        D3a <- matrix(0,nrow=0,ncol=(alevel*blevel*clevel))  
        if(ord==FALSE){
            for(i in 1:(alevel-1)){
                for(j  in (i+1):alevel){
                    base.v <- rep(0,alevel)
                    base.v[i] <- -1
                    base.v[j] <- 1
                    base.mat <- diag((blevel*clevel)) %x% t(base.v)
                    D3a <- rbind(D3a,base.mat)
                }
            }
        }else if(ord==TRUE){
            for(i in 1:(alevel-1)){
                base.v <- rep(0,alevel)
                base.v[i] <- -1
                base.v[(i+1)] <- 1
                base.mat <- diag((blevel*clevel)) %x% t(base.v)
                D3a <- rbind(D3a,base.mat)
            }        
        }
        return(D3a)
    }

    ## (B): Make the difference in two-way interactions 
    D3Bfunction <- function(alevel,blevel,clevel,ord=FALSE){    
        if(ord==FALSE){
            D3b <- matrix(0,nrow=0,ncol=(alevel*blevel*clevel))
            for(i in 1:(blevel-1)){     
                left.temp <- matrix(0,ncol=((i-1)*(alevel)),nrow=(alevel))
                left.main <- -diag((alevel))
                left.final <- cbind(left.temp,left.main)
                for(j  in (i+1):blevel){
                    middle.temp <-  matrix(0,ncol=((j-i-1)*(alevel)),nrow=(alevel))
                    right.main <- diag((alevel))
                    right.temp <- matrix(0,ncol=((blevel-j)*(alevel)),nrow=(alevel))
                    right.final <- cbind(middle.temp,right.main,right.temp)
                    D2b <- cbind(left.final,right.final) ## The base
                    D3base <- diag(clevel) %x% D2b
                    D3b <- rbind(D3b, D3base)
                }
            }
        }else if (ord==TRUE){
            D3b <- matrix(0,nrow=0,ncol=(alevel*blevel*clevel))
            for(i in 1:(blevel-1)){     
                left <- matrix(0,ncol=((i-1)*(alevel)),nrow=(alevel))
                Main <- cbind(-diag((alevel)),diag((alevel)))
                right <- matrix(0,ncol=((blevel-(i+1))*(alevel)),nrow=(alevel))
                D2b <- cbind(left,Main,right)
                D3base <- diag(clevel) %x% D2b
                D3b <- rbind(D3b, D3base)
            }
        }
        return(D3b)
        ## psy_B(2,3), psy_B(2,4), psy_B(2,5) 
    }

    D3Cfunction <- function(alevel,blevel,clevel,ord=FALSE){
        if(ord==FALSE){
            D3c <- matrix(0,nrow=0,ncol=(alevel*blevel*clevel))
            for(i in 1:(clevel-1)){     
                left.temp <- matrix(0,ncol=((i-1)*(alevel*blevel)),nrow=(alevel*blevel))
                left.main <- -diag((alevel*blevel))
                left.final <- cbind(left.temp,left.main)
                for(j  in (i+1):clevel){
                    middle.temp <-  matrix(0,ncol=((j-i-1)*(alevel*blevel)),nrow=(alevel*blevel))
                    right.main <- diag((alevel*blevel))
                    right.temp <- matrix(0,ncol=((clevel-j)*(alevel*blevel)),
                                         nrow=(alevel*blevel))
                    right.final <- cbind(middle.temp,right.main,right.temp)
                    D3c <- rbind(D3c,cbind(left.final,right.final))
                }
            }
        }else if (ord==TRUE){
            Left <- cbind(-diag((alevel*blevel)*(clevel-1)),
                          matrix(0,ncol=(alevel*blevel),nrow=(alevel*blevel)*(clevel-1)))
            Right <- cbind(matrix(0,ncol=(alevel*blevel),nrow=(alevel*blevel)*(clevel-1)),
                           diag((alevel*blevel)*(clevel-1)))
            D3c <- Left + Right
        }
        return(D3c)
    }


    ## From this version, I use D2 rather than D3 
    ## Make the difference in two-way interactions 
    D2function <- function(alevel,blevel,type,ordAB){
        if(type=="A"){
            D2 <- D2Afunction(alevel=alevel,blevel=blevel,ord=ordAB[1])
        }else if(type=="B"){        
            D2 <- D2Bfunction(alevel=alevel,blevel=blevel,ord=ordAB[2])
        }
        return(D2)
        ## D2 matrix
        ## nrow = differences between two-way interactions fixing the sublevel
        ## ncol = alevel*blevel (coefficients for two-ways)
    }

    ## (A): Make the difference in two-way interactions 
    D2Afunction <- function(alevel,blevel,ord=FALSE){
        D2a <- matrix(0,nrow=0,ncol=(alevel*blevel))
        if(ord==FALSE){
            for(i in 1:(alevel-1)){
                for(j  in (i+1):alevel){
                    base.v <- rep(0,alevel)
                    base.v[i] <- -1
                    base.v[j] <- 1
                    base.mat <- diag(blevel) %x% t(base.v)
                    D2a <- rbind(D2a,base.mat)
                }
            }
        }else if(ord==TRUE){
            for(i in 1:(alevel-1)){
                base.v <- rep(0,alevel)
                base.v[i] <- -1
                base.v[(i+1)] <- 1
                base.mat <- diag(blevel) %x% t(base.v)
                D2a <- rbind(D2a,base.mat)
            }        
        }
        return(D2a)
    }

    ## (B): Make the difference in two-way interactions 
    D2Bfunction <- function(alevel,blevel,ord=FALSE){    
        if(ord==FALSE){
            D2b <- matrix(0,nrow=0,ncol=(alevel*blevel))
            for(i in 1:(blevel-1)){     
                left.temp <- matrix(0,ncol=((i-1)*(alevel)),nrow=(alevel))
                left.main <- -diag((alevel))
                left.final <- cbind(left.temp,left.main)
                for(j  in (i+1):blevel){
                    middle.temp <-  matrix(0,ncol=((j-i-1)*(alevel)),nrow=(alevel))
                    right.main <- diag((alevel))
                    right.temp <- matrix(0,ncol=((blevel-j)*(alevel)),nrow=(alevel))
                    right.final <- cbind(middle.temp,right.main,right.temp)
                    D2b <- rbind(D2b,cbind(left.final,right.final))
                }
            }
        }else if (ord==TRUE){
            Left <- cbind(-diag((alevel)*(blevel-1)),matrix(0,ncol=alevel,nrow=(alevel)*(blevel-1)))
            Right <- cbind(matrix(0,ncol=alevel,nrow=(alevel)*(blevel-1)), diag((alevel)*(blevel-1)))
            D2b <- Left + Right
        }
        return(D2b)
        ## psy_B(2,3), psy_B(2,4), psy_B(2,5) 
    }

    ## (Main): Making the difference between all levels
    D1function <- function(nlevel,ord=FALSE){
        if(ord==FALSE){
            qp = (nlevel)*(nlevel-1)/2
            D1.mat <- matrix(0,ncol=(nlevel),nrow=0)    
            if (qp==1) D1.mat <- matrix(c(-1,1),ncol=2,nrow=1)
            if (qp>1){
                for(i in 1 : (nlevel-1)){
                    w.diag <- diag((nlevel-i))
                    left.v <- c(rep(0,(i-1)),-1)
                    left.mat <- matrix(rep(left.v,(nlevel-i)),nrow=(nlevel-i),byrow=TRUE)
                    w.mat <- cbind(left.mat,w.diag)
                    D1.mat <- rbind(D1.mat,w.mat)
                }
            }
        }else if(ord==TRUE){
            if(nlevel==2){
                D1.mat <- matrix(c(-1,1),ncol=2,nrow=1)
            }else{
                D1.mat <- cbind(0,diag(nlevel-1)) + cbind(-diag(nlevel-1),0)
            }
            
        }    
        return(D1.mat)
        ## nrow = differences between main effects
        ## ncol = mainlevel (coefficients for the main factor)
    }

    L1function <- function(nlevel,ord=FALSE){
        ## Creating L for the main effects
        D1 <- D1function(nlevel,ord=ord)
        if(ord==FALSE){
            qp <- (nlevel)*(nlevel-1)/2
        }else if(ord==TRUE){
            qp <- nlevel-1
        }
        L1 <- cbind(D1, -D1, -diag(qp), diag(qp))
        return(L1)
    }

    ## From this version, we call this L2 function.

    L2function <- function(alevel,blevel,ordAB=c(FALSE,FALSE)){
        ## Creating L for interaction effects
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
        D2a <- D2function(alevel=alevel,blevel=blevel,type="A",ordAB=ordAB)
        D2b <- D2function(alevel=alevel,blevel=blevel,type="B",ordAB=ordAB)    
        zero.r1 <- matrix(0,ncol=r1,nrow=r2)
        zero.r2 <- matrix(0,ncol=r2,nrow=r1)
        L2top <- cbind(D2a, -D2a, -diag(r1), diag(r1), zero.r2, zero.r2)
        L2bottom <- cbind(D2b, -D2b, zero.r1, zero.r1, -diag(r2), diag(r2))
        L2 <- rbind(L2top,L2bottom)
        return(L2)
    }

    L3function <- function(alevel,blevel,clevel,ordABC=c(FALSE,FALSE,FALSE)){
        ## Creating L for interaction effects
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
        D3a <- D3function(alevel=alevel,blevel=blevel,clevel=clevel,type="A",ordABC=ordABC)
        D3b <- D3function(alevel=alevel,blevel=blevel,clevel=clevel,type="B",ordABC=ordABC)
        D3c <- D3function(alevel=alevel,blevel=blevel,clevel=clevel,type="C",ordABC=ordABC)
        zero.1 <- matrix(0,ncol=(2*(r2+r3)),nrow=r1)
        zero.2l <- matrix(0,ncol=(2*r1),nrow=r2)
        zero.2r <- matrix(0,ncol=(2*r3),nrow=r2)
        zero.3 <- matrix(0,ncol=(2*(r1+r2)),nrow=r3)
        L3.1 <- cbind(D3a, -D3a, -diag(r1), diag(r1), zero.1)
        L3.2 <- cbind(D3b, -D3b, zero.2l, -diag(r2), diag(r2), zero.2r)
        L3.3 <- cbind(D3c, -D3c, zero.3, -diag(r3), diag(r3))
        L3 <- rbind(L3.1,L3.2,L3.3)
        return(L3)
    }



    ## fac.level <- c(4,7,4)
    ## ord <- c(TRUE,TRUE,TRUE)
    
    levelIndex <- CreatelevelIndex(fac.level=fac.level,ord.fac=ord.fac,Gorder=Gorder)
    use.ind <- (levelIndex$plus==1) * (levelIndex$dif==0)
    Index.use <- levelIndex[use.ind==1,]
    Fac.index <- Index.use[,regexpr("Fac",colnames(Index.use))>0]
    coef.length <- sum(levelIndex$length)


    L.list <- list()
    
    for(z in 1:nrow(Index.use)){
        if(Index.use$order[z]==1){
            level.main <- fac.level[which(Fac.index[z,]==1)]
            ord.main <- ord.fac[which(Fac.index[z,]==1)]
            L.list[[z]] <- L1function(nlevel=level.main,ord=ord.main)
        }else if(Index.use$order[z]==2){
            level.main <- fac.level[which(Fac.index[z,]==1)]
            ord.main <- ord.fac[which(Fac.index[z,]==1)]
            L.list[[z]] <- L2function(alevel=level.main[1],
                                      blevel=level.main[2],ordAB=ord.main)
        }else if(Index.use$order[z]==3){
            level.main <- fac.level[which(Fac.index[z,]==1)]
            ord.main <- ord.fac[which(Fac.index[z,]==1)]
            L.list[[z]] <- L3function(alevel=level.main[1],blevel=level.main[2],
                                      clevel=level.main[3],
                                      ordABC=ord.main)
        }
    }   
    
    Lm <- matrix(0,nrow=0,ncol=(coef.length))
    ## Need to make this generalizable
    for(l in 1:length(L.list)){
        L.minilist <- list()        
        for(ind in 1:length(L.list)){
            if(ind!=l){
                L.minilist[[ind]] <- matrix(0,nrow=nrow(L.list[[l]]),
                                            ncol=ncol(L.list[[ind]]))
            }else if(ind==l){
                L.minilist[[ind]] <- L.list[[l]]
            }            
        }
        Lb <- do.call(cbind,L.minilist)
        Lm <- rbind(Lm,Lb)
    }

    l.slack <- lengthSlack(fac.level=fac.level,ord.fac=ord.fac,Gorder=Gorder,
                           facCons=facCons)$l.slack
    Lm <- cbind(0,Lm)
    Lm <- rbind(0, Lm)
    Lslack <- matrix(0,nrow=nrow(Lm),ncol=l.slack)
    
    L <- cbind(Lm,Lslack)

    return(L)
}

