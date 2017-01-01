##################################################
## Naoki Egami
## 2016/08/22
##################################################

CreateWeights <- function(formula, Gorder=3,dif=TRUE,facCons=FALSE,data,side,type.ind,fac.level,ord.fac){

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

    if(dif==TRUE){
        data1 <- data[side==1,]
        data2 <- data[side==0,]
        model.frame  <- model.frame(formula,data=data1)
        contr <- rep(list("contr.sum"), ncol(model.frame) - 1)
        names(contr) <- colnames(model.frame)[-1]
        y <- model.frame[,1]
        x1 <- model.matrix(formula, data=data1,contrast=contr)[,-1]
        x2 <- model.matrix(formula, data=data2,contrast=contr)[,-1]
        X <- x1-x2
    }else if (dif==FALSE){
        model.frame  <- model.frame(formula,data=data)
        contr <- rep(list("contr.sum"), ncol(model.frame) - 1)
        names(contr) <- colnames(model.frame)[-1]
        y <- model.frame[,1]
        X <- model.matrix(formula, data=data,contrast=contr)[,-1]
    }
    lm.fit <- lm(y ~ X)
    coef.use <- coef(lm.fit)[-1]
    
    levelIndex <- CreatelevelIndex(fac.level=(fac.level-1),ord.fac=ord.fac, Gorder=Gorder)
    use.ind <-  (levelIndex$plus==1)*(levelIndex$dif==0)
    Index.use <- levelIndex[use.ind==1,]
    Fac.index <- levelIndex[,regexpr("Fac",colnames(levelIndex))>0]
    Fac.Ind.use <- Fac.index[use.ind==1,]
    Fac.level.useC <- (rep(1,nrow(Fac.Ind.use)) %x% t((fac.level-1))) * Fac.Ind.use
    Fac.level.use <- (rep(1,nrow(Fac.Ind.use)) %x% t((fac.level))) * Fac.Ind.use
    Index.use$start <- c(1,cumsum(Index.use$length)[-nrow(Index.use)]+1)
    Index.use$end <-  cumsum(Index.use$length)

    if(Gorder==2){
        Index.use$seq.order <- c(seq(1:sum(Index.use$order==1)),
                                 seq(1:sum(Index.use$order==2)))
    }else if(Gorder>2){
        Index.use$seq.order <- c(seq(1:sum(Index.use$order==1)),
                                 seq(1:sum(Index.use$order==2)),
                                 seq(1:sum(Index.use$order==3)))
    }
    
    Coef.list <- CreateCoef(coef.use=coef.use, 
                            Index.use=Index.use,
                            Fac.level.use=Fac.level.useC)

    if(facCons==TRUE){
        n.fac <- length(fac.level)
        MainDif <- list()
        IntDif <- list()
        ThreeDif <- list()
        for(z in 1:n.fac){
            Maincoef <- Coef.list$One.coef[[z]]
            MainD1 <- D1function(nlevel=fac.level[z],ord=ord.fac[z])
            MainDif[[z]] <- t(MainD1 %*% Maincoef)
        }
        
        ind.use.Two <- (Index.use$order==2)
        Index.useTwo <- Index.use[ind.use.Two==1,]
        Fac.useTwo <- Fac.Ind.use[ind.use.Two==1,]
        Fac.level.useTwo <- Fac.level.use[ind.use.Two==1,]
        ## IntDif <- matrix(0,nrow=0,ncol=qp)
        for(z in 1:nrow(Index.useTwo)){
            level.use <- as.numeric(Fac.level.useTwo[z,Fac.useTwo[z, ]==1])
            ord.use <- ord.fac[Fac.useTwo[z, ]==1]
            coef.int <- c(Coef.list$Two.coef[[z]])
            D2A  <- D2function(alevel=level.use[1],blevel=level.use[2],
                               type="A",ordAB=ord.use)
            D2B  <- D2function(alevel=level.use[1],blevel=level.use[2],
                               type="B",ordAB=ord.use)
            IntDif.A <- D2A %*% coef.int
            IntDif.B <- D2B %*% coef.int
            IntDif[[z]] <- c(IntDif.A, IntDif.B)

            if(Gorder>2){
                ## Three-ways   
                ind.use.Three <- (Index.use$order==3)
                Index.useThree <- Index.use[ind.use.Three==1,]
                Fac.useThree <- Fac.Ind.use[ind.use.Three==1,]
                Fac.level.useThree <- Fac.level.use[ind.use.Three==1,]
                ThreeDif <- list()
                
                for(z in 1:nrow(Index.useThree)){
                    level.use <- as.numeric(Fac.level.useThree[z,Fac.useThree[z, ]==1])
                    ord.use <- ord.fac[Fac.useThree[z, ]==1]                    
                    D3A  <- D3function(alevel=level.use[1],blevel=level.use[2],clevel=level.use[3],
                                       type="A",ordABC=ord.use)
                    D3B  <- D3function(alevel=level.use[1],blevel=level.use[2],clevel=level.use[3],
                                       type="B",ordABC=ord.use)
                    D3C  <- D3function(alevel=level.use[1],blevel=level.use[2],clevel=level.use[3],
                                       type="C",ordABC=ord.use)
                    coef.three <- c(Coef.list$Three.coef[[z]])
                    ThreeDif.A <- D3A %*% coef.three
                    ThreeDif.B <- D3B %*% coef.three
                    ThreeDif.C <- D3C %*% coef.three
                    ThreeDif[[z]] <- c(ThreeDif.A, ThreeDif.B, ThreeDif.C)
                }
            }            
        }
        
        Main.max <- unlist(lapply(1:length(MainDif),function(x) 1/max(abs(MainDif[[x]]))))
        Int.max <- unlist(lapply(1:length(IntDif),function(x) 1/max(abs(IntDif[[x]]))))        
        if(Gorder > 2){
            Three.max <- unlist(lapply(1:length(ThreeDif),function(x) 1/max(abs(ThreeDif[[x]]))))
            weight <- c(Main.max,Int.max, Three.max)
        }else if(Gorder==2){
            weight <- c(Main.max,Int.max)
        }
    }
    else if(facCons==FALSE){
        ## Main Effects 
        Maincoef <- Coef.list$One.coef[[type.ind]]
        MainD1 <- D1function(nlevel=fac.level[type.ind],ord=ord.fac[type.ind])
        MainDif <- t(MainD1 %*% Maincoef)
        qp <- ncol(MainDif)
        
        ## Two-ways
        ind.use.Two <- (Index.use$order==2)*(Fac.Ind.use[type.ind]==1)
        Index.useTwo <- Index.use[ind.use.Two==1,]
        Fac.useTwo <- Fac.Ind.use[ind.use.Two==1,]
        Fac.level.useTwo <- Fac.level.use[ind.use.Two==1,]
        IntDif <- matrix(0,nrow=0,ncol=qp)
        for(z in 1:nrow(Index.useTwo)){
            level.use <- as.numeric(Fac.level.useTwo[z,Fac.useTwo[z, ]==1])
            ord.use <- ord.fac[Fac.useTwo[z, ]==1]
            if(min(which(Fac.useTwo[z,]==1))==type.ind){
                type.d2 <- "A"}else{ type.d2 <- "B" }
            D2  <- D2function(alevel=level.use[1],blevel=level.use[2],
                              type=type.d2,ordAB=ord.use)
            coef.int <- c(Coef.list$Two.coef[[Index.useTwo[z,"seq.order"]]])
            IntDif.t <- D2 %*% coef.int
            IntDif.mat <- matrix(IntDif.t,ncol=qp)
            IntDif <- rbind(IntDif, IntDif.mat)
        }

        if(Gorder>2){
            ## Three-ways   
            ind.use.Three <- (Index.use$order==3)*(Fac.Ind.use[type.ind]==1)
            Index.useThree <- Index.use[ind.use.Three==1,]
            Fac.useThree <- Fac.Ind.use[ind.use.Three==1,]
            Fac.level.useThree <- Fac.level.use[ind.use.Three==1,]
            ThreeDif <- matrix(0,nrow=0,ncol=qp)
            
            for(z in 1:nrow(Index.useThree)){
                level.use <- as.numeric(Fac.level.useThree[z,Fac.useThree[z, ]==1])
                ord.use <- ord.fac[Fac.useThree[z, ]==1]
                if(min(which(Fac.useThree[z,]==1))==type.ind){
                    type.d3 <- "A"
                }else if (max(which(Fac.useThree[z,]==1))==type.ind){
                    type.d3 <- "C"
                }else{ type.d3 <- "B"}
                D3  <- D3function(alevel=level.use[1],blevel=level.use[2],clevel=level.use[3],
                                  type=type.d3,ordABC=ord.use)
                coef.three <- c(Coef.list$Three.coef[[Index.useThree[z,"seq.order"]]])
                ThreeDif.t <- D3 %*% coef.three
                ThreeDif.mat <- matrix(ThreeDif.t,ncol=qp)
                ThreeDif <- rbind(ThreeDif, ThreeDif.mat)
            }
        }

        if(Gorder == 2){
            Dif <- abs(rbind(MainDif,IntDif))
        }else if(Gorder > 2){
            Dif <- abs(rbind(MainDif,IntDif,ThreeDif))
        }
        weight <- apply(Dif,2,function(x) 1/max(x))
    }

    levelIndex <- CreatelevelIndex(fac.level=(fac.level),ord.fac=ord.fac, Gorder=Gorder)
    use.ind <-  (levelIndex$plus==1)*(levelIndex$dif==1) * (levelIndex$order==1)
    Index.use <- levelIndex[use.ind==1,]
    temp.w <- 1/(sqrt(fac.level[type.ind])*(fac.level[type.ind]+1))
    weight.u <- rep(temp.w,times=Index.use$length[type.ind])

    ## Combine two weights. 
    output <- list("weight.ols"=weight, "weight.fac"=weight.u)
    return(output)
}





