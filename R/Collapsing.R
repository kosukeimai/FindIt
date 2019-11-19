Collapsing <- function(fit){
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

    fac.level <- fit$fac.level
    ord.fac  <- fit$ord.fac
    AME <- fit$AME
    AMIE2 <- fit$AMIE2
    AMIE3 <- fit$AMIE3
    eps <- fit$eps   
    fac.name <- all.vars(fit$formula)[-1]
    n.fac <- length(fac.level)
    indTwo <- fit$indTwo
    indThree <- fit$indThree
    Gorder <- fit$Gorder
    
    levelIndex <- CreatelevelIndex(fac.level=fac.level,ord.fac=ord.fac,Gorder=Gorder,
                                   indTwo=indTwo, indThree=indThree)
    use.ind <-  (levelIndex$plus==1)*(levelIndex$dif==0)
    Index.use <- levelIndex[use.ind==1,]
    Fac.index <- levelIndex[,regexpr("Fac",colnames(levelIndex))>0]
    Fac.Ind.use <- Fac.index[use.ind==1,]

    order.f <- attr(terms(fit$formula, data=fit$data), "order")
    Fac.Ind.use1 <- Fac.Ind.use[order.f==1,]
    if(any(order.f==2)) Fac.Ind.use2 <- Fac.Ind.use[order.f==2,]
    if(any(order.f==3)) Fac.Ind.use3 <- Fac.Ind.use[order.f==3,]
    
       
    ## Do Collapsing for each factor
    Collapse <- list()
    for(z in 1:n.fac){
        ## First Order 
        ## type.ind <- z    
        MainD1 <- D1function(nlevel=fac.level[z],ord=ord.fac[z])
        MainDif <- t(MainD1 %*% AME[[z]])
        
        ## Two-ways
        onlyOne <- sum(Fac.Ind.use[,z])==1
        if(any(order.f==2)==TRUE) yesTwo <-  sum(Fac.Ind.use2[,z]==1)>0 else yesTwo <- FALSE
        if(any(order.f==3)==TRUE) yesThree <- sum(Fac.Ind.use3[,z]==1)>0 else yesThree <- FALSE
        if(onlyOne==TRUE){
            ## No Interaction
            Dif <- abs(MainDif)
            Collapse[[z]] <- apply(Dif <= eps, 2, all)
            ## print(Dif)
        }else if(yesTwo==TRUE){

            int2.index <- which(Fac.Ind.use2[,z]==1)
            D2.mat <- matrix(NA, ncol=ncol(MainDif), nrow=0)
            for(w in int2.index){
                ## Setup the coefficients
                AMIE2.u <- AMIE2[[w]]

                ## Construct D matrix
                Fac1 <- min(which(Fac.Ind.use2[w,]==1))
                Fac2 <- max(which(Fac.Ind.use2[w,]==1))
                if(z == Fac1) type.d2 <- "A"
                if(z == Fac2) type.d2 <- "B"
                D2  <- D2function(alevel=fac.level[Fac1],blevel=fac.level[Fac2],
                                  type=type.d2, ordAB=ord.fac[c(Fac1,Fac2)])
                D2.mat0 <- D2 %*% AMIE2.u
                D2.mat01  <- matrix(D2.mat0, ncol=ncol(MainDif))
                D2.mat <- rbind(D2.mat, D2.mat01)
            }
            D2.mat <- abs(D2.mat)
            if(yesThree==FALSE){                
                Dif <- abs(rbind(MainDif, D2.mat))
                Collapse[[z]] <- apply(Dif <= eps, 2, all)
            }else if(yesThree==TRUE){
                int3.index <- which(Fac.Ind.use3[,z]==1)
                D3.mat <- matrix(NA, ncol=ncol(MainDif), nrow=0)

                for(w in int3.index){
                    ## Setup the coefficients
                    AMIE3.u <- AMIE3[[w]]
                    
                    ## Construct D matrix
                    Fac1 <- which(Fac.Ind.use3[w,]==1)[1]
                    Fac2 <- which(Fac.Ind.use3[w,]==1)[2]
                    Fac3 <- which(Fac.Ind.use3[w,]==1)[3]
                    if(z == Fac1) type.d3 <- "A"
                    if(z == Fac2) type.d3 <- "B"
                    if(z == Fac3) type.d3 <- "C"
                    D3  <- D3function(alevel=fac.level[Fac1], blevel=fac.level[Fac2], clevel=fac.level[Fac3],
                                      type=type.d3, ordAB=ord.fac[c(Fac1,Fac2, Fac3)])
                    D3.mat0 <- D3 %*% c(AMIE3.u)
                    D3.mat01  <- matrix(D3.mat0, ncol=ncol(MainDif))
                    D3.mat <- rbind(D3.mat, D3.mat01)
                }
                D3.mat <- abs(D3.mat)         
                Dif <- abs(rbind(MainDif, D2.mat, D3.mat))
                Collapse[[z]] <- apply(Dif <= eps, 2, all)
            }
        }
    }
    
    ## I got Collapsing Index, then decide which levels will be collapsed.
    collapse.level <- list()
    for(z in 1:n.fac){
        adj <- matrix(0, ncol=fac.level[z], nrow=fac.level[z])
        if(ord.fac[z]==TRUE){
            for(i in 1:length(Collapse[[z]])){
                adj[i, (i+1)] <- as.numeric(Collapse[[z]][i])
            }
            adj <- adj + t(adj)
        }else if(ord.fac[z]==FALSE){
            ref <- combn(seq(1:fac.level[z]),2)
            for(i in 1:length(Collapse[[z]])){
                adj[ref[1,i], ref[2,i]] <- as.numeric(Collapse[[z]][i])
            }
            adj <- adj + t(adj)
        }        
        g <- graph_from_adjacency_matrix(adj, mode="undirected")
        collapse.level[[z]] <- components(g)$membership
    }
    names(collapse.level) <- fac.name
    
    ## Combine two weights. 
    output <- list("Collapse.Index"=Collapse, "collapse.level"=collapse.level)
    return(output)
}





