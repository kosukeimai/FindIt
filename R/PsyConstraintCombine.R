##################################################
## Naoki Egami
## 2016/08/22
##################################################
PsyConstraintCombine <- function(fac.level,ord.fac, Gorder){

    PsyIJfunction <- function(i,j,type.ind,plus,fac.level,ord.fac, Gorder){
        plus.t <- ifelse(plus==TRUE,1,0)
        
        ## Actually, level.vec will do the most work to keep track of names
        levelIndex <- CreatelevelIndex(fac.level=fac.level,ord.fac=ord.fac, Gorder=Gorder)
        Fac.index <- levelIndex[,regexpr("Fac",colnames(levelIndex))>0]
        use.ind <-  (levelIndex$plus==plus.t)*(levelIndex$dif==1)*(Fac.index[,type.ind]==2)
        Index.use <- levelIndex[use.ind==1,]
        Fac.Ind.use <- Fac.index[use.ind==1,]
        Index.ind <- seq(1:nrow(levelIndex))[use.ind==1]

        mainlevel <- fac.level[type.ind]
        ordMain <- ord.fac[type.ind]
        
        ## Create the Basline Vectors
        if(ordMain==FALSE){
            level.comb <- combn(seq(from=1,to=mainlevel),2)
            ind.c <- which(((level.comb[1,]==i) * (level.comb[2,]==j))==1)
            base <- rep(0,times=((mainlevel)*(mainlevel-1)/2))
            base[ind.c] <- 1
        }else if(ordMain==TRUE){
            base <- rep(0,times=(mainlevel-1))
            base[i] <- 1
        }

        ## Create Sublevel Matrix
        Fac.sub <-t(t(Fac.Ind.use==1) * fac.level)
        sublevel <- rep(0,times=nrow(Fac.Ind.use))

        for(z in 2:nrow(Index.use)){
            sublevel[z] <- prod(Fac.sub[z,][Fac.sub[z,]!=0])
        }
        
        ## The base matrix for each set of coefficients
        Base.P.list <- list()
        for(z in 1:nrow(Index.use)){
            if(z>1){
                Base.P.list[[z]] <- t(base) %x% diag(sublevel[z])
            }else{
                Base.P.list[[z]] <- base
            }
        }

        ## First Row 
        left.n.col1 <- sum(levelIndex$length[1:(Index.ind[1]-1)])
        right.n.col1 <- sum(levelIndex$length[(Index.ind[1]+1):nrow(levelIndex)])
        PsyIJ <- c(rep(0,left.n.col1),base,rep(0,right.n.col1))

        ## Need to make this generalizable
        ## I need to refer to the LevelIndex 
        for(p in 2:length(Base.P.list)){
            left.n.col <- sum(levelIndex$length[1:(Index.ind[p]-1)])
            Pb.left <- matrix(0,ncol=left.n.col,nrow=nrow(Base.P.list[[p]]))
            if((nrow(levelIndex)-Index.ind[p])>=2){
                right.n.col <- sum(levelIndex$length[(Index.ind[p]+1):nrow(levelIndex)])         
                Pb.right <- matrix(0,ncol=right.n.col,nrow=nrow(Base.P.list[[p]]))
                Pb <- cbind(Pb.left, Base.P.list[[p]],Pb.right)
            }else if(nrow(levelIndex)==(1+Index.ind[p])){
                Pb <- cbind(Pb.left, Base.P.list[[p]],
                            matrix(0,ncol=ncol(Base.P.list[[p]]),nrow=nrow(Base.P.list[[p]])))
            }else if(nrow(levelIndex)==(Index.ind[p])){
                Pb <- cbind(Pb.left, Base.P.list[[p]])
            }
            PsyIJ <- rbind(PsyIJ,Pb)
        }
        
        ## Add a column for Mu.
        PsyIJ <- cbind(0, PsyIJ)
        return(PsyIJ)
    }

    ## a <- PsyConstraint(type.ind=3,fac.level,ord.fac)



    ## Turn PhyMatrix Intro Constraint Form
    PsyConstraintFunction <- function(type.ind,fac.level,ord.fac,Gorder){
        ## Main Job is to combine PsyIJ and Matrix for slack variables
        
        levelIndex <- CreatelevelIndex(fac.level,ord.fac, Gorder)
        coef.length <- sum(levelIndex$length)
        coef.length.mu <- coef.length + 1
        
        l.slack <- lengthSlack(fac.level=fac.level,ord.fac=ord.fac)$l.slack
        slack.full <- lengthSlack(fac.level=fac.level,ord.fac=ord.fac)$length.full

        mainlevel <- fac.level[type.ind]
        ordMain <- ord.fac[type.ind]
        if(type.ind>1){
            jump <- cumsum(slack.full)[(type.ind-1)]
        }else{
            jump <- 0
        }
        
        PsyFinal <- matrix(0,nrow=0,ncol=(coef.length.mu+l.slack))

        if(ordMain==FALSE){
            for(i in 1:(mainlevel-1)){    
                for(j in (i+1):mainlevel){
                    PijPlus <- PsyIJfunction(i=i,j=j,type.ind=type.ind,plus=TRUE,
                                             fac.level=fac.level,ord.fac=ord.fac, Gorder=Gorder)
                    PijMinus <- PsyIJfunction(i=i,j=j,type.ind=type.ind,plus=FALSE,
                                              fac.level=fac.level,ord.fac=ord.fac,
                                              Gorder=Gorder)                
                    Pij <- PijPlus + PijMinus
                    
                    ## Create Matrix for Slack variables
                    base <- rep(0,times=l.slack)
                    level.comb <- combn(seq(from=1,to=mainlevel),2)
                    ind <- which(((level.comb[1,]==i) * (level.comb[2,]==j))==1)
                    ind.use <- ind + jump
                    base[ind.use] <- 1
                    PsySlackIJ <- t(base %x% t(rep(1,times=nrow(Pij))))
                    PsyFinalIJ <- cbind(Pij, -PsySlackIJ)
                    PsyFinal <- rbind(PsyFinal,PsyFinalIJ)
                }
            }
        }else if(ordMain==TRUE){
            for(i in 1:(mainlevel-1)){    
                PijPlus <- PsyIJfunction(i=i,j=(i+1),type.ind=type.ind,plus=TRUE,
                                         fac.level=fac.level,ord.fac=ord.fac,
                                         Gorder=Gorder)
                PijMinus <- PsyIJfunction(i=i,j=(i+1),type.ind=type.ind,plus=FALSE,
                                          fac.level=fac.level,ord.fac=ord.fac,
                                          Gorder=Gorder)
                Pij <- PijPlus + PijMinus
                
                ## Create Matrix for Slack variables
                base <- rep(0,times=l.slack)
                ind <- i
                ind.use <- ind + jump
                base[ind.use] <- 1
                PsySlackIJ <- t(base %x% t(rep(1,times=nrow(Pij))))                      
                PsyFinalIJ <- cbind(Pij, -PsySlackIJ)
                PsyFinal <- rbind(PsyFinal,PsyFinalIJ)
            }        
        }    
        return(PsyFinal)
    }
    
    
    n.fac <- length(fac.level)
    levelIndex <- CreatelevelIndex(fac.level,ord.fac,Gorder)
    coef.length <- sum(levelIndex$length)
    coef.length.mu <- coef.length + 1    
    l.slack <- lengthSlack(fac.level=fac.level,ord.fac=ord.fac)$l.slack
    PsyConstraintF <- matrix(0, ncol=(coef.length.mu+l.slack),nrow=0)
    for(psy in 1:n.fac){
        Psy.temp <- PsyConstraintFunction(type.ind=psy,fac.level=fac.level,ord.fac=ord.fac,
                                          Gorder=Gorder)
        PsyConstraintF <- rbind(PsyConstraintF, Psy.temp)
    }
    return(PsyConstraintF)
}






