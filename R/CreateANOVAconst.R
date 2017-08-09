##################################################
## Naoki Egami
## 2016/08/25
##################################################
CreateANOVAconst <- function(fac.level,ord.fac,Gorder,facCons, indTwo=NULL, indThree=NULL){
    ## Two-way Interactions
    TwoWayConst <- function(alevel,blevel,type="alpha"){
        ## type determines which to fix
        if(type=="beta"){
            ## Fix j (beta for AB)
            int.constAB <- diag(blevel) %x% t(rep(1,times=alevel))
        }else if(type=="alpha"){
            ## Fix i (alpha for AB)
            int.constAB <- do.call(cbind, replicate(blevel, diag(alevel), simplify=FALSE))
        }
        return(int.constAB)
    }

    ThreeWayConst <- function(alevel,blevel,clevel,type="alpha"){
        ## type determines which to average over
        if(type=="gamma"){
            ## Fix j (beta for AB)
            int.constABC <- do.call(cbind, replicate(clevel, diag(alevel*blevel), 
                                                     simplify=FALSE))
        }else if(type=="beta"){
            int.constBase <- do.call(cbind, replicate(blevel, diag(alevel), 
                                                      simplify=FALSE))
            ## Fix i (alpha for AB)
            int.constABC <- diag(clevel) %x% int.constBase
        }else if(type=="alpha"){
            int.constABC <- diag(blevel*clevel) %x% t(rep(1,times=(alevel)))
        }
        return(int.constABC)
    }


    ## I introduce Slack Matrix later
    ## I will add intercept later

    ## Even with ordered Categorical varialbe, almost the same

    ## We can fix this. 
    levelIndex <- CreatelevelIndex(fac.level=fac.level,ord.fac=ord.fac,Gorder=Gorder,
                                   indTwo=indTwo, indThree=indThree)
    Fac.index <- levelIndex[,regexpr("Fac",colnames(levelIndex))>0]
    use.ind <-  (levelIndex$plus==1)*(levelIndex$dif==0)
    Index.use <- levelIndex[use.ind==1,]
    Fac.Ind.use <- Fac.index[use.ind==1,]
    Index.ind <- seq(1:nrow(levelIndex))[use.ind==1]
    coef.length <- sum(levelIndex$length)
    l.slack <- lengthSlack(fac.level=fac.level,ord.fac=ord.fac,Gorder=Gorder,
                           facCons=facCons)$l.slack       
    
    ## The base matrix for each set of coefficients
    Base.C.list <- list()
    for(z in 1:nrow(Index.use)){
        if(Index.use$order[z]==1){
            mainlevel <- fac.level[z]
            Base.C.list[[z]] <- rep(1,mainlevel)
        }else if(Index.use$order[z]==2){
            level.temp <- (Fac.Ind.use[z,]==1) * fac.level
            mainlevel <- level.temp[level.temp!=0]
            ## for the same interaction coefficients, so we can take rbind here. 
            Base.C.list[[z]] <- rbind(TwoWayConst(alevel=mainlevel[1],blevel=mainlevel[2],type="alpha"),
                                      TwoWayConst(alevel=mainlevel[1],blevel=mainlevel[2],type="beta"))
            
        }else if(Index.use$order[z]==3){
            level.temp <- (Fac.Ind.use[z,]==1) * fac.level
            mainlevel <- level.temp[level.temp!=0]           
            Base.C.list[[z]] <- rbind(ThreeWayConst(alevel=mainlevel[1],blevel=mainlevel[2],
                                                    clevel=mainlevel[3],type="alpha"),
                                      ThreeWayConst(alevel=mainlevel[1],blevel=mainlevel[2],
                                                    clevel=mainlevel[3],type="beta"),
                                      ThreeWayConst(alevel=mainlevel[1],blevel=mainlevel[2],
                                                    clevel=mainlevel[3],type="gamma"))
        }
    }

    ## Combine into the larger Matrix using the index 
    Const <- matrix(0,nrow=0,ncol=coef.length)
    for(z in 1:nrow(Index.use)){
        if(Index.use$order[z]==1){
            right.n.col1 <- sum(levelIndex$length[(Index.ind[z]+2):nrow(levelIndex)])
            if(Index.ind[z]>1){
                left.n.col1 <- sum(levelIndex$length[1:(Index.ind[z]-1)])
                Const.temp <- c(rep(0,left.n.col1),Base.C.list[[z]],-Base.C.list[[z]], rep(0,right.n.col1))
            }else{
                Const.temp <- c(Base.C.list[[z]],-Base.C.list[[z]], rep(0,right.n.col1))
            }
            Const <- rbind(Const, Const.temp)
        }else if(Index.use$order[z]>1){
            left.n.col <- sum(levelIndex$length[1:(Index.ind[z]-1)])
            right.n.col <- sum(levelIndex$length[(Index.ind[z]+2):nrow(levelIndex)])
            Cb.left <- matrix(0,ncol=left.n.col,nrow=nrow(Base.C.list[[z]]))
            Cb.right <- matrix(0,ncol=right.n.col,nrow=nrow(Base.C.list[[z]]))
            Const.temp <- cbind(Cb.left, Base.C.list[[z]], -Base.C.list[[z]], Cb.right)
            Const <- rbind(Const, Const.temp)
        }
    }
    
    ## Add intercept and slack variables
    ConstFinal <- cbind(0, Const, matrix(0,nrow=nrow(Const),ncol=l.slack))
    return(ConstFinal)
}
