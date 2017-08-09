##################################################
## Naoki Egami
## 2016/08/22
## Updated on 2017/07/20
##################################################


## Global function: (it does change depending on the formula)
CreatelevelIndex <- function(fac.level,ord.fac, indTwo=NULL, indThree=NULL, Gorder){
    if(missing(Gorder)){
        Gorder <- 3 
    }
    n.fac <- length(fac.level)
    
    dif.level <- rep(0,times=n.fac)
    for(i in 1:n.fac){
        if(ord.fac[i]==TRUE){
            dif.level[i] <- fac.level[i] - 1
        }else{
            dif.level[i] <- fac.level[i]*(fac.level[i] - 1)/2
        }
    }
    
    fac.name <- paste("Fac",seq(1:n.fac),sep=".")

    ## Setup
    if(is.null(indTwo)==TRUE){
        indTwo <- combn(seq(1:n.fac),2)
    }
    if(is.null(indThree)==TRUE & Gorder==3){
        indThree <- combn(seq(1:n.fac),3)
    }

    ## %%%% IndTwo can be imported. 

    ## For one-way
    Index.one <- matrix(0,nrow=0,ncol=(5 + n.fac))
    BaseOne <- as.data.frame(cbind(rep(c(1,0),2),rep(c(0,1),each=2)))
    BaseOne <- cbind(rep(1,4),BaseOne)
    for(one in 1:n.fac){
        FacOne  <- matrix(0,ncol=n.fac,nrow=4)
        FacOne[,one] <- rep(c(1,2),each=2)
        length <- c(rep(fac.level[one],each=2),rep(dif.level[one],each=2))
        Index.one <- rbind(Index.one, cbind(length,BaseOne,FacOne))
    }
    colnames(Index.one) <- c("length","order","plus","dif",fac.name)
    ## Index.one

    ## For two-way
    if(Gorder>=2){
        Index.two <- matrix(0,nrow=0,ncol=(5 + n.fac))
        BaseTwo <- as.data.frame(cbind(rep(c(1,0),3),c(rep(0,2),rep(1,4))))
        BaseTwo <- cbind(rep(2,6),BaseTwo)
        for(z in 1:ncol(indTwo)){
            FacTwo  <- matrix(0,ncol=n.fac,nrow=6)
            FacTwo[,indTwo[,z]] <- cbind(rep(c(1,2,1),each=2),rep(c(1,1,2),each=2))
            fac.int.level <- fac.level[indTwo[1,z]]*fac.level[indTwo[2,z]]
            fac.int1 <- dif.level[indTwo[1,z]]*fac.level[indTwo[2,z]]
            fac.int2 <- fac.level[indTwo[1,z]]*dif.level[indTwo[2,z]]
            length <- c(rep(fac.int.level,each=2),rep(fac.int1,each=2),rep(fac.int2,each=2))
            Index.two <- rbind(Index.two, cbind(length,BaseTwo,FacTwo))
        }
        colnames(Index.two) <- c("length","order","plus","dif",fac.name)
    }
    ## Index.two
    
    if(Gorder == 3){
        ## For three-way
        Index.three <- matrix(0,nrow=0,ncol=(5 + n.fac))
        BaseThree <- as.data.frame(cbind(rep(c(1,0),4),c(rep(0,2),rep(1,6))))
        BaseThree <- cbind(rep(3,8),BaseThree)
        for(z in 1:ncol(indThree)){
            FacThree  <- matrix(0,ncol=n.fac,nrow=8)
            FacThree[,indThree[,z]] <- cbind(rep(c(1,2,1,1),each=2),
                                             rep(c(1,1,2,1),each=2),
                                             rep(c(1,1,1,2),each=2))
            fac.int.level <- fac.level[indThree[1,z]]*fac.level[indThree[2,z]]*fac.level[indThree[3,z]]
            fac.int1 <- dif.level[indThree[1,z]]*fac.level[indThree[2,z]]*fac.level[indThree[3,z]]
            fac.int2 <- fac.level[indThree[1,z]]*dif.level[indThree[2,z]]*fac.level[indThree[3,z]]
            fac.int3 <- fac.level[indThree[1,z]]*fac.level[indThree[2,z]]*dif.level[indThree[3,z]]
            length <- c(rep(fac.int.level,each=2),rep(fac.int1,each=2),
                        rep(fac.int2,each=2), rep(fac.int3,each=2))
            Index.three <- rbind(Index.three, cbind(length,BaseThree,FacThree))
        }
        colnames(Index.three) <- c("length","order","plus","dif",fac.name)
    }
    
    ## Add Start Index
    if(Gorder==1){
        Index.Mat <- Index.one
    }else if(Gorder ==2){
        Index.Mat <- rbind(Index.one, Index.two)
    }else if(Gorder ==3){
        Index.Mat <- rbind(Index.one, Index.two,Index.three)
    }
    start <- c(1,cumsum(Index.Mat$length)[-nrow(Index.Mat)]+1)
    
    Index.Mat <- cbind(start,Index.Mat)
    
    return(Index.Mat)
}


