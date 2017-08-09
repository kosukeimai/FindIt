########################################
## Naoki Egami
## 2016/09/03
########################################

CreateCoef <- function(coef.use, Index.use, Fac.level.use, Gorder){

    OneANOVA <- function(Ovec){
        OvecF<- c(Ovec,-sum(Ovec))
        return(OvecF)
    }
    TwoANOVA <- function(Tmatrix){
        TmatrixF <- cbind(Tmatrix,-apply(Tmatrix,1,sum))
        TmatrixF2 <- rbind(TmatrixF,-apply(TmatrixF,2,sum))
        return(TmatrixF2)
    }
    ThreeANOVA <- function(Tarray,array.dim){
        TarrayF <- array(NA,(array.dim+1))
        ## Fill everything
        for(k in 1:(array.dim[3])){
            if(is.numeric(Tarray[,,k])){
                Tarray.k <- matrix(Tarray[,,k],nrow=array.dim[1],ncol=array.dim[2])
            }else{
                Tarray.k <- Tarray[,,k]
            }
            TarrayF[,,k] <- TwoANOVA(Tarray.k)
        }
        ## Fill the last array
        coef.mat.sum <- matrix(0,nrow=nrow(TarrayF[,,1]),ncol=ncol(TarrayF[,,1]))
        for(k in 1:(array.dim[3])){
            coef.mat.sum <- coef.mat.sum + TarrayF[,,k]
        }        
        TarrayF[,,(array.dim[3]+1)] <- - coef.mat.sum
        return(TarrayF)
    }

    
    One.coef <- list()
    Two.coef <- list()
    Three.coef <- list()
    one <- two <- three <- 0
    for(z in 1:nrow(Index.use)){
        if(Index.use$order[z]==1){
            one <- one + 1
            One.coef[[one]] <- OneANOVA(coef.use[Index.use$start[z]:Index.use$end[z]])
        }
        if(Index.use$order[z]==2){
            two <- two + 1
            two.level <- Fac.level.use[z,][Fac.level.use[z,]!=0]
            two.mat <- matrix(coef.use[Index.use$start[z]:Index.use$end[z]],
                              nrow=two.level[1],ncol=two.level[2])
            two.mat <- matrix(two.mat,nrow=two.level[1],ncol=two.level[2])
            Two.coef[[two]] <- TwoANOVA(Tmatrix=two.mat)
        }
        if(Index.use$order[z]==3){
            three <- three + 1
            three.level <- Fac.level.use[z,][Fac.level.use[z,]!=0]
            three.array <- array(coef.use[Index.use$start[z]:Index.use$end[z]],
                                 dim=three.level)
            Three.coef[[three]] <- ThreeANOVA(three.array,three.level)
        }
    }
    if(Gorder ==1){
        full <- One.coef
    }else if(Gorder==2){
        full <- append(One.coef, Two.coef)
    }else if(Gorder==3){
        full <- append(One.coef, Two.coef)
        full <- append(full, Three.coef)
    }
    output <- list("One.coef"=One.coef,
                   "Two.coef"=Two.coef,
                   "Three.coef"=Three.coef,
                   "full"=full)
    return(output)

}


