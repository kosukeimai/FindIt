##################################################
## Naoki Egami
## 2016/08/22
##################################################
lengthSlack <- function(fac.level,ord.fac,Gorder=2, facCons=FALSE){
    ## This does not depend on Gorder.
    if(facCons==FALSE){
        n.fac <- length(fac.level)
        length.slack <- rep(0,times=n.fac)
        
        for(i in 1:n.fac){
            if(ord.fac[i]==FALSE){
                length.slack[i] <- fac.level[i]*(fac.level[i]-1)/2
            }else{
                length.slack[i] <- (fac.level[i]-1)
            }
        }
        l.slack <- sum(length.slack)
        output <- list("l.slack"=l.slack,
                       "length.full"=length.slack)
    }else if(facCons==TRUE){
        n.fac <- length(fac.level)
        first <- n.fac
        second <- choose(n.fac,2)
        if(Gorder==3){
            third <- choose(n.fac,3)
            l.slack <- first + second + third
        }else if(Gorder==2){
            l.slack <- first + second 
        }
        output <- list("l.slack"=l.slack)
    }
    return(output)
}
