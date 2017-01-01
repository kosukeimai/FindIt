##############################
## stab.GashANOVA.R
## 2016/11/28
## Naoki Egami
############################## 
stab.CausalANOVA <- function(object,block.id=NULL,boot=500,seed=1234){

    if(class(object)[1] != "CausalANOVA" & class(object)[1] != "CausalANOVAFit"){
        warning("This function takes the output from CausalANOVA")
    }


    data <- object$data.orig
    pair.id <- object$pair.id

    if(is.null(block.id)==FALSE){
        data$block.id <- block.id    
    }
    
    formula <- object$formula
    Gorder <- object$nway
    diff <- object$diff
    fac.level <- object$fac.level
    ord.fac <- object$ord.fac
    cost <- object$cost
    eps <- object$eps
    
    if(diff==TRUE){
        data <- data[order(pair.id),]
        if(is.null(block.id)==FALSE) {block.id <- data$block.id}
        side <- rep(c(1,0),times=nrow(data)/2)
        data1 <- data[side==1,]
        data2 <- data[side==0,]
        X1 <- model.matrixBayes(formula,data=data1)
        X2 <- model.matrixBayes(formula,data=data2)
        X <- X1 - X2 
        y <- model.frame(formula,data=data1)[,1]
    }else if(diff==FALSE){
        X <- model.matrixBayes(formula,data=data)
        y <- model.frame(formula,data=data)[,1]
    }
    
    if(is.null(block.id)==FALSE){
        data$block.id <- NULL
    }
    
    n.stab <- nrow(X)
    index.stab <- seq(1:n.stab)
    
    coef.stab.v <- c()
    range.stab.v <- c()
    AME.stab.v <- c()
    AMIE.stab.v <- c()
    
    set.seed(seed)
    for(i in 1:boot){
        if(is.null(block.id)==FALSE){         
            ind.block <- sample(unique(block.id),size=length(unique(block.id)),replace=TRUE)
            data.block <- data
            if(diff==TRUE) data.block$side.block <- side
            data.block$block.id <- block.id
            data.block <- do.call("rbind", lapply(seq(1:length(ind.block)), 
                                                  FUN= function(x) subset(data.block, 
                                                      block.id==ind.block[x])))
            if(diff==TRUE){
                data.stab <- rbind(data.block[data.block$side.block==1, ],
                                   data.block[data.block$side.block==0, ])
                pair.stab <- rep(seq(1:nrow(data.block[data.block$side.block==1, ])),2)
            }else{
                data.stab <- data.block
                pair.stab <- NULL
            }

            data.stab$side.block <- NULL
            data.stab$block.id <- NULL         
            
        }else{
            ## No block
            ind.chosen <- sample(index.stab,size=length(index.stab),replace=TRUE)
            if(diff==TRUE){
                data.stab <- rbind(data1[ind.chosen, ],data2[ind.chosen, ])
                pair.stab <- rep(seq(1:nrow(data1[ind.chosen,])),2)
            }else{
                data.stab <- data[ind.chosen, ]
                pair.stab <- NULL
            }
        }
        gash.stab  <- CausalANOVAFit(formula,data=data.stab,pair.id=pair.stab,diff=diff,
                                     fac.level=fac.level,eps=eps,
                                     ord.fac=ord.fac,cost=cost, nway=Gorder,verbose=FALSE)
        
        coef.stab.v <- cbind(coef.stab.v,c(gash.stab$int, unlist(gash.stab$coef)))
        sum.gash <- rangeCausalANOVAFit(gash.stab,verbose=FALSE)
        range.stab.v <- cbind(range.stab.v,sum.gash$range)
        AME.stab.v <- cbind(AME.stab.v,sum.gash$AME)
        AMIE.stab.v <- cbind(AMIE.stab.v,sum.gash$AMIE)
        
        if(i %% 50 ==0){
            print(paste(round(i*(100/boot)),"% done.",sep=""))
        }
    }
    rownames(range.stab.v) <- rangeCausalANOVAFit(gash.stab,verbose=FALSE)$range.name

    StabBase <- rangeCausalANOVAFit(object,verbose=FALSE)
    ## Stability
    ## Range
    digit <- nchar(eps) - 1
    range.stab <- cbind(round(StabBase$range,digits=digit),
                        apply(apply(sign(round(range.stab.v,digits=digit)),
                                    2,function(x) x==1),1,mean))
    colnames(range.stab) <- c("Range","stability")

    ## ## AME/AMIE
    ## AMIE.main <- StabBase$AMIE.na
    ## AMIE.stab <- cbind(StabBase$AMIE,
    ##                    apply(apply(sign(AMIE.stab.v),2,function(x) x==sign(AMIE.main)),1,mean))
    ## colnames(AMIE.stab) <- c("AME/AMIE","stability")

    ## AME
    AME.main <- StabBase$AME
    AME.stab <- cbind(round(StabBase$AME,digits=digit),
                      apply(apply(sign(round(AME.stab.v,digits=digit)),2,
                                  function(x) x==sign(round(AME.main,digits=digit))),1,mean))
    stab.t <- AME.stab[round(AME.main,digits=digit)==0,2]
    AME.stab[round(AME.main,digits=digit)==0,2] <- 1 - stab.t
    colnames(AME.stab) <- c("AME",
                            "stability") 
    
    output <- list("range.stab"=range.stab,"AME.stab"=AME.stab,
                   "coef.stab.v"=coef.stab.v,
                   "range.stab.v"=range.stab.v,
                   "AME.stab.v"=AME.stab.v,"AMIE.stab.v"=AMIE.stab.v)    
}
    
