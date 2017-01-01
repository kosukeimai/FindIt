########################################
## cv.GashANOVA.R
## 2016/11/27
## Naoki Egami
########################################
cv.CausalANOVA <- function(formula,data,cv.cost=c(0.1,0.3,0.5,0.7,1.0),type="bin",
                           pair.id=NULL,nway=2,diff=TRUE,eps=1e-5,nfolds=10,seed=1234){

    if(diff==TRUE & is.null(pair.id)) warning("Set pair.id when diff=TRUE.")
    if(any(cv.cost<= 0)){
        warning("cv.cost should be greater than 0.")
    }
    if(any(cv.cost > 1)){
        warning("cv.cost should be smaller than or equal to 1. 1 corresponds to no regularization.")
    }
    
    if(type=="bin"){
        print("binary outcomes are used.")
    }else if (type=="cont"){
        print("continuous outcomes are used.")
    }

    Gorder <- nway
    
    ## We convert the data into the data only with necessary variables. 
    data.main <- model.frame(formula,data=data)
    data.x <- model.frame(formula,data=data)[,-1]
    
    all.fac <- all(unlist(lapply(data.x,FUN=function(x) is.element("factor",class(x)))))
    if(all.fac==FALSE) warning("Design matrix should contain only factors.")
    rm(all.fac)
    
    ## Extract information. 
    fac.level <- unlist(lapply(data.x,FUN=function(x) length(levels(x))))
    ord.fac <- unlist(lapply(data.x,FUN=function(x) is.ordered(x)))    
    for(i in 1:ncol(data.x)){
        if(ord.fac[i]==TRUE){
            data.main[,(i+1)] <- factor(data.main[,(i+1)],ordered=FALSE)
        }
    }
    ord.fac[fac.level==2] <- FALSE
    
    ## Print the status
    print.fac <- cbind(fac.level,as.data.frame(ord.fac))
    colnames(print.fac) <- c("levels","ordered")
    print(print.fac)
    rm(print.fac)
    rm(data.x)
    
    if(Gorder==2){
        formula <- as.formula(paste(all.vars(formula)[1]," ~ .*.", sep=""))
    }else if(Gorder==3){
        formula <- as.formula(paste(all.vars(formula)[1]," ~ .*.*.", sep=""))
    }
    data <- model.frame(formula,data=data.main)
    rm(data.main)
    
    if(diff==TRUE){
        data <- data[order(pair.id),]
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

    set.seed(seed)
    foldid <- sample(rep(seq(nfolds), length = nrow(X)))
    cv.result <- as.list(seq(nfolds))
    for (i in seq(nfolds)) {
        which = foldid == i
        
        if(diff==TRUE){
            data.cv <- rbind(data1[!which, ],data2[!which, ])
            pair.cv <- rep(seq(1:nrow(data1[!which,])),2)
        }else{
            data.cv <- data[!which, ]
            pair.cv <- NULL
        }
        coef.cv <- c()        
        for(z in 1:length(cv.cost)){
            gash.cv  <- CausalANOVAFit(formula=formula,data=data.cv,pair.id=pair.cv,diff=diff,
                                       fac.level=fac.level,eps=eps,
                                       ord.fac=ord.fac,cost=cv.cost[z],nway=Gorder,
                                       verbose=FALSE)
            coef.cv <- cbind(coef.cv,c(gash.cv$int, unlist(gash.cv$coef)))
        }
        
        ## Compute Missclassification 
        X.left <- X[which,]
        if (is.matrix(y)) y.left = y[which,] else y.left = y[which]
        y.pred <- cbind(1,X.left) %*% coef.cv
        y.left.mat <- matrix(rep(y.left,length(cv.cost)),ncol=length(cv.cost))
        if(type=="bin"){
            y.pred.bin <- y.pred > 0.5
            cv.result[[i]] <- y.left.mat - y.pred.bin
        }else if (type=="cont"){
            cv.result[[i]] <- y.left.mat - y.pred
        }        
        print(paste(round(i*(100/nfolds)),"% done.",sep=""))
    }

    cv.result.mat <- do.call(rbind,cv.result)
    if(type=="bin"){
        cv.error  <- apply(cv.result.mat,2,function(x) 1-mean(x==0))
        cv.each.list <- lapply(X=cv.result,FUN=function(x) apply(x,2,function(x) 1-mean(x==0)))
    }else if(type=="cont"){
        cv.error  <- apply(cv.result.mat,2,function(x) mean(x^2))
        cv.each.list <- lapply(X=cv.result,FUN=function(x) apply(x,2,function(x) mean(x^2)))
    }
    names(cv.error) <- cv.cost
    cv.min <- cv.cost[which.min(cv.error)]
    
    cv.each.mat <- do.call(rbind,cv.each.list)
    colnames(cv.each.mat) <- cv.cost

    cv.sd.each <- apply(cv.each.mat,2,sd)
    cv.sd1.value <- min(cv.error) + cv.sd.each[which.min(cv.error)]
    cv.sd1 <- min(cv.cost[cv.error < cv.sd1.value])
    
    output <- list("cv.min"=cv.min,"cv.sd1"=cv.sd1,
                   "cv.error"=cv.error,
                   "cv.each.mat"=cv.each.mat,
                   "cv.cost"=cv.cost)
    class(output) <- "cv.CausalANOVA"
    return(output)
}


