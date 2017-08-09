########################################
## cv.GashANOVA.R
## 2016/11/27
## Naoki Egami
########################################
cv.CausalANOVA <- function(formula, int2.formula=NULL, int3.formula=NULL,
                           data, nway=1, pair.id=NULL, diff=FALSE,
                           cv.collapse.cost=c(0.1,0.3,0.7), nfolds=5,
                           screen=TRUE, screen.type="fixed", screen.num.int=3,                                                
                           family="binomial", cluster=NULL,                   
                           maxIter=50, eps=1e-5, seed=1234,         
                           fac.level=NULL,ord.fac=NULL, verbose=TRUE){

    ## ############################ 
    ## HouseKeeping (the same as CausalANOVA function)
    ## ############################
    cv.cost <- cv.collapse.cost
    if(missing(fac.level)) fac.level <- NULL
    if(missing(ord.fac)) ord.fac <- NULL    
    if(any(cv.cost > 1) | any(cv.cost < 0) ){
        stop("Specify 'cost.cv' between 0 and 1")
    }
    if(diff==TRUE & is.null(pair.id)==TRUE){
        stop("When 'diff=TRUE', specify 'pair.id'.")
    }
    if((nway %in% c(1, 2,3))==FALSE){
        stop("'nway' should be 1, 2 or 3.")
    }
    if(screen==TRUE & nway==1){
        warning("When 'nway=1', no screening is needed ('screen=FALSE').")
        screen <- FALSE
    }
    if((family %in% c("gaussian","binomial"))==FALSE){
        stop("'family' should be 'gaussian' or 'binomial'.")
    }
    y <- model.frame(formula,data=data)[,1]
    if(family=="gaussian" & is.numeric(y)==FALSE){
        stop("When 'family=gaussian', outcomes should be 'numeric'.") 
    }
    if(family=="binomial" & all(y %in% c(0,1))==FALSE){
        stop("When 'family=binomial', outcomes should be 0 or 1.") 
    }
    rm(y)
    if((screen.type %in% c("fixed","cv.min","cv.1Std"))==FALSE){
        stop("'screen.type' should be 'fixed', 'cv.min' or 'cv.1Std'.")
    }    
    if(is.null(int3.formula)==FALSE & nway!=3){
        warning("When 'int3.formula' is specified, 'nway=3'.")
        nway <- 3
    }
    if(is.null(int2.formula)==TRUE & is.null(int3.formula)==FALSE){
        stop("Cannot specify 'int3.formula' without specifying 'int2.formula'.")
    }
    if(is.null(int2.formula)==FALSE & screen==TRUE){
        warning("When 'int2.formula' is specified, 'screen' should be FALSE.")
        screen <- FALSE
    }
    if(nway==3 & screen==TRUE & screen.type=="fixed" & screen.num.int < 3){
        stop("'screen.num.int >=3' when 'nway=3'.")
    }

    ## Factors only
    all.fac <- all(unlist(lapply(model.frame(formula,data=data)[,-1],
                                 FUN=function(x) is.element("factor",class(x)))))
    if(all.fac==FALSE) stop("Design matrix should contain only factors.")
    rm(all.fac);
    
    ## when int.formual != NULL
    if(is.null(int2.formula)==TRUE){
        internal.int <- TRUE
    }else if(is.null(int2.formula)==FALSE){
        internal.int <- FALSE
        if(is.null(int3.formula)==TRUE){
            formula <- as.formula(paste(as.character(formula)[2], "~",
                                        as.character(formula)[3], "+", as.character(int2.formula)[2], sep=""))
        }else if(is.null(int3.formula)==FALSE){
            formula <- as.formula(paste(as.character(formula)[2], "~",
                                        as.character(formula)[3],
                                        "+", as.character(int2.formula)[2],
                                        "+", as.character(int3.formula)[2],
                                        sep=""))
        }
    }
    ## ######################################
    ## #####################################
    ## ###################################### 

    Gorder <- nway
    formula.orig <- formula
    ## We convert the data into the data only with necessary variables. 
    data.main <- model.frame(formula,data=data)
    data.x <- model.frame(formula,data=data)[,-1]
    
    ## Extract information.
    if(is.null(fac.level)){
        fac.level <- unlist(lapply(data.x,FUN=function(x) length(levels(x))))
    }
    n.fac <- length(fac.level)
    
    if(is.null(ord.fac)){
        ord.fac <- unlist(lapply(data.x,FUN=function(x) is.ordered(x)))
    }
    for(i in 1:ncol(data.x)){
        if(ord.fac[i]==TRUE){
            data.main[,(i+1)] <- factor(data.main[,(i+1)],ordered=FALSE)
        }
    }

    if(any(ord.fac[fac.level==2]==TRUE)){
        if(verbose==TRUE) cat("Note: binary factors are recoded to unordered.")
        ord.fac[fac.level==2] <- FALSE
    }
    
    ## Print the status
    if(verbose==TRUE){
        print.fac <- cbind(fac.level,as.data.frame(ord.fac))
        colnames(print.fac) <- c("levels","ordered")
        cat("\nCheck: the number of levels for factors and whether they are ordered.\n")
        print(print.fac); rm(print.fac)
    }



    ## We need Formula (when user specified), internal.int
    rm(data.x);
    
    data <- model.frame(formula,data=data.main); rm(data.main)    
    data.orig <- data
    

    ## ############################
    ## (0.3) Take Differences for Paired Data
    ## ############################ 
    if(diff==TRUE){
        data <- data[order(pair.id),]
        side <- rep(c(1,0),times=nrow(data)/2)
    }else{
        side <- NULL
    }
    
    if(diff==TRUE){
        data1 <- data[side==1,]
        data2 <- data[side==0,]
        y <- model.frame(formula,data=data1)[,1]
    }else if(diff==FALSE){
        y <- model.frame(formula,data=data)[,1]
    }

    ## We always Regularize, so no need for cost==1

    ## ############################
    ## (2.0) CV Setting
    ## ############################ 
    set.seed(seed)
    foldid <- sample(rep(seq(nfolds), length = length(y)))
    cv.result <- as.list(seq(nfolds))
    for (i in seq(nfolds)) {
        which = foldid == i
        
        if(diff==TRUE){
            data.cv <- rbind(data1[!which, ],data2[!which, ])
            pair.cv <- rep(seq(1:nrow(data1[!which,])),2)
            data1.cv.u <- data1[which,]
            data2.cv.u <- data2[which,]
        }else{
            data.cv <- data[!which, ]
            data.cv.u <- data[which, ]
            pair.cv <- NULL
        }
        ## By screening, formula can be different!
        coef.cv <- list()
        X.l <- list()
        for(z in 1:length(cv.cost)){
            ## Special is data=data.cv,cost=cv.cost[z],pair.id=pair.cv,
            ## For different cost, I would get different data. 
            gash.cv  <- CausalANOVAFit(formula=formula, internal.int=internal.int,
                                       data=data.cv,cost=cv.cost[z],pair.id=pair.cv,
                                       nway=Gorder, diff=diff, eps=eps, family=family,
                                       verbose=FALSE, cluster=NULL,
                                       screen=screen, screen.type=screen.type,
                                       screen.num.int=screen.num.int,
                                       maxIter=maxIter)
            coef.cv[[z]] <- c(gash.cv$intercept, unlist(gash.cv$coefs))
            if(diff==TRUE){
                ## Model.matrixBayes expand the matrix, even if there is no variation!!
                ## This is why we do not need to see errors here!
                X1.l <- model.matrixBayes(gash.cv$formula, data=data1.cv.u)
                X2.l <- model.matrixBayes(gash.cv$formula, data=data2.cv.u)
                X.l[[z]] <- X1.l - X2.l 
            }else if(diff==FALSE){
                X.l[[z]] <- model.matrixBayes(gash.cv$formula,data=data.cv.u)
            }
        }
        
        ## Compute Missclassification        
        if (is.matrix(y)) y.left = y[which,] else y.left = y[which]
        y.pred <- matrix(NA, ncol=length(cv.cost), nrow=sum(which))
        for(z in 1:length(cv.cost)){
            y.pred[1:sum(which),z] <- cbind(1,X.l[[z]]) %*% coef.cv[[z]]
        }
        y.left.mat <- matrix(rep(y.left,length(cv.cost)),ncol=length(cv.cost)) ## True outcomes
        if(family=="binomial"){
            y.pred.bin <- y.pred > 0.5
            cv.result[[i]] <- y.left.mat - y.pred.bin
        }else if (family=="gaussian"){
            cv.result[[i]] <- y.left.mat - y.pred
        }        
        print(paste(round(i*(100/nfolds)),"% done.",sep=""))
    }

    cv.result.mat <- do.call(rbind,cv.result)
    if(family=="binomial"){
        cv.error  <- apply(cv.result.mat,2,function(x) 1-mean(x==0))
        cv.each.list <- lapply(X=cv.result,FUN=function(x) apply(x,2,function(x) 1-mean(x==0)))
    }else if(family=="gaussian"){
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
    
    output <- list("cv.min"=cv.min,"cv.1Std"=cv.sd1,
                   "cv.error"=cv.error,
                   "cv.each.mat"=cv.each.mat,
                   "cv.cost"=cv.cost)
    class(output) <- "cv.CausalANOVA"
    return(output)
}


