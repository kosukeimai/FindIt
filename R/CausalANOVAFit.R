CausalANOVAFit <- function(formula,data,cost,pair.id,
                           nway=2,diff=TRUE,eps=1e-5,
                           fac.level=NULL,ord.fac=NULL,
                           verbose=TRUE){
    
    if(diff==TRUE & is.null(pair.id)) warning("Set pair.id when diff=TRUE.")
    if(cost <= 0) warning("cost should be greater than 0.")
    if(cost > 1) warning("cost should be smaller than or equal to 1. 1 corresponds to no regularization.")   
    
    Gorder <- nway
    ##
    formula.orig <- formula
    ## We convert the data into the data only with necessary variables. 
    data.main <- model.frame(formula,data=data)
    data.x <- model.frame(formula,data=data)[,-1]

    all.fac <- all(unlist(lapply(data.x,FUN=function(x) is.element("factor",class(x)))))
    if(all.fac==FALSE) warning("Design matrix should contain only factors.")
    rm(all.fac)
    
    ## Extract information.
    if(is.null(fac.level)){
        fac.level <- unlist(lapply(data.x,FUN=function(x) length(levels(x))))
    }
    if(is.null(ord.fac)){
        ord.fac <- unlist(lapply(data.x,FUN=function(x) is.ordered(x)))
    }
    for(i in 1:ncol(data.x)){
        if(ord.fac[i]==TRUE){
            data.main[,(i+1)] <- factor(data.main[,(i+1)],ordered=FALSE)
        }
    }

    if(any(ord.fac[fac.level==2]==TRUE)){
        print("Note that binary factors should always be unordered.")
        ord.fac[fac.level==2] <- FALSE
    }
    
    ## Print the status
    print.fac <- cbind(fac.level,as.data.frame(ord.fac))
    colnames(print.fac) <- c("levels","ordered")
    if(verbose==TRUE){
        cat("\nCheck: the number of levels for factors and whether they are ordered.\n")
        print(print.fac)
    }
    rm(print.fac)
    rm(data.x)
    
    if(Gorder==2){
        formula <- as.formula(paste(all.vars(formula)[1]," ~ .*.", sep=""))
    }else if(Gorder==3){
        formula <- as.formula(paste(all.vars(formula)[1]," ~ .*.*.", sep=""))
    }
    data <- model.frame(formula,data=data.main)
    rm(data.main)
    
    
    data.orig <- data
    if(diff==TRUE){
        data <- data[order(pair.id),]
        side <- rep(c(1,0),times=nrow(data)/2)
    }else{
        side <- NULL
    }
    
    if(diff==TRUE){
        data1 <- data[side==1,]
        data2 <- data[side==0,]
        X1 <- model.matrixBayes(formula,data=data1)
        X2 <- model.matrixBayes(formula,data=data2)
        X <- X1 - X2 
        y <- model.frame(formula,data=data1)[,1]
        base.name <- colnames(X1)
    }else if(diff==FALSE){
        X <- model.matrixBayes(formula,data=data)
        y <- model.frame(formula,data=data)[,1]
        base.name <- colnames(X)
    }
      
    n.fac <- length(fac.level)
    levelIndex <- CreatelevelIndex(fac.level=fac.level,ord.fac=ord.fac, 
                                   Gorder=Gorder)   
    coef.length <- sum(levelIndex$length)
    coef.length.mu <- coef.length + 1

    ## #########################
    ## Making the basic matrices
    ## Note: these do not change depending on X. Only depends on fac.level and ord.fac
    ## ######################### 
    ## (1) Making L matrix    
    L <- Lcombinefunction(fac.level=fac.level,ord.fac=ord.fac,facCons=FALSE,Gorder=Gorder)    
    ## (2) Making Psy matrix
    PsyConstraintMat <- PsyConstraintCombine(fac.level=fac.level,ord.fac=ord.fac,
                                                 Gorder=Gorder)

    ## (3): Making the weights
    weight.fac <- c()
    weight.ols <- c()
    for(w in 1:n.fac){
        weight.list <- CreateWeights(formula,data=data,Gorder=Gorder,
                                     facCons=FALSE,
                                     side=side,fac.level=fac.level,dif=diff,
                                     ord.fac=ord.fac,type.ind=w)
        weight.fac <- c(weight.fac, weight.list$weight.fac)
        weight.ols <- c(weight.ols, weight.list$weight.ols)        
    }
    weight.fac <- (weight.fac/sum(weight.fac))*length(weight.fac)
    weight <- weight.fac * weight.ols

    ## (4): Create Slack variables
    l.slack <- lengthSlack(fac.level=fac.level,ord.fac=ord.fac,Gorder=Gorder,
                           facCons=FALSE)$l.slack

    ## (5) : Normalize sum constraints so that 0<t<1
    sum.constraint <- c(rep(0,times=coef.length.mu),weight)
    t.sum <- cost*length(weight)
    
    ## (6): Positivity constraint
    Positive <- diag(coef.length.mu+l.slack)
    
    ## (7): Create zero-sum constraint
    Constraint <- CreateANOVAconst(fac.level=fac.level,ord.fac=ord.fac,
                                   Gorder=Gorder,facCons=FALSE)

    ## (8) : Create Z matrix 
    Zcombine <- Zcombinefunction(X=X,fac.level=fac.level,ord.fac=ord.fac, 
                                 Gorder=Gorder)
    ## Adjust Z matrix for slack variables
    ZcombineF <- cbind(Zcombine, matrix(0,ncol=l.slack,nrow=nrow(Zcombine)))
    
    ## (9) : Combine all inequality constraints    
    G <- rbind(-PsyConstraintMat, -sum.constraint, Positive)
    H <- c(rep(0,times=nrow(PsyConstraintMat)),-t.sum,rep(0,times=nrow(Positive)))
    LargeG <- rbind(L, -L, Constraint, -Constraint, G)
    LargeH <- c(rep(-eps, times=nrow(L)),rep(-eps, times=nrow(L)),
                rep(-eps, times=nrow(Constraint)),rep(-eps, times=nrow(Constraint)),
                H)

    solve.fit <- Glsei(A=ZcombineF,B=y,G=LargeG,H=LargeH,verbose=TRUE,
                      tol=5*eps,type="2")
    ## solve.fit <- lsei(A=ZcombineF,B=y,G=G,H=H, E=E, F=F, tol=eps, verbose=TRUE,type="1")
    
    
    BaseCoef <- solve.fit$X[-1]
    coefs <- CoefExtract(BaseCoef,base.name=base.name, 
                         fac.level=fac.level, ord.fac=ord.fac,
                         Gorder=Gorder)
    if(verbose==TRUE){
        cat(paste("\nCheck: cost=",cost,"\n",sep=""))
        ##         cat("\nThe smaller is the cost, the stronger is the regularization.
        ## \nNote that cost=1 corresponds to no regularization.\n")
    }
    output <- list("solve.fit"=solve.fit,
                   "intercept"=solve.fit$X[1],
                   "coefs"=coefs,"data"=data,
                   "nway"=nway,"formula"=formula,
                   "cost"=cost,"fac.level"=fac.level,
                   "ord.fac"=ord.fac,"diff"=diff,
                   "eps"=eps,"side"=side,
                   "data.orig"=data.orig,"formula.orig"=formula.orig,
                   "pair.id"=pair.id)
    
    class(output) <- c("CausalANOVAFit","list")
    return(output)
}
