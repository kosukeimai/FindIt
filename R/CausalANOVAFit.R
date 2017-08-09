CausalANOVAFit <- function(formula, internal.int,
                           data,cost,pair.id=NULL,
                           nway=2, diff=TRUE, eps=1e-5, family="binomial",
                           fac.level=NULL,ord.fac=NULL, cluster=cluster,
                           verbose=TRUE,
                           screen=FALSE, screen.type="fixed", screen.num.int=5,
                           collapse=FALSE, collapse.type="fixed", collapse.cost=0.3,
                           maxIter=100){

    ## ############################
    ## (0.1) Setup
    ## ############################ 
    Gorder <- nway
    formula.orig <- formula
    ## We convert the data into the data only with necessary variables. 
    data.main <- model.frame(formula, data=data)
    data.x <- model.frame(formula, data=data)[,-1]
    
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
        print(print.fac);rm(print.fac)
    }



    ## ############################
    ## (0.2) Screening
    ## ############################
    if(Gorder!=1){
        screen.out <- ScreenINT(formula=formula, data.main=data.main,
                                data.x=data.x, screen=screen, screen.type=screen.type,
                                family=family, fac.level=fac.level, screen.num.int=screen.num.int,
                                Gorder=Gorder, maxIter=maxIter, internal.int=internal.int,
                                verbose=verbose)
        
        indTwo <-  screen.out$indTwo
        indThree <- screen.out$indThree
        formula <-  screen.out$formula
        Gorder <- screen.out$Gorder
        
        if(screen==FALSE){
            formula.orig <- formula
        }
    }else{
        indTwo <-  NULL
        indThree <- NULL
    }
    
    rm(data.x); 
    data <- model.frame(formula,data=data.main); rm(data.main)
    level.list <- lapply(data[,-1], levels)
    data.orig <- data
    order.f <- attr(terms(formula,data=data), "order")

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
        cluster <- cluster[side==1]
        X1 <- model.matrixBayes(formula,data=data1)
        X2 <- model.matrixBayes(formula,data=data2)
        X <- X1 - X2 
        y <- model.frame(formula,data=data1)[,1]
        base.name <- colnames(X1)
        if(cost==1){
            model.frameC  <- model.frame(formula,data=data1)
            contr <- rep(list("contr.sum"), ncol(model.frameC) - 1)
            names(contr) <- colnames(model.frameC)[-1]
            x1 <- model.matrix(formula, data=data1,contrast=contr)[,-1]
            x2 <- model.matrix(formula, data=data2,contrast=contr)[,-1]
            X.no <- x1-x2; rm(model.frameC)
        }
    }else if(diff==FALSE){
        X <- model.matrixBayes(formula,data=data)
        y <- model.frame(formula,data=data)[,1]
        base.name <- colnames(X)
        if(cost==1){
            model.frameC  <- model.frame(formula,data=data)
            contr <- rep(list("contr.sum"), ncol(model.frameC) - 1)
            names(contr) <- colnames(model.frameC)[-1]
            X.no <- model.matrix(formula, data=data,contrast=contr)[,-1]
            rm(model.frameC)
        }
    }

    ## ############################ 
    ## (0.4) If No Regularization
    ## ############################ 
    if(cost==1 & screen==FALSE){
        ## if(verbose==TRUE) cat("\nNo Collapsing when 'collapse=FALSE'.\n")      
        NoReg <- NoRegularization(y=y, X.no=X.no, base.name=base.name, fac.level=fac.level,
                                  ord.fac=ord.fac, Gorder=Gorder, indTwo=indTwo,indThree=indThree,
                                  cluster=cluster)                
        intercept <- NoReg$intercept
        coefs <- NoReg$coefs
        vcov  <- NoReg$vcov

        ## Index
        levelIndex <- CreatelevelIndex(fac.level=fac.level, ord.fac=ord.fac, Gorder=Gorder,
                                       indTwo=indTwo, indThree=indThree)
        use.ind <-  (levelIndex$plus==1)*(levelIndex$dif==0)
        Index.use <- levelIndex[use.ind==1,]
        Index.use$start <- c(1,cumsum(Index.use$length)[-nrow(Index.use)]+1)
        Index.use$end <-  cumsum(Index.use$length)
        
        order.f <- attr(terms(formula, data=data), "order")
        
        CI.tab <- list()
        for(z in 1:length(coefs)){
            seq.u <- Index.use$start[z]:Index.use$end[z]
            effect <- coefs[[z]] - coefs[[z]][length(coefs[[z]])]
            vcov.p <- c()
            for(i in 1:(length(seq.u)-1)){
                vcov.p[i] <- VarEffect(ind1=seq.u[i],ind2=seq.u[length(seq.u)],vcov.full=vcov)
            }
            sd <- c(sqrt(round(vcov.p,digits=10)),0)
            table <- as.data.frame(cbind(effect, sd, effect - 1.96*sd, effect + 1.96*sd))
            if(Index.use$order[z]==1){
                colnames(table) <- c("AME", "sd", "2.5% CI", "97.5% CI")
            }else{
                colnames(table) <- c("AMIE", "sd", "2.5% CI", "97.5% CI")
            }
            CI.tab[[z]] <- table
        }
        
        AME <- coefs[order.f==1]
        AME <- lapply(AME, function(x) round(x/(eps*10))*(eps*10))
        if(Gorder>=2){
            AMIE2 <- coefs[order.f==2]
            AMIE2 <- lapply(AMIE2, function(x) round(x/(eps*10))*(eps*10))
        }else{
            AMIE2 <- NULL
        }
        if(Gorder==3){
            AMIE3 <- coefs[order.f==3]
            AMIE3 <- lapply(AMIE3, function(x) round(x/(eps*10))*(eps*10))
        }else{
            AMIE3 <- NULL
        }
        
        output <- list("intercept"=intercept,
                       "coefs"=coefs,"data"=data,
                       "vcov"=vcov, "CI.table"=CI.tab,
                       "AME"=AME, "AMIE2"=AMIE2,"AMIE3"=AMIE3,
                       "nway"=nway,"formula"=formula,
                       "cost"=cost,"fac.level"=fac.level,
                       "ord.fac"=ord.fac,
                       "level.list"=level.list,
                       "diff"=diff,
                       "eps"=eps,"side"=side,
                       "data.orig"=data.orig,
                       "formula.orig"=formula.orig,
                       "pair.id"=pair.id,
                       "indTwo"=indTwo,
                       "indThree"=indThree,
                       "family"=family,
                       "internal.int"=internal.int,
                       "Gorder"=Gorder,
                       "inference"=TRUE,
                       "screen"=screen,
                       "screen.type"=screen.type,
                       "screen.num.int"=screen.num.int,
                       "collapse"=collapse,
                       "collapse.type"=collapse.type,
                       "collapse.cost"=collapse.cost)
        
        class(output) <- c("CausalANOVAFit","list")
        return(output)
    }
    
    ## ##################
    ## (1.0) Create the basic Index
    ## ################## 
    levelIndex <- CreatelevelIndex(fac.level=fac.level,ord.fac=ord.fac, 
                                   Gorder=Gorder, indTwo=indTwo, indThree=indThree)   
    coef.length <- sum(levelIndex$length)
    coef.length.mu <- coef.length + 1


    ## #########################
    ## Making the basic matrices
    ## Note: these do not change depending on X. Only depends on fac.level, ord.fac, indTwo, indThree
    ## ######################### 
    ## (1.1) Making L matrix    
    L <- Lcombinefunction(fac.level=fac.level,ord.fac=ord.fac,facCons=FALSE,Gorder=Gorder,
                          indTwo=indTwo,indThree=indThree)

    
    ## (1.2) Making Psy matrix
    PsyConstraintMat <- PsyConstraintCombine(fac.level=fac.level,ord.fac=ord.fac,
                                             Gorder=Gorder, 
                                             indTwo=indTwo, indThree=indThree)

    ## (1.3): Making the weights
    weight.fac <- c()
    weight.ols <- c()
    for(w in 1:n.fac){
        weight.list <- CreateWeights(formula,data=data,Gorder=Gorder,
                                     facCons=FALSE,
                                     side=side,fac.level=fac.level,dif=diff,
                                     ord.fac=ord.fac,type.ind=w,
                                     indTwo=indTwo, indThree=indThree)
        weight.fac <- c(weight.fac, weight.list$weight.fac)
        weight.ols <- c(weight.ols, weight.list$weight.ols)        
    }
    ## weight.fac <- (weight.fac/sum(weight.fac))*length(weight.fac)
    ## weight.fac <- (weight.fac/sum(weight.fac))*length(weight.fac)
    weight <- weight.fac * weight.ols

    ## (1.4): Create Slack variables
    l.slack <- lengthSlack(fac.level=fac.level,ord.fac=ord.fac,Gorder=Gorder,
                           facCons=FALSE)$l.slack

    ## (1.5) : Normalize sum constraints so that 0<t<1
    sum.constraint <- c(rep(0,times=coef.length.mu), weight)
    t.sum <- cost*sum(weight.fac)

    ## (1.6): Positivity constraint
    Positive <- diag(coef.length.mu+l.slack)
    
    ## (1.7): Create zero-sum constraint
    Constraint <- CreateANOVAconst(fac.level=fac.level,ord.fac=ord.fac,
                                   Gorder=Gorder,facCons=FALSE,
                                   indTwo=indTwo, indThree=indThree)

    ## (1.8) : Create Z matrix (Expand version of Design Matrix)
    Zcombine <- Zcombinefunction(X=X,fac.level=fac.level,ord.fac=ord.fac, 
                                 Gorder=Gorder,
                                 indTwo=indTwo, indThree=indThree)
    ## Adjust Z matrix for slack variables
    ZcombineF <- cbind(Zcombine, matrix(0,ncol=l.slack,nrow=nrow(Zcombine)))
    
    ## (1.9) : Combine all inequality constraints    
    G <- rbind(-PsyConstraintMat, -sum.constraint, Positive)
    H <- c(rep(0,times=nrow(PsyConstraintMat)),-t.sum,rep(0,times=nrow(Positive)))
    LargeG <- rbind(L, -L, Constraint, -Constraint, G)
    LargeH <- c(rep(-eps, times=nrow(L)),rep(-eps, times=nrow(L)),
                rep(-eps, times=nrow(Constraint)),rep(-eps, times=nrow(Constraint)),
                H)
    
    if(collapse==TRUE){if(verbose==TRUE) cat("\n ***** Collapsing ***** \n")}

    ## ############################
    ## (2.0) Fit the Main Model
    ## ############################ 
    solve.fit <- Glsei(A=ZcombineF,B=y,G=LargeG,H=LargeH,verbose=TRUE,
                       tol=5*eps,type="2")
        
    BaseCoef <- solve.fit$X[-1]
    coefs <- CoefExtract(BaseCoef,base.name=base.name, 
                         fac.level=fac.level, ord.fac=ord.fac,
                         Gorder=Gorder,indTwo=indTwo,indThree=indThree)

    AME <- coefs[order.f==1]
    AME <- lapply(AME, function(x) round(x/(eps*10))*(eps*10))
    if(Gorder>=2){
        AMIE2 <- coefs[order.f==2]
        AMIE2 <- lapply(AMIE2, function(x) round(x/(eps*10))*(eps*10))
    }else{
        AMIE2 <- NULL
    }
    if(Gorder==3){
        AMIE3 <- coefs[order.f==3]
        AMIE3 <- lapply(AMIE3, function(x) round(x/(eps*10))*(eps*10))
    }else{
        AMIE3 <- NULL
    }
    
    output <- list("solve.fit"=solve.fit,
                   "intercept"=solve.fit$X[1],
                   "coefs"=coefs,"data"=data,
                   "AME"=AME, "AMIE2"=AMIE2,"AMIE3"=AMIE3,
                   "nway"=nway,"formula"=formula,
                   "cost"=cost,"fac.level"=fac.level,
                   "level.list"=level.list,
                   "ord.fac"=ord.fac,"diff"=diff,
                   "eps"=eps,"side"=side,
                   "data.orig"=data.orig,
                   "formula.orig"=formula.orig,
                   "pair.id"=pair.id,
                   "indTwo"=indTwo,
                   "indThree"=indThree,
                   "family"=family,
                   "screen"=screen,
                   "internal.int"=internal.int,
                   "screen.type"=screen.type,
                   "screen.num.int"=screen.num.int,
                   "Gorder"=Gorder,
                   "inference"=FALSE,
                   "collapse"=collapse,
                   "collapse.type"=collapse.type,
                   "collapse.cost"=collapse.cost)
    
    class(output) <- c("CausalANOVAFit","list")
    return(output)
}
