CausalANOVA <- function(formula, int2.formula=NULL, int3.formula=NULL,
                        data, nway=1,
                        pair.id=NULL, diff=FALSE,
                        screen=FALSE, screen.type="fixed", screen.num.int=3,
                        collapse=FALSE, collapse.type="fixed", collapse.cost=0.3,
                        family="binomial", cluster=NULL,
                        maxIter=50, eps=1e-5,
                        fac.level=NULL, ord.fac=NULL,
                        select.prob=FALSE, boot=100, seed=1234,
                        verbose=TRUE){    

    ## House Keeping
    if(collapse==FALSE){
        cost <- 1
    }else if(collapse.type=="fixed"){
        cost <- collapse.cost
        if(cost > 1 | cost < 0 ){
            stop("Specify 'collapse.cost' between 0 and 1")
        }
    }
    if(missing(fac.level)) fac.level <- NULL
    if(missing(ord.fac)) ord.fac <- NULL        
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
    if(collapse.type=="cv.min" | collapse.type=="cv.1Std"){
        cat("\nRunning Cross Validation might take time...\n")
        cv <- cv.CausalANOVA(formula = formula,
                             int2.formula=int2.formula, int3.formula=int3.formula,
                             data=data,
                             cv.collapse.cost=c(0.1,0.3,0.7),
                             nfolds=5,
                             pair.id=pair.id, diff=diff,
                             nway=nway, family=family,
                             seed=seed, cluster=cluster,
                             screen=screen, screen.type=screen.type,
                             screen.num.int=screen.num.int,
                             maxIter=maxIter, eps=eps,                                                         
                             fac.level=fac.level,ord.fac=ord.fac, verbose=verbose)                            
        if(collapse.type=="cv.min") cost <- cv$cv.min
        if(collapse.type=="cv.1Std") cost <- cv$cv.1Std
        print(paste("Selected Cost parameter=", cost, sep=""))
    }
    if(cost==1 & select.prob==TRUE){
        select.prob <- FALSE
    }

    ## Factors only
    all.fac <- all(unlist(lapply(model.frame(formula,data=data)[,-1],
                                 FUN=function(x) is.element("factor",class(x)))))
    if(all.fac==FALSE) stop("Design matrix should contain only factors.")
    rm(all.fac);

    main.formula <- formula
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
    
    ## ############################
    ## Main Function
    ## ############################ 
    fit <- CausalANOVAFit(formula=formula,
                          internal.int=internal.int,
                          data=data,cost=cost,pair.id=pair.id,
                          family=family, cluster=cluster,
                          screen=screen,screen.type=screen.type,
                          screen.num.int=screen.num.int,
                          maxIter=maxIter,
                          nway=nway,diff=diff,eps=eps,
                          collapse=collapse, collapse.type=collapse.type,
                          collapse.cost=collapse.cost,
                          fac.level=fac.level,ord.fac=ord.fac,
                          verbose=verbose)

    fit <- c(fit, main.formula=main.formula, int2.formula=int2.formula, int3.formula=int3.formula)
    
    if(select.prob==TRUE){
        stab.fit <- stab.CausalANOVA(fit,cluster=cluster,boot=boot)
        output <- list("fit"=fit, "stab.fit"=stab.fit)
        class(output) <- c("CausalANOVA","stab","list")
    }else{
        output <- fit
        class(output) <- c("CausalANOVA","list")
    }        
    return(output)
}
