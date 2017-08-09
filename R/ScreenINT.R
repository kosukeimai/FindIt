ScreenINT <- function(formula, data.main, data.x, screen, screen.type,
                      family, fac.level,screen.num.int, Gorder, maxIter,
                      internal.int, verbose=verbose){

    n.fac <- length(fac.level)
    ## verbose=FALSE, maxIter=maxIte
    term.f <- terms(formula, data=data.main)
    order.f <- attr(term.f, "order")
    all.vars0 <- c(all.vars(term.f)[1], attr(term.f, "term.labels")[order.f==1])
    if(screen==FALSE){        
        if(internal.int==FALSE){
            ## ############################ 
            ## Specified by Users.
            ## ############################ 
            ## Two-ways 
            if(any(order.f==2)==TRUE){
                indTwo.n <- strsplit(attr(term.f, "term.labels")[attr(term.f, "order")==2], "\\:")
                if(verbose==TRUE){
                    select.int <- as.data.frame(do.call("rbind", indTwo.n))
                    colnames(select.int) <- c("Fac.1","Fac.2")
                    cat("\nSpecified Two-way Interactions.\n")
                    print(select.int); rm(select.int)
                }
                indTwo <- matrix(NA, ncol=length(indTwo.n), nrow=2)
                for(j in 1:length(indTwo.n)){
                    indTwo[1,j] <- which(indTwo.n[[j]][1]==all.vars0) - 1
                    indTwo[2,j] <- which(indTwo.n[[j]][2]==all.vars0) - 1
                }
                ## indTwo <- apply(indTwo, 2, function(x) x[order(x)])
                ## indTwo <- indTwo[,order(indTwo[1,])]
                ## indTwo <- indTwo[,order(indTwo[2,])]
            }
            ## Three-ways 
            if(any(order.f==3)==TRUE){
                indThree.n <- strsplit(attr(term.f, "term.labels")[attr(term.f, "order")==3], "\\:")
                indThree <- matrix(NA, ncol=length(indThree.n), nrow=3)
                for(j in 1:length(indThree.n)){
                    indThree[1,j] <- which(indThree.n[[j]][1]==all.vars0) - 1
                    indThree[2,j] <- which(indThree.n[[j]][2]==all.vars0) - 1
                    indThree[3,j] <- which(indThree.n[[j]][3]==all.vars0) - 1
                }
                ## CHECK
                indThreeCheck <- indTwo2Three(indTwo=indTwo, n.fac=n.fac)
                if(is.null(indThreeCheck)){
                    check.three <- FALSE
                }else{
                    check.three <- c()
                    for(z in 1:ncol(indThree)){
                        check.three[z] <- any(apply(indThreeCheck == indThree[,j],2,all))
                    }
                }
                if(any(check.three==FALSE)){
                    stop("\nEvery pair of terms in Three-way Interactions should be specified in 'int2.formula' to satisfy strong hierarchy.\n")
                }
                rm(check.three); rm(indThreeCheck)

                if(verbose==TRUE){
                    select.int3 <- as.data.frame(do.call("rbind", indThree.n))
                    colnames(select.int3) <- c("Fac.1","Fac.2", "Fac.3")
                    cat("\nSpecified Three-way Interactions.\n")
                    print(select.int3); rm(select.int3)
                }
                Gorder <- 3
            }else{ indThree <- NULL}
        }else if(internal.int==TRUE & Gorder==2){
            ## ############################ 
            ## Analyze all Two-ways
            ## ############################ 
            if(verbose==TRUE) cat("\nAnalyzing all two-way interactions...\n")
            formula <- as.formula(paste(all.vars0[1]," ~ .*.", sep=""))
            indTwo <- combn(seq(1:n.fac),2)
            indThree <- NULL
        }else if(internal.int==TRUE & Gorder==3){
            ## ############################ 
            ## Analyze all Three-ways
            ## ############################ 
            if(verbose==TRUE) cat("\nAnalyzing all two-way and three-way interactions...\n")
            formula <- as.formula(paste(all.vars0[1]," ~ .*.*.", sep=""))
            indTwo <- combn(seq(1:n.fac),2)
            indThree <- combn(seq(1:n.fac),3)
        }
    }else if(screen==TRUE){
        ## ############################ 
        ## Data-driven Screening 
        ## ############################ 
        ## prepare the data for screening
        data.s <- matrix(NA, ncol=ncol(data.x),nrow=nrow(data.x))
        for(i in 1:ncol(data.s)){
            data.s[,i] <- as.numeric(data.x[,i]) - 1            
        }
        ## Do the screening here and modify formula.
        main.for <- paste(all.vars0[1],"~",paste(all.vars0[-1],collapse="+"),sep="")
        if(verbose==TRUE) cat("\n ***** Screening ***** \n")
        
        if(screen.type=="cv.min" | screen.type=="cv.1Std"){
            gl.screen <- glinternet.cv(X=data.s, Y=data.main[,1], family=family, numLevels=fac.level,
                                       nLambda=5, nFolds=10,verbose=FALSE, maxIter=maxIter)            
            if(screen.type=="cv.min"){
                if(is.null(gl.screen$activeSet[[1]]$catcat)==TRUE) indTwo <- NULL
                else indTwo <- t(gl.screen$activeSet[[1]]$catcat)
            }else if(screen.type=="cv.1Std"){
                if(is.null(gl.screen$activeSet1Std[[1]]$catcat)==TRUE) indTwo <- NULL
                else indTwo <- t(gl.screen$activeSet1Std[[1]]$catcat)
            }
            
        }else if (screen.type=="fixed"){
            gl.screen <- glinternet(X=data.s, Y=data.main[,1], family=family, numLevels=fac.level,
                                    numToFind=screen.num.int,
                                    verbose=FALSE, maxIter=maxIter)
            check.indTwo <- gl.screen$activeSet[[length(gl.screen$activeSet)]]$catcat
            if(is.null(check.indTwo)==FALSE){
                indTwo <- t(gl.screen$activeSet[[length(gl.screen$activeSet)]]$catcat)
            }else{
                indTwo <- NULL
            }
        }
        if(is.null(indTwo)==FALSE){
            if(verbose==TRUE){
                select.int <- as.data.frame(cbind(all.vars0[(indTwo[1,]+1)], all.vars0[(indTwo[2,]+1)]))
                colnames(select.int) <- c("Fac.1","Fac.2")
                cat("\nSelected Two-way Interactions.\n")
                print(select.int); rm(select.int)
            }
            
            int.for <- c()
            for(j in 1:ncol(indTwo)){
                int.for[j] <- paste(all.vars0[indTwo[1,j]+1], all.vars0[indTwo[2,j]+1], sep=":")
            }
            int.for <- paste(int.for, collapse="+")
        }else{
            if(verbose==TRUE) cat("\nNo Two-way Interaction is selected.\n")
            formula <- as.formula(main.for)
            indThree <- NULL
            Gorder <- 1
        }
        
        if(Gorder==3 & is.null(indTwo)==FALSE){
            indThree <- indTwo2Three(indTwo=indTwo, n.fac=n.fac)
            if(is.null(indThree)==FALSE){
                if(verbose==TRUE){
                    select.int3 <- as.data.frame(cbind(all.vars0[(indThree[1,]+1)],
                                                       all.vars0[(indThree[2,]+1)],
                                                       all.vars0[(indThree[3,]+1)]))
                    colnames(select.int3) <- c("Fac.1","Fac.2", "Fac.3")
                    cat("\nSelected Three-way Interactions.\n")
                    print(select.int3); rm(select.int3)
                }
                Gorder <- 3
                ## Formula
                int3.for <- c()
                for(j in 1:ncol(indThree)){
                    int3.for[j] <- paste(all.vars0[indThree[1,j]+1],
                                         all.vars0[indThree[2,j]+1],
                                         all.vars0[indThree[3,j]+1],sep=":")
                }
                int3.for <- paste(int3.for, collapse="+")
                formula <- as.formula(paste(main.for, int.for, int3.for, sep="+"))
            }else if(is.null(indTwo)==FALSE & is.null(indThree)==TRUE){
                if(verbose==TRUE) cat("\nNo Three-way Interaction is selected.\n")
                formula <- as.formula(paste(main.for, int.for, sep="+"))
                Gorder <- 2
            }
        }else if(Gorder==2 & is.null(indTwo)==FALSE){
            formula <- as.formula(paste(main.for, int.for, sep="+"))
            indThree <- NULL
            Gorder <- 2
        }
    }    
    output <- list("formula"=formula, "indTwo"=indTwo, "indThree"=indThree, "Gorder"=Gorder)
    return(output)
}
