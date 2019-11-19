#' Summarizing CausalANOVA output
#' @param object An object from \code{CausalANOVA}
#' @param digit the number of digits
#' @param verbose report additional summary results
#' @param verbose.full report full summary results
#' @param ... Other parameters
#' @method summary CausalANOVA
#' @export
summary.CausalANOVA <-function(object, digit=4, verbose=TRUE, verbose.full=TRUE, ...){

    if(class(object)[2]=="stab"){
        fit <- object$fit
        stab.fit <- object$stab.fit
    }else{
        fit <- object
    }
    
    data <- fit$data
    Gorder <- fit$Gorder
    formula <- fit$formula
    formula.orig <- fit$formula.orig
    cost <- fit$cost
    n.fac <- length(fit$fac.level)
    fac.level <- fit$fac.level
    terms.f <- terms(formula, data=data)
    order.f <- attr(terms.f, "order")
    Name <- attr(terms.f, "term.labels")
    all.vars0 <- Name[order.f==1]
    inference <- fit$inference
    Stab <- as.numeric(class(object)[2]=="stab")
    level.list <- fit$level.list
    indTwo <- fit$indTwo
    indThree <- fit$indThree
    collapse <- fit$collapse
    screening <- fit$screen
        
    ## ################# 
    ## Compute Range
    ## ################# 
    range.reg <- c()
    for(z in 1:length(fit$coef)){
        range.reg[z] <-  max(fit$coef[[z]]) - min(fit$coef[[z]])
    }
    range <- round(range.reg,digits=3)

    Range <-as.data.frame(range)
    rownames(Range) <- Name

    ## Name Preparation 
    Factor.ame <- rep(Name[order.f==1], fac.level)
    Level.ame <- unlist(level.list)
    Base.ame <- unlist(lapply(fac.level, function(x) c(rep("",(x-1)), "***")))
    index.ame<- data.frame(cbind(Factor.ame, Level.ame, Base.ame))

    if(Gorder >=2){
        Name.amie2 <- Name[order.f==2]
        Factor.amie2 <- list()
        Level.amie2 <- list()
        Base.amie2 <- c()
        for(z in 1:sum(order.f==2)){
            size <- fac.level[[indTwo[1,z]]]*fac.level[[indTwo[2,z]]]
            Factor.amie2[[z]] <- rep(Name.amie2[z], times=size)
            Level.amie2[[z]] <- expand.grid(Var1=level.list[[indTwo[1,z]]], Var2=level.list[[indTwo[2,z]]])
            Base.amie2 <- c(Base.amie2, c(rep("", size-1), "***"))
        }
        Factor.amie2 <- unlist(Factor.amie2)
        Level.amie2 <- do.call("rbind",Level.amie2)
        index.amie2<- data.frame(cbind(Factor.amie2, Level.amie2, Base.amie2))

        
        if(Gorder == 3){
           Name.amie3 <- Name[order.f==3]
           Factor.amie3 <- list()
           Level.amie3 <- list()
           Base.amie3 <- c()
           for(z in 1:sum(order.f==3)){
               size <- fac.level[[indThree[1,z]]]*fac.level[[indThree[2,z]]]*fac.level[[indThree[3,z]]]
               Factor.amie3[[z]] <- rep(Name.amie3[z], times=size)
               Level.amie3[[z]] <- expand.grid(Var1=level.list[[indThree[1,z]]],
                                               Var2=level.list[[indThree[2,z]]],
                                               Var3=level.list[[indThree[3,z]]])
               Base.amie3 <- c(Base.amie3, c(rep("", size-1), "***"))
           }
           Factor.amie3 <- unlist(Factor.amie3)
           Level.amie3 <- do.call("rbind",Level.amie3)
           index.amie3<- data.frame(cbind(Factor.amie3, Level.amie3, Base.amie3))
        }                        
    }
    ## ################
    ## AME
    ## ################
    if(inference==TRUE){
        AME.CI <- fit$CI.table[order.f==1]
        AME.out <- do.call("rbind", AME.CI)
        AME.table <- round(AME.out, digit)
        AME.table2 <- data.frame(AME.table)
                       
        AME.table.pr <- cbind(index.ame, AME.table2)
        colnames(AME.table.pr) <- c("Factor","Levels","base","AME","Std.Err","2.5%CI","97.5%CI")

        if(Gorder>=2){                        
            AMIE2.CI <- fit$CI.table[order.f==2]
            AMIE2.out <- do.call("rbind", AMIE2.CI)
            AMIE2.table <- round(AMIE2.out, digit)

            AMIE2.table.pr <- cbind(index.amie2, data.frame(AMIE2.table))
            colnames(AMIE2.table.pr) <- c("Factor","Level1","Level2", "base","AMIE","Std.Err","2.5%CI","97.5%CI")

            if(Gorder == 3){
                AMIE3.CI <- fit$CI.table[order.f==3]
                AMIE3.out <- do.call("rbind", AMIE3.CI)
                AMIE3.table <- round(AMIE3.out, digit)

                AMIE3.table.pr <- cbind(index.amie3, data.frame(AMIE3.table))
                colnames(AMIE3.table.pr) <- c("Factor","Level1","Level2","Level3",
                                              "base","AMIE","Std.Err","2.5%CI","97.5%CI")
            }else{
                AMIE3.table.pr <- NULL
            }
        }else{
            AMIE2.table.pr <- NULL
            AMIE3.table.pr <- NULL
        }
        
    }else if(inference==FALSE){
        AME.CI <- fit$coefs[order.f==1]
        AME.out0 <- list()
        for(z in 1:n.fac){
            AME.out0[[z]] <-  AME.CI[[z]] - AME.CI[[z]][length(AME.CI[[z]])]
        }
        AME.table <- data.frame(matrix(round(unlist(AME.out0),digit), ncol=1))

        AME.table.pr <- cbind(index.ame, AME.table)
        colnames(AME.table.pr) <- c("Factor","Levels","base","AME")

        if(Gorder>=2){
            AMIE2.CI <- fit$coefs[order.f==2]
            AMIE2.out0 <- list()
            for(z in 1:length(AMIE2.CI)){
                AMIE2.out0[[z]] <-  AMIE2.CI[[z]] - AMIE2.CI[[z]][length(AMIE2.CI[[z]])]
            }
            AMIE2.table <- data.frame(matrix(round(unlist(AMIE2.out0),digit), ncol=1))
            AMIE2.table.pr <- cbind(index.amie2, data.frame(AMIE2.table))
            colnames(AMIE2.table.pr) <- c("Factor","Level1","Level2", "base","AMIE")
            
            if(Gorder>=3){
                AMIE3.CI <- fit$coefs[order.f==3]
                AMIE3.out0 <- list()
                for(z in 1:length(AMIE3.CI)){
                    AMIE3.out0[[z]] <-  AMIE3.CI[[z]] - AMIE3.CI[[z]][length(AMIE3.CI[[z]])]
                }
                AMIE3.table <- data.frame(matrix(round(unlist(AMIE3.out0),digit), ncol=1))
                AMIE3.table.pr <- cbind(index.amie3, data.frame(AMIE3.table))
                colnames(AMIE3.table.pr) <- c("Factor","Level1","Level2","Level3", "base","AMIE")
            }else{
                AMIE3.table.pr <- NULL
            }           
        }else{
            AMIE2.table.pr <- NULL
            AMIE3.table.pr <- NULL
        }

        if(Stab==TRUE){
            ## Range 
            Range$select.prob <- round(stab.fit$range.stab[,2], digits=digit)
            colnames(Range)[2] <- "select.prob"

            ## AME
            AME.stab <- stab.fit$AME.stab[,2]      
            AME.table.pr$select.prob <- AME.stab            
        }
    }


    Range.pr <- as.data.frame(Range)

    ## Collapsing Results
    if(collapse==TRUE){
        collapse.u  <- Collapsing(fit)$collapse.level
        drop.fac <- unlist(lapply(collapse.u, function(x) all(x==1)))                        
        ## New Level Names
        collapse.name <- list()
        for(i in 1:n.fac){
            original.level <- level.list[[i]]
            new.name <- c()
            for(j in 1:length(unique(collapse.u[[i]]))){
                new.name[j] <- paste(original.level[collapse.u[[i]]==j],collapse="/")
            }
            collapse.name[[i]] <- unique(new.name[collapse.u[[i]]])
        }
        names(collapse.name) <- Name[order.f==1]
    }
    if(screening==TRUE){
        if(is.null(indTwo)==FALSE){
            select.int2 <- as.data.frame(cbind(all.vars0[(indTwo[1,])], all.vars0[(indTwo[2,])]))
            colnames(select.int2) <- c("Fac.1","Fac.2")
        }
        if(is.null(indTwo)==FALSE){
            select.int3 <- as.data.frame(cbind(all.vars0[(indThree[1,])], all.vars0[(indThree[2,])], all.vars0[(indThree[3,])]))
            colnames(select.int3) <- c("Fac.1","Fac.2", "Fac.3")
        }
    }
    
    
    if(verbose==TRUE){
        cat(" \nModel:\n ")
        print(formula)               
        
        cat("******\n")
        cat("Range:\n")
        cat("******\n")
        print(Range.pr)
        
        cat("*******************************\n")
        cat("Average Marginal Effects (AME):\n")
        cat("*******************************\n")
        print(AME.table.pr, row.names=F)
        cat("\n")
        
        if(verbose.full==TRUE & is.null(AMIE2.table.pr)==FALSE){
            cat("********************************************\n")
            cat("Average Marginal Interaction Effects (AMIE):\n")
            cat("********************************************\n")
            cat("Two-way Interactions:\n")
            print(AMIE2.table.pr, row.names=FALSE)

            if(is.null(AMIE3.table.pr)==FALSE){
                cat("\nThree-way Interactions:\n")
                print(AMIE3.table.pr, row.names=FALSE)
            }
        }

        if(screening==TRUE){
            cat("--------------------------------------------\n")
            cat("Selected Interactions; Results of Screening:\n")
            cat("--------------------------------------------\n")
            if(is.null(indTwo)==FALSE){
                print(select.int2)
            }
            if(is.null(indThree)==FALSE){
                print(select.int3)
            }
        }
            
        if(collapse==TRUE){
            cat("-----------------------------\n")
            cat("Results of Collapsing Levels:\n")
            cat("-----------------------------\n")
            print(collapse.name)           
        }

        if(inference==FALSE & Stab==FALSE){
            cat("===============================================================\n")
            cat(": To get confidence intervals, Use 'test.CausalANOVA' function:\n")
            cat("===============================================================\n")
        }
        
    }    
    output <- list("range"=Range.pr,"range.name"=Name,"AME"=AME.table.pr, "AMIE2"=AMIE2.table.pr,"AMIE3"=AMIE3.table.pr)
}