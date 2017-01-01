AMIE <-function(object,fac.name,level.name,base.name=NULL,verbose=TRUE,...){

    if(class(object)[2]=="stab"){
        fit <- object$fit
        stab.fit <- object$stab.fit
        stability <- TRUE
    }else{
        fit <- object
    }

    if(fit$nway!=2){
        warning("AMIE currently works only for nway=2.")
    }

    digit <- nchar(fit$eps) - 1
    
    output <- AMIEFit(object=object,fac.name=fac.name,
                      level.name=level.name,
                      base.name=base.name,
                      verbose=FALSE)

    
    if(class(object)[2]!="stab"){        
        baseline <- output$baseline
        decom <- output$decom
    }
    
    if(class(object)[2]=="stab"){
        data <- fit$data        
        formula <- fit$formula        

        AME <- stab.fit$AME.stab.v
        AMIE <- stab.fit$AMIE.stab.v
        
        data.x <- model.frame(formula,data=data)[,-1]
        fac.ind <- which(colnames(data.x) %in% fac.name)
        
        ## AMIE
        ind.fac1 <- as.numeric(regexpr(fac.name[1],rownames(AMIE))>0)
        ind.fac2 <- as.numeric(regexpr(fac.name[2],rownames(AMIE))>0)
        ## AME
        ind.AME.fac1 <- as.numeric(regexpr(fac.name[1],rownames(AME))>0)
        ind.AME.fac2 <- as.numeric(regexpr(fac.name[2],rownames(AME))>0)
        
        if(is.null(base.name)){              
            baseline <- c(tail(levels(data.x[,fac.ind[1]]),1),tail(levels(data.x[,fac.ind[2]]),1))
            AMIE.base <- 0
            AME.base1 <- 0
            AME.base2 <- 0
            names(baseline) <- colnames(data.x)[fac.ind]
        }else{
            ## AMIE
            ind.base1 <- as.numeric(regexpr(base.name[1],rownames(AMIE))>0)
            ind.base2 <- as.numeric(regexpr(base.name[2],rownames(AMIE))>0)
            ind.base <- ind.fac1*ind.fac2*ind.base1*ind.base2 
            AMIE.base <- AMIE[ind.base==1,]
            
            ## AME
            ind.AME.base1t <- as.numeric(regexpr(base.name[1],rownames(AME))>0)
            ind.AME.base2t <- as.numeric(regexpr(base.name[2],rownames(AME))>0)
            ind.AME.base1 <- ind.AME.fac1*ind.AME.base1t
            ind.AME.base2 <- ind.AME.fac2*ind.AME.base2t
            AME.base1 <- AME[ind.AME.base1==1,]
            AME.base2 <- AME[ind.AME.base2==1,]
            
            baseline <- base.name
            names(baseline) <- fac.name
        }
        
        ## Name match
        ind.lev1 <- as.numeric(regexpr(level.name[1],rownames(AMIE))>0)
        ind.lev2 <- as.numeric(regexpr(level.name[2],rownames(AMIE))>0)
        ind <- ind.fac1*ind.fac2*ind.lev1*ind.lev2
        AMIE.main.s <- AMIE[ind==1,] - AMIE.base
        
        ind.AME.lev1 <- as.numeric(regexpr(level.name[1],rownames(AME))>0)
        ind.AME.lev2 <- as.numeric(regexpr(level.name[2],rownames(AME))>0)
        ind.AME1 <- ind.AME.fac1*ind.AME.lev1
        ind.AME2 <- ind.AME.fac2*ind.AME.lev2
        AME.main1.s <- AME[ind.AME1==1,] - AME.base1
        AME.main2.s <- AME[ind.AME2==1,] - AME.base2
     
        ACE <- AMIE.main.s + AME.main1.s + AME.main2.s

        ACE <- round(ACE,digits=digit)
        AMIE.main.s <- round(AMIE.main.s,digits=digit)
        AME.main1.s <- round(AME.main1.s,digits=digit)
        AME.main2.s <- round(AME.main2.s,digits=digit)
        
        AME1.name <- paste("AME (",fac.name[1],")",sep="")
        AME2.name <- paste("AME (",fac.name[2],")",sep="")
        prop.AME1.name <- paste("Prop.",AME1.name,sep="")
        prop.AME2.name <- paste("Prop.",AME2.name,sep="")

        ## Compute Stability
        ACE.main <- round(output$decom[1,1],digits=digit)
        ACE.stab.t <- mean(sign(ACE) ==sign(ACE.main))
        if(ACE.main==0) ACE.stab.t <- 1 - ACE.stab.t
        ACE.stab <- c(ACE.main,ACE.stab.t)
        ## ACE.stab <- c(ACE.main,mean(sign(round(ACE,digits=digit)) 
        ##                             ==sign(round(ACE.main,digits=digit))))
        

        AMIE.main.orig <- round(output$AMIE.main,digits=digit)
        AMIE.stab.t <- mean(sign(AMIE.main.s) ==sign(AMIE.main.orig))
        if(AMIE.main.orig==0) AMIE.stab.t <- 1 - AMIE.stab.t
        AMIE.stab <- c(AMIE.main.orig,AMIE.stab.t)
        
        if(any(is.element(baseline,level.name)==TRUE)){
            if(all(AME.main1.s==0)){
                AME.main2.orig <- round(output$decom[3,1],digits=digit)
                AME2.stab.t <- mean(sign(AME.main2.s) ==sign(AME.main2.orig))
                if(AME.main2.orig==0) AME2.stab.t <- 1 - AME2.stab.t
                AME2.stab <- c(AME.main2.orig,AME2.stab.t)
                stab.prob <- c(ACE.stab[2],AMIE.stab[2],AME2.stab[2],".",".")
                estimate <- c(ACE.main,AMIE.main.orig,AME.main2.orig,
                              AMIE.main.orig/ACE.main, AME.main2.orig/ACE.main)
                estimate <- round(estimate,digits=3)
                decom <- cbind(estimate,as.data.frame(stab.prob))
                decom <- as.data.frame(decom)
                rownames(decom) <- c("Conditional Effect","AMIE",AME2.name,
                                     "Prop.AMIE",prop.AME2.name)
            }else{              
                AME.main1.orig <- round(output$decom[3,1],digits=digit)
                AME1.stab.t <- mean(sign(AME.main1.s) ==sign(AME.main1.orig))
                if(AME.main1.orig==0) AME1.stab.t <- 1 - AME1.stab.t
                AME1.stab <- c(AME.main1.orig,AME1.stab.t)
                ## AME1.stab <- c(AME.main1.orig,mean(sign(AME.main1.s) ==sign(AME.main1.orig)))
                stab.prob <- c(ACE.stab[2],AMIE.stab[2],AME1.stab[2],".",".")
                estimate <- c(ACE.main,AMIE.main.orig,AME.main1.orig,
                              AMIE.main.orig/ACE.main, AME.main1.orig/ACE.main)
                                estimate <- round(estimate,digits=3)
                decom <- cbind(estimate,as.data.frame(stab.prob))
                rownames(decom) <- c("Conditional Effect","AMIE",AME1.name,
                                     "Prop.AMIE",prop.AME1.name)
            }
        }else{
            AME.main1.orig <- round(output$decom[3,1],digits=digit)
            AME.main2.orig <- round(output$decom[4,1],digits=digit)
            AME2.stab.t <- mean(sign(AME.main2.s) ==sign(AME.main2.orig))
            if(AME.main2.orig==0) AME2.stab.t <- 1 - AME2.stab.t
            AME2.stab <- c(AME.main2.orig,AME2.stab.t)
            AME1.stab.t <- mean(sign(AME.main1.s) ==sign(AME.main1.orig))
            if(AME.main1.orig==0) AME1.stab.t <- 1 - AME1.stab.t
            AME1.stab <- c(AME.main1.orig,AME1.stab.t)
            ## AME1.stab <- c(AME.main1.orig,mean(sign(AME.main1.s) ==sign(AME.main1.orig)))
            ## AME2.stab <- c(AME.main2.orig,mean(sign(AME.main2.s) ==sign(AME.main2.orig)))
            stab.prob <- c(ACE.stab[2],AMIE.stab[2],AME1.stab[2],AME2.stab[2],".",".",".")
            estimate <- c(ACE.main,AMIE.main.orig,AME.main1.orig,AME.main2.orig,
                          AMIE.main.orig/ACE.main, AME.main1.orig/ACE.main,
                          AME.main2.orig/ACE.main)
            estimate <- round(estimate,digits=3)
            decom <- as.data.frame(cbind(estimate,as.data.frame(stab.prob)))
            rownames(decom) <- c("Combination Effect","AMIE",AME1.name,AME2.name,
                                 "Prop.AMIE",prop.AME1.name,prop.AME2.name)
        }
        colnames(decom)  <- c("Estimate","select.prob")
        ## decom <- round(decom,digits=3)
    }

    AMIE.main <- output$AMIE.main
    
    cat("\nAMIE:\n")
    print(AMIE.main)
    
    cat("\nBaseline:\n")
    print(baseline)
    
    cat("\nDecomposition:\n")
    print(decom)
    
    output <- list("AMIE.main"=AMIE.main,"baseline"=baseline,"decompose"=decom)
}
