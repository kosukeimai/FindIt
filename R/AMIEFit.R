AMIEFit <-function(object,fac.name,level.name,base.name,verbose=TRUE,...){
    
    if(class(object)[2]=="stab"){
        fit <- object$fit
    }else{
        fit <- object
    }
    
    data <- fit$data
    Gorder <- fit$nway
    formula <- fit$formula
    cost <- fit$cost
    n.fac <- length(fit$fac.level)

    ## AME
    AME <- c()
    for(z in 1:n.fac){
        AME.z <-  fit$coef[[z]] - fit$coef[[z]][length(fit$coef[[z]])]
        AME <- c(AME,AME.z)
    }
    AME.table <- as.data.frame(AME)
    colnames(AME.table) <- "AME"

    ## #################
    ## Compute each AMIEs by specifying the baseline
    ## #################
    
    AMIE <- c()
    for(i in 1:(length(fit$coef)-n.fac)){
        z <- i + n.fac
        AMIE.z <-  fit$coef[[z]] - fit$coef[[z]][length(fit$coef[[z]])]
        AMIE <- c(AMIE,AMIE.z)
    }
  
    data.x <- model.frame(formula,data=data)[,-1]
    fac.ind <- which(colnames(data.x) %in% fac.name)

    ## AMIE
    ind.fac1 <- as.numeric(regexpr(fac.name[1],names(AMIE))>0)
    ind.fac2 <- as.numeric(regexpr(fac.name[2],names(AMIE))>0)
    ## AME
    ind.AME.fac1 <- as.numeric(regexpr(fac.name[1],names(AME))>0)
    ind.AME.fac2 <- as.numeric(regexpr(fac.name[2],names(AME))>0)
    
    if(missing(base.name)){              
        baseline <- c(tail(levels(data.x[,fac.ind[1]]),1),tail(levels(data.x[,fac.ind[2]]),1))
        AMIE.base <- 0
        AME.base1 <- 0
        AME.base2 <- 0
        names(baseline) <- colnames(data.x)[fac.ind]
    }else{
        ## AMIE
        ind.base1 <- as.numeric(regexpr(base.name[1],names(AMIE))>0)
        ind.base2 <- as.numeric(regexpr(base.name[2],names(AMIE))>0)
        ind.base <- ind.fac1*ind.fac2*ind.base1*ind.base2 
        AMIE.base <- AMIE[ind.base==1]

        ## AME
        ind.AME.base1t <- as.numeric(regexpr(base.name[1],names(AME))>0)
        ind.AME.base2t <- as.numeric(regexpr(base.name[2],names(AME))>0)
        ind.AME.base1 <- ind.AME.fac1*ind.AME.base1t
        ind.AME.base2 <- ind.AME.fac2*ind.AME.base2t
        AME.base1 <- AME[ind.AME.base1==1]
        AME.base2 <- AME[ind.AME.base2==1]
        
        baseline <- base.name
        names(baseline) <- fac.name
    }

    ## Name match
    ind.lev1 <- as.numeric(regexpr(level.name[1],names(AMIE))>0)
    ind.lev2 <- as.numeric(regexpr(level.name[2],names(AMIE))>0)
    ind <- ind.fac1*ind.fac2*ind.lev1*ind.lev2
    AMIE.main <- AMIE[ind==1] - AMIE.base

    ind.AME.lev1 <- as.numeric(regexpr(level.name[1],names(AME))>0)
    ind.AME.lev2 <- as.numeric(regexpr(level.name[2],names(AME))>0)
    ind.AME1 <- ind.AME.fac1*ind.AME.lev1
    ind.AME2 <- ind.AME.fac2*ind.AME.lev2
    AME.main1 <- AME[ind.AME1==1] - AME.base1
    AME.main2 <- AME[ind.AME2==1] - AME.base2
    
    ACE <- AMIE.main + AME.main1 + AME.main2
    
    AME1.name <- paste("AME (",fac.name[1],")",sep="")
    AME2.name <- paste("AME (",fac.name[2],")",sep="")
    prop.AME1.name <- paste("Prop.",AME1.name,sep="")
    prop.AME2.name <- paste("Prop.",AME2.name,sep="")

    if(any(is.element(baseline,level.name)==TRUE)){
        if(AME.main1==0){
            decom <- as.data.frame(c(ACE,AMIE.main,AME.main2,AMIE.main/ACE,
                                     AME.main2/ACE))
            rownames(decom) <- c("Conditional Effect","AMIE",AME2.name,
                                 "Prop.AMIE",prop.AME2.name)
        }else{
            decom <- as.data.frame(c(ACE,AMIE.main,AME.main1,AMIE.main/ACE,
                                     AME.main1/ACE))
            rownames(decom) <- c("Conditional Effect","AMIE",AME1.name,
                                 "Prop.AMIE",prop.AME1.name)
        }
    }else{
        decom <- as.data.frame(c(ACE,AMIE.main,AME.main1,AME.main2,AMIE.main/ACE,
                                 AME.main1/ACE,AME.main2/ACE))
        rownames(decom) <- c("Combination Effect","AMIE",AME1.name,AME2.name,
                             "Prop.AMIE",prop.AME1.name,prop.AME2.name)
    }
    colnames(decom)  <- "Estimate"
    decom <- round(decom,digits=3)
    
    ## cat("\nAMIE:\n")
    ## print(AMIE.main)

    ## cat("\nBasline:\n")
    ## print(baseline)
    
    ## cat("\nDecomposition:\n")
    ## print(decom)
   
    output <- list("AMIE.main"=AMIE.main,"baseline"=baseline,"decompose"=decom)
}
