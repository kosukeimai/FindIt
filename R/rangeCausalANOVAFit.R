rangeCausalANOVAFit <-function(object,verbose=TRUE,...){
    
    data <- object$data
    Gorder <- object$nway
    formula <- object$formula
    formula.orig <- object$formula.orig
    cost <- object$cost
    n.fac <- length(object$fac.level)
    
    ## ################## 
    ## Name
    ## ################# 
    data.name <- colnames(model.frame(formula= ~ .*., data=data)[1,-1])
    ## Gorder =2 
    combTwo <- combn(length(data.name),2)
    NameInt <- c()
    for(i in 1:ncol(combTwo)){
        NameInt[i] <- paste(data.name[combTwo[1,i]],data.name[combTwo[2,i]],sep=":")
    }
    ## Gorder = 3
    if(Gorder==3){
        combThree <- combn(length(data.name),3)
        threename <- c()
        for(i in 1:ncol(combThree)){
            threename[i] <- paste(data.name[combThree[1,i]],
                                  data.name[combThree[2,i]],
                                  data.name[combThree[2,i]],
                                  sep=":")
        }
        NameInt <- c(NameInt,threename)
    }
    NAME <- c(data.name,NameInt)
    
    ## ################# 
    ## Compute Range
    ## ################# 
    range.reg <- c()
    for(z in 1:length(object$coef)){
        range.reg[z] <-  max(object$coef[[z]]) - min(object$coef[[z]])
    }
    ## range <- round(range.reg,digits=3)
    range <- range.reg

    Range <-as.data.frame(range)
    rownames(Range) <- NAME

    ## ################
    ## AMIE
    ## ################
    AMIE <- c()
    for(i in 1:(length(object$coef)-n.fac)){
        z <- i + n.fac
        AMIE.z <-  object$coef[[z]] - object$coef[[z]][length(object$coef[[z]])]
        AMIE.z[length(AMIE.z)] <- 0
        AMIE <- c(AMIE,AMIE.z)       
    }
    AMIE.table <- as.data.frame(AMIE)        

    ## ## #################
    ## ## Compute each AMIEs by specifying the baseline
    ## ## #################

    AME <- c()
    for(z in 1:n.fac){
        AME.z <-  object$coef[[z]] - object$coef[[z]][length(object$coef[[z]])]
        AME.z[length(AME.z)] <- 0
        AME <- c(AME,AME.z)
    }
    ## AME <- round(AME,digits=3)
    ## AME.na <- round(AME.na,digits=3)
    AME.table <- as.data.frame(AME)
    ## colnames(AME.table) <- "AME/AMIE"

    ## if(verbose==TRUE){
    ##     cat("\nCall:\n")
    ##     cat(" \nModel:\n ")
    ##     print(formula.orig)
    ##     cat(paste("\nInteraction: ",Gorder,"-way Interaction\n",sep=""))
        
    ##     cat(paste("\nCost parameter : ",cost,"\n",sep=""))
        
    ##     cat("\nRange:\n")
    ##     print(Range)
        
        
    ##     cat("\nAME:\n")
    ##     print(AME.table)
    ## }

    output <- list("range"=range,"range.name"=NAME,
                   "AME"=AME,"AMIE"=AMIE)
}
