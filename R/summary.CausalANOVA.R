summary.CausalANOVA <-function(object,verbose=TRUE,...){

    if(class(object)[2]=="stab"){
        fit <- object$fit
        stab.fit <- object$stab.fit
    }else{
        fit <- object
    }
    
    data <- fit$data
    Gorder <- fit$nway
    formula <- fit$formula
    formula.orig <- fit$formula.orig
    cost <- fit$cost
    n.fac <- length(fit$fac.level)
    fac.level <- fit$fac.level
    
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
    for(z in 1:length(fit$coef)){
        range.reg[z] <-  max(fit$coef[[z]]) - min(fit$coef[[z]])
    }
    range <- round(range.reg,digits=3)

    Range <-as.data.frame(range)
    rownames(Range) <- NAME

    ## ################
    ## AME
    ## ################
    AME <- c()
    AME.name <- c()
    for(z in 1:n.fac){
        AME.z <-  fit$coef[[z]] - fit$coef[[z]][length(fit$coef[[z]])]
        AME.z <- round(AME.z,digits=3)
        AME.z <- c("",AME.z)
        AME.z[length(AME.z)] <- "base"
        names(AME.z)[1] <- paste("Fac: ", NAME[z],sep="")
        AME.z <- c(AME.z,"")
        names(AME.z)[length(AME.z)] <- paste(rep(" ",z),collapse="")
        AME <- c(AME,AME.z)
    }
    AME.table <- as.data.frame(AME)
    colnames(AME.table) <- "AME"

    if(class(object)[2]=="stab"){
        Range$stability <- stab.fit$range.stab[,2]
        AME.stab <- stab.fit$AME.stab[,2]

        AME.stab.tab <- c()
        start <- c(1,c(cumsum(fac.level)[-n.fac] + 1))
        end <- cumsum(fac.level)
        for(z in 1:n.fac){
            stab.t <- round(AME.stab[start[z]:end[z]],digits=3)
            stab.t[length(stab.t)] <- "."
            AME.stab.tab <- c(AME.stab.tab,c("", stab.t,""))
        }
        AME.table$stability <- AME.stab.tab

        Range.pr <- Range
        AME.table.pr <- AME.table
        AME.table.pr$select.prob <- AME.table.pr$stability
        AME.table.pr$stability <- NULL
        colnames(Range.pr)[2] <- colnames(AME.table.pr)[2] <- "select.prob"
        Range.pr$select.prob <- round(Range.pr$select.prob,digits=3)
        ## AME.table.pr$select.prob[is.numeric(AME.table.pr$select.prob)] <- 
        ##     round(AME.table.pr$select.prob[is.numeric(AME.table.pr$select.prob)],digits=3)
    }else{
        Range.pr <- Range
        AME.table.pr <- AME.table
    }

    ## ## #################
    ## ## Compute each AMIEs by specifying the baseline
    ## ## #################

    ## AMIE <- c()
    ## AMIE.na <- c()
    ## for(z in 1:length(fit$coef)){
    ##     AMIE.z <-  fit$coef[[z]] - fit$coef[[z]][length(fit$coef[[z]])]
    ##     AMIE.z.na <- AMIE.z
    ##     AMIE.z.na[length(AMIE.z.na)] <- NA
    ##     AMIE.na <- c(AMIE.na,AMIE.z.na)
    ##     AMIE <- c(AMIE,AMIE.z)
    ## }
    ## AMIE <- round(AMIE,digits=3)
    ## AMIE.na <- round(AMIE.na,digits=3)
    ## AMIE.table <- as.data.frame(AMIE)
    ## colnames(AMIE.table) <- "AME/AMIE"

   

    if(verbose==TRUE){
        cat("\nCall:\n")
        cat(" \nModel:\n ")
        print(formula.orig)
        cat(paste("\nInteraction: ",Gorder,"-way Interaction\n",sep=""))
        
        cat(paste("\nCost parameter : ",cost,"\n",sep=""))
        
        cat("\nRange:\n")
        print(Range.pr)
        
        
        cat("\nAME:\n")
        print(AME.table.pr)
    }

    output <- list("range"=range,"range.name"=NAME,"AME"=AME.table)
}
