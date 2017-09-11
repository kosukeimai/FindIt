rangeCausalANOVAFit <-function(object,verbose=TRUE,...){
    
    data <- object$data
    Gorder <- object$Gorder
    formula <- object$formula
    formula.orig <- object$formula.orig
    cost <- object$cost
    n.fac <- length(object$fac.level)
    terms.f <- terms(formula, data=data)
    order.f <- attr(terms.f, "order")
    Name <- attr(terms.f, "term.labels")
    
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
    rownames(Range) <- Name

    ## ################
    ## AMIE
    ## ################
    if(Gorder>=2){
        AMIE <- c()
        for(i in 1:(length(object$coef) - n.fac)){
            z <- i + n.fac
            AMIE.z <-  object$coef[[z]] - object$coef[[z]][length(object$coef[[z]])]
            AMIE.z[length(AMIE.z)] <- 0
            AMIE <- c(AMIE,AMIE.z)       
        }
        AMIE.table <- as.data.frame(AMIE)
    }else{
        AMIE <- NULL
    }

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
    

    output <- list("range"=range,"range.name"=Name,
                   "AME"=AME,"AMIE"=AMIE)
}
