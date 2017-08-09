
test.CausalANOVA <- function(fit, newdata, collapse.level=TRUE, diff=FALSE, pair.id=NULL,cluster=NULL){
    ## I need to install newdata as well as new information about pair.id
    
    eps <- fit$eps
    formula <- fit$formula
    main.formula <- fit$main.formula
    data <- model.frame(main.formula, fit$data)
    terms.f <- terms(formula, data=data)
    order.f <- attr(terms.f, "order")
    var.name <- attr(terms.f, "term.labels")[order.f==1]
    fac.level <- fit$fac.level
    n.fac <- length(fac.level)
    data <- fit$data
    family <- fit$family
    
    data.main0 <- model.frame(main.formula, data=newdata)
    
    if(collapse.level==TRUE){
        collapse.u <- collapse <- Collapsing(fit)$collapse.level
        drop.fac <- unlist(lapply(collapse, function(x) all(x==1)))
        
        ## Recover Collapsed Factor
        for(z in which(drop.fac)){
            n1 <- ceiling(length(collapse.u[[z]])/2)
            n2 <- length(collapse.u[[z]]) - n1
            collapse.u[[z]] <- c(rep(1, n1), rep(2,n2))
        }
        
        ## New Level Names
        collapse.u2 <- list()
        for(i in 1:n.fac){
            original.level <- levels(data[,(i+1)])
            new.name <- c()
            for(j in 1:length(unique(collapse.u[[i]]))){
                new.name[j] <- paste(original.level[collapse.u[[i]]==j],collapse="/")
            }
            collapse.u2[[i]] <- new.name[collapse.u[[i]]]
        }
        
        ## ############################ 
        ## Reconstruct Formula
        ## ############################ 
        
        fac.table <- attr(terms.f, "factors")[-1,]
        d.table <- fac.table[drop.fac==TRUE,]

        if(sum(drop.fac)==0) drop.ind <- 0
        if(sum(drop.fac)==1) drop.ind <- (d.table > 0)*(order.f>1)
        if(sum(drop.fac)>1){
            drop.ind <- (apply(d.table,2,sum) > 0)*(order.f>1)
        }
        
        term.labels <- attr(terms.f, "term.labels")
        
        if(all(drop.fac==FALSE) | all(drop.ind==0)==TRUE){
            new.formula <- formula
        }else{
            new.formula <-
                as.formula(paste(as.character(formula)[2], "~", paste(term.labels[drop.ind==0], collapse="+"),sep=""))
        }
    }else{
        new.formula <- formula
    }
    ## ############################
    ## Prepare New Data 
    ## ############################
    data.main <- model.frame(new.formula, data=data.main0)
    data.x <- model.frame(new.formula, data=data.main0)[,-1]
    terms.f.new <- terms(new.formula, data=data.main0)
    nway.new <- max(attr(terms.f.new, "order"))
    ord.fac <- rep(FALSE, n.fac)

    ## Releveling New Data
    if(collapse.level==TRUE){
        for(i in 1:n.fac){
            levels(data.main[,(i+1)]) <- collapse.u2[[i]]
            data.main[,(i+1)] <- factor(data.main[,(i+1)],ordered=FALSE)
        }
    }
    
    fit.new <- CausalANOVAFit(formula = new.formula, internal.int=FALSE,
                              pair.id=pair.id, diff=diff, 
                              data=data.main, cost=1, cluster=cluster,
                              family=family, nway=nway.new,
                              eps=eps, screen=FALSE, ord.fac=NULL,
                              collapse=FALSE, verbose=FALSE)
    
    coefs <- fit.new$coefs
    vcov  <- fit.new$vcov
    fit.new$inference <- TRUE
    
    ## Index
    levelIndex <- CreatelevelIndex(fac.level=fit.new$fac.level, ord.fac=fit.new$ord.fac, Gorder=fit.new$Gorder,
                                   indTwo=fit.new$indTwo, indThree=fit.new$indThree)
    use.ind <-  (levelIndex$plus==1)*(levelIndex$dif==0)
    Index.use <- levelIndex[use.ind==1,]
    Index.use$start <- c(1,cumsum(Index.use$length)[-nrow(Index.use)]+1)
    Index.use$end <-  cumsum(Index.use$length)
    
    order.f <- attr(terms.f.new, "order")

    result.tab <- list()
    for(z in 1:length(coefs)){
        seq.u <- Index.use$start[z]:Index.use$end[z]
        effect <- coefs[[z]] - coefs[[z]][length(coefs[[z]])]
        vcov.p <- c()
        for(i in 1:(length(seq.u)-1)){
            vcov.p[i] <- VarEffect(ind1=seq.u[i],ind2=seq.u[length(seq.u)],vcov.full=vcov)
        }
        sd <- c(sqrt(vcov.p),0)
        table <- as.data.frame(cbind(effect, sd, effect - 1.96*sd, effect + 1.96*sd))
        if(Index.use$order[z]==1){
            colnames(table) <- c("AME", "sd", "2.5% CI", "97.5% CI")
        }else{
            colnames(table) <- c("AMIE", "sd", "2.5% CI", "97.5% CI")
        }
        result.tab[[z]] <- table
    }

    fit.new[[length(fit.new)+1]] <- result.tab
    names(fit.new)[length(fit.new)] <- "CI.table"
    ## output <- c(fit.new, "result.tab"=result.tab)

    output <- fit.new
    class(output) <- c("CausalANOVA","list")
    return(output)
}
