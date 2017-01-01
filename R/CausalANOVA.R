CausalANOVA <- function(formula,data,cost,pair.id=NULL,
                        nway=2,diff=TRUE,
                        select.prob=FALSE,boot=100,block.id=NULL,
                        seed=1234,
                        eps=1e-5,
                        fac.level=NULL,ord.fac=NULL,
                        verbose=TRUE){    
    
    if(missing(fac.level)) fac.level <- NULL
    if(missing(ord.fac)) ord.fac <- NULL
    
    fit <- CausalANOVAFit(formula=formula,data=data,cost=cost,pair.id=pair.id,
                          nway=nway,diff=diff,eps=eps,
                          fac.level=fac.level,ord.fac=ord.fac,
                          verbose=verbose)
    if(select.prob==TRUE){
        stab.fit <- stab.CausalANOVA(fit,block.id=block.id,boot=boot)
        output <- list("fit"=fit, "stab.fit"=stab.fit)
        class(output) <- c("CausalANOVA","stab","list")
    }else{
        output <- fit
        class(output) <- c("CausalANOVA","list")
    }        
    return(output)
}
