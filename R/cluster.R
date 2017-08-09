cluster_se_glm <- function(model, cluster){
    
    #  Drop unused cluster indicators, if cluster var is a factor
    if (class(cluster) == "factor") {
        cluster <- droplevels(cluster)
    }
    
    if (nrow(model.matrix(model)) != length(cluster)) {
        stop("check your data: cluster variable has different N than model - you may have observations with missing data") 
    }
    
    M <- length(unique(cluster))
    N <- length(cluster)           
    K <- model$rank
    
    ## if(M<50) {
    ##     warning("Fewer than 50 clusters, variances may be unreliable (could try block bootstrap instead).")
    ## }
   
    dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
    uj  <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum));
    rcse.cov <- dfc * sandwich(model, meat. = crossprod(uj)/N)
    return(rcse.cov)
}
