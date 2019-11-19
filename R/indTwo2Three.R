indTwo2Three <- function(indTwo, n.fac){
    Ad <- matrix(0, ncol=n.fac, nrow=n.fac)
    for(j in 1:ncol(indTwo)){
        Ad[indTwo[1,j],indTwo[2,j]] <- 1
    }
    Ad <- Ad + t(Ad)
    indThree0 <- matrix(triangles(graph_from_adjacency_matrix(Ad)), nrow=3)
    indThree <- apply(indThree0, 2, function(x) x[order(x)])
    if(length(indThree)==0) indThree <- NULL
    return(indThree)
}
