makeallway<-function(X,threshold=0.999999,deletion=TRUE,
                     make.reference=TRUE,
                     sparse.use=FALSE,
		     nway){
    if(!is.vector(X) & !is.matrix(X)){
        warning("X should be a vector or matrix.")
    }
    X2<-NULL
    ## Make indicator variables for all columns.
    if(is.matrix(X)){
        matrix <- 1
        X <- data.frame(X)
        for(i in 1:dim(X)[2]) {
            X2[[i]]<-cbind(sapply(sort(unique(X[,i])),FUN=function(j) 1*(X[,i]==j)))
            colnames(X2[[i]])<-c(paste(names(X[i]), sort(unique(X[,i])),sep="_"))
        }
        n <- ncol(X)
    }
    ## NEW data sets
    ##One-way 
    one.way.data <- NULL
    one.way.data <- as.data.frame(X2)
    if(sparse.use==TRUE){
        one.way.data <- as.matrix(one.way.data)
        one.way.data <- Matrix(one.way.data,sparse=TRUE)
    }
    ## print(head(one.way.data))
    
    
    if(matrix == 1){
        ##Two-way
        if(ncol(X)>=2){
            two.way.data <- NULL
            two.way.name.w <- NULL
            for(j in 1:choose(n,2)){
                formula  <- ~ X2[[combn(n,2)[1,j]]]:X2[[combn(n,2)[2,j]]]
                two.way.data[[j]] <- model.matrix(formula)[,-1]
                if(sparse.use==TRUE){
                    two.way.data[[j]] <- Matrix(two.way.data[[j]],sparse=TRUE)
                }
                two.way.name.w[[j]] <- expand.grid(colnames(X2[[combn(n,2)[1,j]]]),
                                                   colnames(X2[[combn(n,2)[2,j]]]))
                name <- c()
                for(i in 1:nrow(two.way.name.w[[j]])){
                    name[i] <- paste(two.way.name.w[[j]][i,1],
                                     two.way.name.w[[j]][i,2],sep=":")
                }
                colnames(two.way.data[[j]]) <- name
            }
            if(sparse.use==TRUE){
                two.way.data <- do.call("cbind",two.way.data)
            }else{
                two.way.data <- as.data.frame(two.way.data)
            }
            ## two.way.data <- as.matrix(two.way.data)
            ## two.way.data <- Matrix(two.way.data,sparse=TRUE)
        }
        
        ## Three-way
         if(nway>=3){
             three.way.data <- NULL
             three.way.name.w <- NULL
             for(j in 1:choose(n,3)){
                 formula  <- ~ X2[[combn(n,3)[1,j]]]:X2[[combn(n,3)[2,j]]]:
                     X2[[combn(n,3)[3,j]]]
                 three.way.data[[j]] <- model.matrix(formula)[,-1]
                 three.way.name.w[[j]] <- expand.grid(colnames(X2[[combn(n,3)[1,j]]]),
                                                      colnames(X2[[combn(n,3)[2,j]]]),
                                                      colnames(X2[[combn(n,3)[3,j]]]))
                 name <- c()
                 for(i in 1:nrow(three.way.name.w[[j]])){
                     name[i] <- paste(three.way.name.w[[j]][i,1],
                                      three.way.name.w[[j]][i,2],
                                      three.way.name.w[[j]][i,3],
                                      sep=":")
                 }
                 colnames(three.way.data[[j]]) <- name
             }
             three.way.data <- as.data.frame(three.way.data)
             ## three.way.data <- as.matrix(three.way.data)
             ## three.way.data <- Matrix(three.way.data,sparse=TRUE)
         }
        
        
        ## Four-way
        if(nway>=4){
            four.way.data <- NULL
            four.way.name.w <- NULL
            for(j in 1:choose(n,4)){
                formula  <- ~ X2[[combn(n,4)[1,j]]]:X2[[combn(n,4)[2,j]]]:
                    X2[[combn(n,4)[3,j]]]:X2[[combn(n,4)[4,j]]]
                four.way.data[[j]] <- model.matrix(formula)[,-1]
                four.way.name.w[[j]] <- expand.grid(colnames(X2[[combn(n,4)[1,j]]]),
                                                    colnames(X2[[combn(n,4)[2,j]]]),
                                                    colnames(X2[[combn(n,4)[3,j]]]),
                                                    colnames(X2[[combn(n,4)[4,j]]])
                                                    )
                name <- c()
                for(i in 1:nrow(four.way.name.w[[j]])){
                    name[i] <- paste(four.way.name.w[[j]][i,1],
                                     four.way.name.w[[j]][i,2],
                                     four.way.name.w[[j]][i,3],
                                     four.way.name.w[[j]][i,4],
                                     sep=":")
                }
                colnames(four.way.data[[j]]) <- name
            }
            four.way.data <- as.data.frame(four.way.data)
            ## four.way.data <- as.matrix(four.way.data)
            ## four.way.data <- Matrix(four.way.data,sparse=TRUE)
        }
    }
    
    ## if(is.vector(X)){FinalData <- as.data.frame(one.way.data)}
    if(ncol(X)==1){FinalData <- one.way.data}
    if(matrix==1){
        if(sparse.use==TRUE){
            if(ncol(X)>=2){FinalData <- cbind(one.way.data,two.way.data)}
        }else{
            if(ncol(X)>=2){FinalData <- cbind(one.way.data,two.way.data)}
        }
        
        if(nway==3){FinalData <- cbind(one.way.data,two.way.data,three.way.data)}
        if(nway==4){FinalData <- cbind(one.way.data,two.way.data,three.way.data,
                   four.way.data)}
    }
    
    ## FinalData <- as.data.frame(FinalData)
    ## print(head(FinalData))
    
    if(deletion==TRUE){
        ## Discard the columns with no variation or the same variation
        ## with other columns,
        ## while keeping the column names
        No.variation <- colnames(FinalData[,apply(FinalData,2,sd)==0])
        T2 <- FinalData[,apply(FinalData,2,sd)>0]
        T2.sparse <- T2
        
        T2 <- as.matrix(T2)
        ## T2 <- t(T2)
        coldrop      <- sapply(1:ncol(T2),
                               FUN=function(x)
                               which(cor(T2[,x] ,T2)^2 > threshold))
        
        col.names.w  <- sapply(coldrop, FUN=function(x) colnames(T2)[x])
        col.names    <- list(Perfect.Correlated.Columns=col.names.w,
                             no.variation.columns=No.variation)
        
        ## it could be three columns are the same. 
        coldrop.cols <- sapply(coldrop,FUN=function(x) x[1])
        Keep.cols <- unique(coldrop.cols)
        ## if there are more than two, the first one is gotten.
        ## so the second or the third one is the discarded one.
        ## if there are only one, then it is the same.
        
        ## colnames(T2)[coldrop[[Keep.cols]][-1]]
        ## This gives the corresponding dropped coefficients.
        T3 <- T2.sparse[,Keep.cols]
        
        if(make.reference==TRUE){
            Corresponding <- list()
            ## Keep Discarded Data.
            Discarded.cols <- seq(1:ncol(T2))[-Keep.cols]
            
            for(i in Keep.cols){
                if(length(coldrop[[i]][-1])==0){
                    Corresponding[[i]] <- "No Match"
                }
                if(length(coldrop[[i]][-1])!=0){
                    Corresponding[[i]] <- colnames(T2)[coldrop[[i]][-1]]
                }
            }
            for(i in Discarded.cols){
                Corresponding[[i]] <- "Discarded"
            }
            
            C <- max(sapply(Corresponding, FUN=function(x) length(x)))
            
            CorrespondingM <- matrix(NA,nrow=ncol(T2), ncol=C)
            Variation <- colnames(T2)
            rownames(CorrespondingM) <- colnames(T2)
            for(j in 1:C){
                for(i in 1:ncol(T2)){
                    if(length(Corresponding[[i]])>=j){
                        CorrespondingM[i,j] <- Corresponding[[i]][j]
                    }
                }
            }
            
            Reference <- matrix(NA,nrow=ncol(FinalData), ncol=C)
            rownames(Reference) <- colnames(FinalData)
            
            for(i in 1:nrow(Reference)){
                if(rownames(Reference)[i] %in% No.variation){
                    Reference[i,1] <- "No Variation"
                }
                if(rownames(Reference)[i] %in% Variation){
                    Reference[i,]  <- CorrespondingM[rownames(CorrespondingM)==
                                                     rownames(Reference)[i]]
                }
            }
            
            names <- c("Matched Variables", rep("Other matched variables",C-1))
            colnames(Reference) <- names
            Reference <- as.data.frame(Reference)
        }
        else{
            Reference <- NULL}
    }
    if(deletion==FALSE){
        T3 <- FinalData
        Reference <- NULL
    }
    
    return(list(FinalData=T3, reference=Reference))
}
