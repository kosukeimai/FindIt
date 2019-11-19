#' Estimating the AMEs and AMIEs after Regularization with the CausalANOVA.
#' 
#' \code{test.CausalANOVA} estimates the AMEs and AMIEs with confidence
#' intervals after regularization with \code{CausalANOVA}function.
#' 
#' See Details in \code{CausalANOVA}.
#' 
#' @param fit The output from \code{CausalANOAV} function.
#' @param newdata A data frame to use for re-estimating the AMEs and AMIEs with
#' confidence intervals.
#' @param collapse.level A logical indicating whether to collapse insignificant
#' levels within factors as suggested by the \code{CausalANOVA} output users
#' provide.
#' @param diff A logical indicating whether the outcome is the choice between a
#' pair. If \code{diff=TRUE}, \code{pair.id} should specify a pair of
#' comparison. Default is \code{FALSE}.
#' @param pair.id (optional).Unique identifiers for each pair of comparison.
#' This option is used when \code{diff=TRUE}.
#' @param cluster Unique identifies with which cluster standard errors are
#' computed.
#' @return \item{fit}{The output of class \code{CausalANOVA}.}
#' @author Naoki Egami and Kosuke Imai.
#' @seealso \link{CausalANOVA}.
#' @references Egami, Naoki and Kosuke Imai. 2019. Causal Interaction in
#' Factorial Experiments: Application to Conjoint Analysis, Journal of the American Statistical Association.
#' \url{http://imai.fas.harvard.edu/research/files/int.pdf}
#' 
#' Lim, M. and Hastie, T. 2015. Learning interactions via hierarchical
#' group-lasso regularization. Journal of Computational and Graphical
#' Statistics 24, 3, 627--654.
#' 
#' Post, J. B. and Bondell, H. D. 2013. ``Factor selection and structural
#' identification in the interaction anova model.'' Biometrics 69, 1, 70--79.
#' @examples
#' ## ####################################### 
#' ## With Screening and Collapsing
#' ## #######################################
#' data(Carlson)
#' ## Specify the order of each factor
#' Carlson$newRecordF<- factor(Carlson$newRecordF,ordered=TRUE,
#'                             levels=c("YesLC", "YesDis","YesMP",
#'                                      "noLC","noDis","noMP","noBusi"))
#' Carlson$promise <- factor(Carlson$promise,ordered=TRUE,levels=c("jobs","clinic","education"))
#' Carlson$coeth_voting <- factor(Carlson$coeth_voting,ordered=FALSE,levels=c("0","1"))
#' Carlson$relevantdegree <- factor(Carlson$relevantdegree,ordered=FALSE,levels=c("0","1"))
#'
#' ## Sample Splitting
#' train.ind <- sample(unique(Carlson$respcodeS), 272, replace=FALSE)
#' test.ind <- setdiff(unique(Carlson$respcodeS), train.ind)
#' Carlson.train <- Carlson[is.element(Carlson$respcodeS,train.ind), ]
#' Carlson.test <- Carlson[is.element(Carlson$respcodeS,test.ind), ]
#'  
#' #################### AMEs and two-way AMIEs ####################
#' fit.r2 <- CausalANOVA(formula=won ~ newRecordF + promise + coeth_voting + relevantdegree,
#'                       data=Carlson.train, pair.id=Carlson.train$contestresp,diff=TRUE,
#'                       screen=TRUE, collapse=TRUE,
#'                       cluster=Carlson.train$respcodeS, nway=2)
#' summary(fit.r2)
#' 
#' ## refit with test.CausalANOVA
#' fit.r2.new <- test.CausalANOVA(fit.r2, newdata=Carlson.test, diff=TRUE,
#'                                pair.id=Carlson.test$contestresp, cluster=Carlson.test$respcodeS)
#' 
#' summary(fit.r2.new)
#' plot(fit.r2.new)
#' plot(fit.r2.new, type="ConditionalEffect", fac.name=c("newRecordF","coeth_voting"))
#' ConditionalEffect(fit.r2.new, treat.fac="newRecordF", cond.fac="coeth_voting")
#' @export
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
