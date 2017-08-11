#' Estimating the AMEs and AMIEs with the CausalANOVA.
#' 
#' \code{CausalANOVA} estimates coefficients of the specified ANOVA with
#' regularization. By taking differences in coefficients, the function recovers
#' the AMEs and AMIEs.
#' 
#' Regularization: \code{screen} and \code{collapse}.
#' 
#' Users can implement regularization in order to reduces false discovery rate
#' and facilitates interpretation. This is particularly useful when analyzing
#' factorial experiments with a large number of factors, each having many
#' levels.  \itemize{ \item When \code{screen=TRUE}, the function selects
#' significant factor interactions with \code{glinternet} (Lim and Hastie 2015)
#' before estimating the AMEs and AMIEs. This option is recommended when there
#' are many factors, e.g., more than 6 factors. Alternatively, users can
#' pre-specify interactions of interest using \code{int2.formula} and
#' \code{int3.formula}. \item When \code{collapse=TRUE}, the function collapses
#' insignificant levels within each factor by GashANOVA (Post and Bondell 2013)
#' before estimating the AMEs and AMIEs. This option is recommended when there
#' are many levels within some factors, e.g., more than 6 levels. }
#' 
#' Inference after Regularization: \itemize{ \item When \code{screen=TRUE} or
#' \code{collapse=TRUE}, in order to make valid inference after regularization,
#' we recommend to use \code{test.CausalANOVA} function. It takes the output
#' from \code{CausalANOVA} function and estimate the AMEs and AMIEs with
#' \code{newdata} and provide confidence intervals. Ideally, users should split
#' samples into two; use a half for regularization with \code{CausalANOVA}
#' function and use the other half for inference with \code{test.CausalANOVA}.
#' \item If users do not need regularization, specify \code{screen=FALSE} and
#' \code{collapse=FALSE}. The function estiamtes the AMEs and AMIEs and compute
#' confidence intervals with the full sample. }
#' 
#' Suggested Workflow: (See Examples below as well) \enumerate{ \item Specify
#' the order of levels within each factor using \code{levels()}.  When
#' \code{collapse=TRUE}, the function places penalties on the differences
#' between adjacent levels when levels are ordered, it is crucial to specify
#' the order of levels within each factor carefully. \item Run
#' \code{CausalANOVA}. \enumerate{ \item Specify \code{formula} to indicate
#' outcomes and treatment variables and \code{nway} to indicate the order of
#' interactions. \item Specify \code{diff=TRUE} and \code{pair.id} if the
#' outcome is the choice between a pair. \item Specify \code{screen}.
#' \code{screen=TRUE} to implement data-driven selection of factor
#' interactions. \code{screen=FALSE} to specify interactions through
#' \code{int2.formula} and \code{int3.formula} by hand. \item Specify
#' \code{collapse}. \code{collapse=TRUE} to implement data-driven collapsing of
#' insignificant levels. \code{collapse=FALSE} to use the original nubmber of
#' levels.}
#' 
#' \item Run \code{test.CausalANOVA} when \code{select=TRUE} or
#' \code{collapse=TRUE}.  \item Run \code{summary} and \code{plot} to explore
#' the AMEs and AMIEs. \item Estimate conditional effects using
#' \code{ConditionalEffect} function and visualize them using \code{plot}
#' function. }
#' 
#' @aliases CausalANOVA summary.CausalANOVA plot.CausalANOVA
#' @param formula A formula that specifies outcome and treatment variables.
#' @param int2.formula (optional). A formula that specifies two-way
#' interactions.
#' @param int3.formula (optional). A formula that specifies three-way
#' interactions.
#' @param data An optional data frame, list or environment (or object coercible
#' by 'as.data.frame' to a data frame) containing the variables in the model.
#' If not found in 'data', the variables are taken from 'environment(formula)',
#' typically the environment from which 'CausalANOVA' is called.
#' @param nway With \code{nway=1}, the function estimates the Average Marginal
#' Effects (AMEs) only. With \code{nway=2}, the function estimates the AMEs and
#' the two-way Average Marginal Interaction Effects (AMIEs). With
#' \code{nway=3}, the function estiamtes the AMEs, the two-way and three-way
#' AMIEs. Default is 1.
#' @param diff A logical indicating whether the outcome is the choice between a
#' pair.  If \code{diff=TRUE}, \code{pair.id} should specify a pair of
#' comparison. Default is \code{FALSE}.
#' @param pair.id (optional).Unique identifiers for each pair of comparison.
#' This option is used when \code{diff=TRUE}.
#' @param screen A logical indicating whether select significant factor
#' interactions with \code{glinternet}.  When users specify interactions using
#' \code{int2.formula} or \code{int3.formula}, this option is ignored.
#' \code{screen} should be used only when users want data-driven selection of
#' factor-interactions. With \code{screen.type}, users can specify how to
#' screen factor interactions. We recommend to use this option when the number
#' of factors is large, e.g., more than 6. Default is \code{FALSE}.
#' @param screen.type Type for screening factor interactions. (1)
#' \code{"fixed"} select the fixed number (specified by \code{screen.num.int})
#' of factor interactions. (2) \code{"cv.min"} selects factor-interactions with
#' the tuning parameter giving the minimum cross-validation error. (3)
#' \code{"cv.1Std"} selects factor-interactions with the tuning parameter
#' giving a cross-validation error that is within 1 standard deviation of the
#' minimum cv error.
#' @param screen.num.int (optional).The number of factor interactions to
#' select. This option is used when and \code{screen=TRUE} and
#' \code{screen.type="fixed"}. Default is 3.
#' @param collapse A logical indicating whether to collapse insignificant
#' levels within factors.  With \code{collapse.type}, users can specify how to
#' collapse levels within factors. We recommend to use this option when the
#' number of levels is large, e.g., more than 6. Default is \code{FALSE}.
#' @param collapse.type Type for collapsing levels within factors. (1)
#' \code{"fixed"} collapses levels with the fixed cost parameter (specified by
#' \code{collapse.cost}). (2) \code{"cv.min"} collapses levels with the cost
#' paramter giving the minimum cross-validation error. This option might take
#' time. (3) \code{"cv.1Std"} collapses with the cost parameter giving a
#' cross-validation error that is within 1 standard deviation of the minimum cv
#' error. This option might take time.
#' @param collapse.cost (optional).A cost parameter ranging from 0 to 1. 1
#' corresponds to no collapsing. The closer to 0, the stronger regularization.
#' Default is 0.3.
#' @param family A family of outcome varialbes. \code{"gaussian"} when
#' continuous outcomes \code{"binomial"} when binary outcomes.  Default is
#' \code{"binomial"}.
#' @param cluster Unique identifies with which cluster standard errors are
#' computed.
#' @param maxIter The number of maximum iteration for \code{glinternet}.
#' @param eps A tolerance parameter in the internal optimization algorithm.
#' @param fac.level (optional). A vector containing the number of levels in
#' each factor. The order of \code{fac.level} should match to the order of
#' columns in the data. For example, when the first and second columns of the
#' design matrix is "Education" and "Race", the first and second element of
#' \code{fac.level} should be the number of levels in "Education" and "Race",
#' respectively.
#' @param ord.fac (optional). Logical vectors indicating whether each factor
#' has ordered (\code{TRUE}) or unordered (\code{FALSE}) levels.  When levels
#' are ordered, the function uses the order given by function \code{levels()}.
#' If levels are ordered, the function places penalties on the differences
#' between adjacent levels.  If levels are unordered, the function places
#' penalties on the differences based on every pairwise comparison.
#' @param select.prob (optional). A logical indicating whether selection
#' probabilities are computed. This option might take time.
#' @param boot The number of bootstrap replicates for \code{select.prob}.
#' Default is 50.
#' @param seed Seed for bootstrap.
#' @param verbose Whether it prints the value of a cost parameter used.
#' @return \item{intercept}{An intercept of the estimated ANOVA model.If
#' \code{diff=TRUE}, this should be close to 0.5.} \item{formula}{The
#' \code{formula} used in the function.} \item{coefs}{A named vector of
#' coefficients of the estimated ANOVA model.} \item{vcov}{The
#' variance-covariance matrix for \code{coefs}. Only when \code{select=FALSE}
#' and \code{collapse=FALSE}.} \item{CI.table}{The summary of AMEs and AMIEs
#' with confidence intervals. Only when \code{select=FALSE} and
#' \code{collapse=FALSE}.} \item{AME}{The estimated AMEs with the grand-mean as
#' baselines.} \item{AMIE2}{The estimated two-way AMIEs with the grand-mean as
#' baselines.} \item{AMIE3}{The estimated three-way AMIEs with the grand-mean
#' as baselines.} \item{...}{arguments passed to the function or arguments only
#' for the internal use.}
#' @author Naoki Egami and Kosuke Imai.
#' @seealso \link{cv.CausalANOVA}
#' @references Egami, Naoki and Kosuke Imai. 2016+. Causal Interaction in
#' Factorial Experiments: Application to Conjoint Analysis. Working paper.
#' \url{http://imai.princeton.edu/research/files/int.pdf}
#' 
#' Lim, M. and Hastie, T. 2015. Learning interactions via hierarchical
#' group-lasso regularization. Journal of Computational and Graphical
#' Statistics 24, 3, 627--654.
#' 
#' Post, J. B. and Bondell, H. D. 2013. Factor selection and structural
#' identification in the interaction anova model. Biometrics 69, 1, 70--79.
#' @examples
#' 
#' data(Carlson)
#' ## Specify the order of each factor
#' Carlson$newRecordF<- factor(Carlson$newRecordF,ordered=TRUE,
#'                          levels=c("YesLC", "YesDis","YesMP",
#'                              "noLC","noDis","noMP","noBusi"))
#' Carlson$promise <- factor(Carlson$promise,ordered=TRUE,levels=c("jobs","clinic","education"))
#' Carlson$coeth_voting <- factor(Carlson$coeth_voting,ordered=FALSE,levels=c("0","1"))
#' Carlson$relevantdegree <- factor(Carlson$relevantdegree,ordered=FALSE,levels=c("0","1"))
#' 
#' ## ####################################### 
#' ## Without Screening and Collapsing
#' ## ####################################### 
#' ## only AMEs 
#' fit1 <- CausalANOVA(formula=won ~ newRecordF + promise + coeth_voting + relevantdegree,
#'                     data=Carlson, pair.id=Carlson$contestresp, diff=TRUE,
#' 		    cluster=Carlson$respcodeS, nway=1)
#' summary(fit1)
#' # plot(fit1)
#' 
#' ## AMEs and two-way AMIEs 
#' fit2 <- CausalANOVA(formula=won ~ newRecordF + promise + coeth_voting + relevantdegree,
#'                     int2.formula = ~ newRecordF:coeth_voting,
#' 		    data=Carlson, pair.id=Carlson$contestresp,diff=TRUE,
#' 		    cluster=Carlson$respcodeS, nway=2)
#' summary(fit2)
#' plot(fit2)
#' plot(fit2, type="ConditionalEffect", fac.name=c("newRecordF","coeth_voting"))
#' ConditionalEffect(fit2, treat.fac="newRecordF", cond.fac="coeth_voting")
#' 
#' \dontrun{
#' ## AMEs and two-way and three-way AMIEs
#' ## Note: All pairs within thee-way interactions should show up in int2.formula (Strong Hierarchy).
#' fit3 <- CausalANOVA(formula=won ~ newRecordF + promise + coeth_voting + relevantdegree,
#'                     int2.formula = ~ newRecordF:promise + newRecordF:coeth_voting
#' 		                         + promise:coeth_voting,
#' 		    int3.formula = ~ newRecordF:promise:coeth_voting,
#' 		    data=Carlson, pair.id=Carlson$contestresp,diff=TRUE,
#' 		    cluster=Carlson$respcodeS, nway=3)
#' summary(fit3)
#' plot(fit3)
#' plot(fit3, type="AMIE", fac.name=c("newRecordF","promise", "coeth_voting"),space=25,adj.p=2.2)
#' ConditionalEffect(fit3, treat.fac="newRecordF", cond.fac="coeth_voting")
#' }
#' 
#' ## ####################################### 
#' ## With Screening and Collapsing
#' ## #######################################
#' ## Sample Splitting
#' train.ind <- sample(unique(Carlson$respcodeS), 272, replace=FALSE)
#' test.ind <- setdiff(unique(Carlson$respcodeS), train.ind)
#' Carlson.train <- Carlson[is.element(Carlson$respcodeS,train.ind), ]
#' Carlson.test <- Carlson[is.element(Carlson$respcodeS,test.ind), ]
#' 
#' ## only AMEs (Note: when nway=1, there is no factor interaction, so screen=FALSE).
#' fit.r1 <- CausalANOVA(formula=won ~ newRecordF + promise + coeth_voting + relevantdegree,
#'                       data=Carlson.train, pair.id=Carlson.train$contestresp, diff=TRUE,
#' 		      collapse=TRUE,
#' 		      cluster=Carlson.train$respcodeS, nway=1)
#' summary(fit.r1)
#' 
#' ## refit with test.CausalANOVA
#' fit.r1.new <- test.CausalANOVA(fit.r1, newdata=Carlson.test, diff=TRUE,
#'                                pair.id=Carlson.test$contestresp, cluster=Carlson.test$respcodeS)
#' summary(fit.r1.new)
#' plot(fit.r1.new)
#' 
#' ## AMEs and two-way AMIEs 
#' fit.r2 <- CausalANOVA(formula=won ~ newRecordF + promise + coeth_voting + relevantdegree,
#'                       data=Carlson.train, pair.id=Carlson.train$contestresp,diff=TRUE,
#' 		      screen=TRUE, collapse=TRUE,
#' 		      cluster=Carlson.train$respcodeS, nway=2)
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
#' 
#' \dontrun{
#' ## AMEs and two-way and three-way AMIEs
#' fit.r3 <- CausalANOVA(formula=won ~ newRecordF + promise + coeth_voting + relevantdegree,
#'                       data=Carlson, pair.id=Carlson$contestresp,diff=TRUE,
#' 		      screen=TRUE, collapse=TRUE,
#' 		      cluster=Carlson$respcodeS, nway=3)
#' summary(fit.r3)
#' 
#' ## refit with test.CausalANOVA
#' fit.r3.new <- test.CausalANOVA(fit.r3, newdata=Carlson.test, diff=TRUE,
#'                                pair.id=Carlson.test$contestresp, cluster=Carlson.test$respcodeS)
#' 
#' summary(fit.r3.new)
#' ConditionalEffect(fit.r3.new, treat.fac="newRecordF", cond.fac="coeth_voting")
#' }
#' 
CausalANOVA <- function(formula, int2.formula=NULL, int3.formula=NULL,
                        data, nway=1,
                        pair.id=NULL, diff=FALSE,
                        screen=FALSE, screen.type="fixed", screen.num.int=3,
                        collapse=FALSE, collapse.type="fixed", collapse.cost=0.3,
                        family="binomial", cluster=NULL,
                        maxIter=50, eps=1e-5,
                        fac.level=NULL, ord.fac=NULL,
                        select.prob=FALSE, boot=100, seed=1234,
                        verbose=TRUE){    

    ## House Keeping
    if(collapse==FALSE){
        cost <- 1
    }else if(collapse.type=="fixed"){
        cost <- collapse.cost
        if(cost > 1 | cost < 0 ){
            stop("Specify 'collapse.cost' between 0 and 1")
        }
    }
    if(missing(fac.level)) fac.level <- NULL
    if(missing(ord.fac)) ord.fac <- NULL        
    if(diff==TRUE & is.null(pair.id)==TRUE){
        stop("When 'diff=TRUE', specify 'pair.id'.")
    }
    if((nway %in% c(1, 2,3))==FALSE){
        stop("'nway' should be 1, 2 or 3.")
    }
    if(screen==TRUE & nway==1){
        warning("When 'nway=1', no screening is needed ('screen=FALSE').")
        screen <- FALSE
    }
    if((family %in% c("gaussian","binomial"))==FALSE){
        stop("'family' should be 'gaussian' or 'binomial'.")
    }
    y <- model.frame(formula,data=data)[,1]
    if(family=="gaussian" & is.numeric(y)==FALSE){
        stop("When 'family=gaussian', outcomes should be 'numeric'.") 
    }
    if(family=="binomial" & all(y %in% c(0,1))==FALSE){
        stop("When 'family=binomial', outcomes should be 0 or 1.") 
    }
    rm(y)
    if((screen.type %in% c("fixed","cv.min","cv.1Std"))==FALSE){
        stop("'screen.type' should be 'fixed', 'cv.min' or 'cv.1Std'.")
    }    
    if(is.null(int3.formula)==FALSE & nway!=3){
        warning("When 'int3.formula' is specified, 'nway=3'.")
        nway <- 3
    }
    if(is.null(int2.formula)==TRUE & is.null(int3.formula)==FALSE){
        stop("Cannot specify 'int3.formula' without specifying 'int2.formula'.")
    }
    if(is.null(int2.formula)==FALSE & screen==TRUE){
        warning("When 'int2.formula' is specified, 'screen' should be FALSE.")
        screen <- FALSE
    }
    if(nway==3 & screen==TRUE & screen.type=="fixed" & screen.num.int < 3){
        stop("'screen.num.int >=3' when 'nway=3'.")
    }    
    if(collapse.type=="cv.min" | collapse.type=="cv.1Std"){
        cat("\nRunning Cross Validation might take time...\n")
        cv <- cv.CausalANOVA(formula = formula,
                             int2.formula=int2.formula, int3.formula=int3.formula,
                             data=data,
                             cv.collapse.cost=c(0.1,0.3,0.7),
                             nfolds=5,
                             pair.id=pair.id, diff=diff,
                             nway=nway, family=family,
                             seed=seed, cluster=cluster,
                             screen=screen, screen.type=screen.type,
                             screen.num.int=screen.num.int,
                             maxIter=maxIter, eps=eps,                                                         
                             fac.level=fac.level,ord.fac=ord.fac, verbose=verbose)                            
        if(collapse.type=="cv.min") cost <- cv$cv.min
        if(collapse.type=="cv.1Std") cost <- cv$cv.1Std
        print(paste("Selected Cost parameter=", cost, sep=""))
    }
    if(cost==1 & select.prob==TRUE){
        select.prob <- FALSE
    }

    ## Factors only
    all.fac <- all(unlist(lapply(model.frame(formula,data=data)[,-1],
                                 FUN=function(x) is.element("factor",class(x)))))
    if(all.fac==FALSE) stop("Design matrix should contain only factors.")
    rm(all.fac);

    main.formula <- formula
    ## when int.formual != NULL
    if(is.null(int2.formula)==TRUE){
        internal.int <- TRUE
    }else if(is.null(int2.formula)==FALSE){
        internal.int <- FALSE
        if(is.null(int3.formula)==TRUE){
            formula <- as.formula(paste(as.character(formula)[2], "~",
                                        as.character(formula)[3], "+", as.character(int2.formula)[2], sep=""))
        }else if(is.null(int3.formula)==FALSE){
            formula <- as.formula(paste(as.character(formula)[2], "~",
                                        as.character(formula)[3],
                                        "+", as.character(int2.formula)[2],
                                        "+", as.character(int3.formula)[2],
                                        sep=""))
        }
    }
    
    ## ############################
    ## Main Function
    ## ############################ 
    fit <- CausalANOVAFit(formula=formula,
                          internal.int=internal.int,
                          data=data,cost=cost,pair.id=pair.id,
                          family=family, cluster=cluster,
                          screen=screen,screen.type=screen.type,
                          screen.num.int=screen.num.int,
                          maxIter=maxIter,
                          nway=nway,diff=diff,eps=eps,
                          collapse=collapse, collapse.type=collapse.type,
                          collapse.cost=collapse.cost,
                          fac.level=fac.level,ord.fac=ord.fac,
                          verbose=verbose)

    fit <- c(fit, main.formula=main.formula, int2.formula=int2.formula, int3.formula=int3.formula)
    
    if(select.prob==TRUE){
        stab.fit <- stab.CausalANOVA(fit,cluster=cluster,boot=boot)
        output <- list("fit"=fit, "stab.fit"=stab.fit)
        class(output) <- c("CausalANOVA","stab","list")
    }else{
        output <- fit
        class(output) <- c("CausalANOVA","list")
    }        
    return(output)
}
