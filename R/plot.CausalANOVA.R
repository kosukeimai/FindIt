#' Plotting CausalANOVA
#' @param x An output from \code{CausalANOVA}
#' @param fac.name Factor names to plot. Length should be 2.
#' @param treat.ind Which factor serves as the main treatment. Should be 1 (the first element of \code{fac.name}) or 2 (the second element of \code{fac.name}).
#' @param type What types of effects to plot. Should be one of \code{AME}, \code{AMIE} and \code{ConditionalEffect}.
#' @param space Space on the left side of the plot. 
#' @param xlim Range for the x-axis
#' @param ... Other graphical parameters
#' @export
plot.CausalANOVA <- function(x, fac.name,
                             treat.ind=1,
                             type = "ConditionalEffect", space = 15, xlim, ...){
    treat.ind <-  1
    adj.p   <- 2.3
    add     <- FALSE
    col <-  "black"
    ## HouseKeeping 
    if(x$Gorder==1 | missing(fac.name)==TRUE){
        type <- "AME"
    }
    
    if(type=="AME"){
        formula <- x$formula
        terms.f <- terms(formula, data=x$data)
        order.f <- attr(terms.f, "order")
        var.name <- attr(terms.f, "term.labels")[order.f==1]
        level.list <- x$level.list        
        inference <- x$inference
        
        if(inference==TRUE){
            AME.CI  <- x$CI.table[which(order.f==1)]
            font <- CI.value <- CI.up <- CI.low <- c()
            Name <- list()
            for(z in 1:length(AME.CI)){
                Name0 <- c(paste("Factor: ", var.name[z], sep=""),
                           paste("    ", level.list[[z]],sep=""))
                Name <- c(Name, Name0)
                font <- c(font, c(2, rep(1, length(level.list[[z]]))))
                AME.base <- AME.CI[[z]]
                AME.base[,3] <- AME.base[,1] -1.96*AME.base[,2]
                AME.base[,4] <- AME.base[,1] +1.96*AME.base[,2]
                CI.value <- c(CI.value, c(NA, AME.base[,1]))
                CI.up <- c(CI.up, c(NA, AME.base[,4]))
                CI.low <- c(CI.low, c(NA, AME.base[,3]))               
            }
            xlim.min <- min(na.omit(CI.low))
            xlim.max <- max(na.omit(CI.up))
        }else{
            AME.CI  <- x$AME
            font <- CI.value <- c()
            Name <- list()
            for(z in 1:length(x$AME)){
                Name0 <- c(paste("Factor: ", var.name[z], sep=""),
                           paste("    ", level.list[[z]],sep=""))
                Name <- c(Name, Name0)
                font <- c(font, c(2, rep(1, length(level.list[[z]]))))
                ame <- unlist(AME.CI[[z]])
                ame <- ame - ame[length(ame)]
                CI.value <- c(CI.value, NA, ame)
            }
            xlim.min <- min(na.omit(CI.value))
            xlim.max <- max(na.omit(CI.value))
        }
        
        if(missing(xlim)==FALSE){
            par(mar=c(4,space,2,2))
            plot(rev(CI.value), seq(from=1, to=length(CI.value)), pch=19, yaxt="n", xlim=xlim,
                 ylim=c(0, length(CI.value)+1), ylab="", main="AME", xlab="Estimated Effects", col=col)
            if(inference==TRUE){
                segments(rev(CI.low), seq(from=1, to=length(CI.value)),
                         rev(CI.up), seq(from=1, to=length(CI.value)), lwd=2,col=col)
            }
            Axis(side=2, at=seq(from=1, to=length(CI.value)), labels=rep("",length(Name)))
            Axis(side=2, at=seq(from=1, to=length(CI.value)), labels=rev(Name),
                 las=2, hadj=0, line=10, tck=0, lwd=0)
            abline(v=0, lty=2)
        }else{
            par(mar=c(4,space,2,2))
            plot(rev(CI.value), seq(from=1, to=length(CI.value)), pch=19, yaxt="n",
                 xlim=c(xlim.min, xlim.max),
                 ylim=c(0, length(CI.value)+1), ylab="", main="AME", xlab="Estimated Effects", col=col)
            if(inference==TRUE){
                segments(rev(CI.low), seq(from=1, to=length(CI.value)),
                         rev(CI.up), seq(from=1, to=length(CI.value)), lwd=2,col=col)
            }
            Axis(side=2, at=seq(from=1, to=length(CI.value)), labels=rep("",length(Name)))
            Axis(side=2, at=seq(from=1, to=length(CI.value)), labels=rev(Name),
                 las=2, hadj=0, line=10, tck=0, lwd=0)
            abline(v=0, lty=2)
        }        
        
    }else if(type=="AMIE"){
        formula <- x$formula
        indTwo <- x$indTwo
        indThree <- x$indThree
        plot.Gorder <- length(fac.name)
        terms.f <- terms(formula,data=x$data)
        order.f <- attr(terms.f, "order")
        var.name <- attr(terms.f, "term.labels")[order.f==1]
        inference <- x$inference
        AMIE2 <- x$AMIE2
        AMIE3 <- x$AMIE3
        
        ## This accomodate Three-ways 
        Fac.ind <- c()
        for(z in 1:length(fac.name)){
            Fac.ind[z] <- which(fac.name[z]==var.name)
        }
        Fac.ind <- Fac.ind[order(Fac.ind)]
        
        ## Find the index
        if(plot.Gorder==2){
            INT.ind <- which(apply(c(Fac.ind[1], Fac.ind[2]) == indTwo, 2, all)==TRUE)
            if(inference==TRUE) AMIE.CI  <- x$CI.table[which(order.f==2)]
            else AMIE.CI  <- AMIE2
        }
        if(plot.Gorder==3){
            INT.ind <- which(apply(c(Fac.ind[1], Fac.ind[2], Fac.ind[3]) == indThree, 2, all)==TRUE)
            if(inference==TRUE) AMIE.CI  <- x$CI.table[which(order.f==3)]
            else AMIE.CI  <- AMIE3
        }

        ## Plot: AMIE (just need fac.name)
        if(inference==TRUE){
            CI.value <- AMIE.CI[[INT.ind]][,1]
            CI.up <- AMIE.CI[[INT.ind]][,4]
            CI.low <-  AMIE.CI[[INT.ind]][,3]
            Name   <- rownames(AMIE.CI[[INT.ind]])
            xlim.min <- min(na.omit(CI.low))
            xlim.max <- max(na.omit(CI.up))
        }else{
            amie <- unlist(AMIE.CI[[INT.ind]])
            CI.value <- amie - amie[length(amie)]
            Name <- names(CI.value)
            xlim.min <- min(na.omit(CI.value))
            xlim.max <- max(na.omit(CI.value))
        }
        
        
        ## plot.value <- plot.value - plot.value[base.ind]
        
        max.nchar <- max(nchar(Name))/adj.p
        ## Create blank plot first
        if(missing(xlim)==FALSE){
            par(mar=c(4,space,2,2))
            plot(rev(CI.value), seq(from=1, to=length(CI.value)), pch=19, yaxt="n", xlim=xlim,
                 ylim=c(0, length(CI.value)+1), ylab="", main="AMIE", xlab="Estimated Effects",col=col)
            if(inference==TRUE){
                segments(rev(CI.low), seq(from=1, to=length(CI.value)),
                         rev(CI.up), seq(from=1, to=length(CI.value)), lwd=2, col=col)
            }
            Axis(side=2, at=seq(from=1, to=length(CI.value)), labels=rep("",length(Name)))
            Axis(side=2, at=seq(from=1, to=length(CI.value)), labels=rev(Name),
                 las=2, hadj=0, line=max.nchar, tck=0, lwd=0)
            abline(v=0, lty=2)                         
        }else{
            par(mar=c(4,space,2,2))
            plot(rev(CI.value), seq(from=1, to=length(CI.value)), pch=19, yaxt="n",
                 xlim=c(xlim.min, xlim.max),
                 ylim=c(0, length(CI.value)+1), ylab="", main="AMIE", xlab="Estimated Effects",col=col)
            if(inference==TRUE){
                segments(rev(CI.low), seq(from=1, to=length(CI.value)),
                         rev(CI.up), seq(from=1, to=length(CI.value)), lwd=2, col=col)
            }
            Axis(side=2, at=seq(from=1, to=length(CI.value)), labels=rep("",length(Name)))
            Axis(side=2, at=seq(from=1, to=length(CI.value)), labels=rev(Name),
                 las=2, hadj=0, line=max.nchar, tck=0, lwd=0)
            abline(v=0, lty=2)
        }
        
    }else if(type=="ConditionalEffect"){      
        treat.fac <- fac.name[treat.ind]
        if(treat.ind==2) cond.ind <- 1
        else if(treat.ind==1) cond.ind <- 2
        cond.fac  <- fac.name[cond.ind]
        inference <- x$inference
        CE <- ConditionalEffect(object=x, treat.fac=treat.fac, inference=inference,
                                cond.fac=cond.fac, base.ind=NULL, verbose=FALSE)
        CE.plot <- CE$Conditional
        cond.level <- CE$cond.level
        treat.level <- CE$treat.level

        if(inference==TRUE){
            font <- CE.value <- CE.up <- CE.low <- c()
            Name <- list()
            for(z in 1:length(CE.plot)){
                Name0 <- c(paste(cond.fac, "=", cond.level[z]),
                           paste("    ", treat.level,sep=""))
                Name <- c(Name, Name0)
                font <- c(font, c(2, rep(1, length(treat.level))))
                CE.value <- c(CE.value, c(NA, CE.plot[[z]][,1]))
                CE.up <- c(CE.up, c(NA, CE.plot[[z]][,4]))
                CE.low <- c(CE.low, c(NA, CE.plot[[z]][,3]))
            }
            xlim.min <- min(na.omit(CE.low))
            xlim.max <- max(na.omit(CE.up))

            if(missing(xlim)==FALSE){
                par(mar=c(4,space,2,2))
                plot(rev(CE.value), seq(from=1, to=length(CE.value)), pch=19, yaxt="n", xlim=xlim,
                     ylim=c(0, length(CE.value)+1), ylab="", xlab = "Estimated Effects", 
                     col=col, main = "Conditional Effects")
                segments(rev(CE.low), seq(from=1, to=length(CE.value)),
                         rev(CE.up), seq(from=1, to=length(CE.value)), lwd=2, col=col)
                Axis(side=2, at=seq(from=1, to=length(CE.value)), labels=rep("",length(Name)))
                Axis(side=2, at=seq(from=1, to=length(CE.value)), labels=rev(Name),
                     las=2, hadj=0, line=10, tck=0, lwd=0)
                abline(v=0, lty=2)
            }else{
                par(mar=c(4,space,2,2))
                plot(rev(CE.value), seq(from=1, to=length(CE.value)), pch=19, yaxt="n",
                     xlim=c(xlim.min, xlim.max),
                     ylim=c(0, length(CE.value)+1), ylab="", 
                     xlab="Estimated Effects", main = "Conditional Effects", col=col)
                segments(rev(CE.low), seq(from=1, to=length(CE.value)),
                         rev(CE.up), seq(from=1, to=length(CE.value)), lwd=2, col=col)
                Axis(side=2, at=seq(from=1, to=length(CE.value)), labels=rep("",length(Name)))
                Axis(side=2, at=seq(from=1, to=length(CE.value)), labels=rev(Name),
                     las=2, hadj=0, line=10, tck=0, lwd=0)
                abline(v=0, lty=2) 
            }
        }else{
            font <- CE.value <- c()
            Name <- list()
            for(z in 1:length(CE.plot)){
                Name0 <- c(paste(cond.fac, "=", cond.level[z]),
                           paste("    ", treat.level,sep=""))
                Name <- c(Name, Name0)
                font <- c(font, c(2, rep(1, length(treat.level))))
                CE.value <- c(CE.value, c(NA, CE.plot[[z]]))
            }
            xlim.min <- min(na.omit(CE.value))
            xlim.max <- max(na.omit(CE.value))

            if(missing(xlim)==FALSE){
                if(add==FALSE){
                    par(mar=c(4,space,2,2))
                }                    
                plot(rev(CE.value), seq(from=1, to=length(CE.value)), pch=19, yaxt="n", xlim=xlim,
                     ylim=c(0, length(CE.value)+1), ylab="", xlab="Estimated Effects", 
                     col=col, main = "Conditional Effects")
                Axis(side=2, at=seq(from=1, to=length(CE.value)), labels=rep("",length(Name)))
                Axis(side=2, at=seq(from=1, to=length(CE.value)), labels=rev(Name),
                     las=2, hadj=0, line=10, tck=0, lwd=0)
                abline(v=0, lty=2)
            }else{
                if(add==FALSE){
                    par(mar=c(4, space, 2,2))
                } 
                plot(rev(CE.value), seq(from=1, to=length(CE.value)), pch=19, yaxt="n",
                     xlim=c(xlim.min, xlim.max),
                     ylim=c(0, length(CE.value)+1), ylab="", 
                     xlab="Estimated Effects", col=col,
                     main = "Conditional Effects")
                Axis(side=2, at=seq(from=1, to=length(CE.value)), labels=rep("",length(Name)))
                Axis(side=2, at=seq(from=1, to=length(CE.value)), labels=rev(Name),
                     las=2, hadj=0, line=10, tck=0, lwd=0)
                abline(v=0, lty=2) 
            }
        }
    }    
}
