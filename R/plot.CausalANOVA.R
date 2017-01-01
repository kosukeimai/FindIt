plot.CausalANOVA <- function(x,fac.name,main=NULL,xlab,ylab, ...){

    theme_minimalN <- function (base_size = 12, size.x=20, size.y=20, size.l=20, size.t=20,base_family = "")
        {
            theme_bw(base_size = base_size, base_family = base_family) %+replace%
            theme (
                axis.text.x = element_text(family = base_family, size = size.x * 0.8, lineheight = 0.9, vjust = 1, angle=0,
                    face="bold"),
                axis.text.y = element_text(family = base_family, size = size.y * 0.8, lineheight = 0.9,
                    hjust = 0, angle=0,
                    face="bold"),
                ## axis.text.y = element_blank(), 
                axis.ticks = element_blank(), 
                axis.title.x = element_text(family = base_family, size = size.t, vjust =0.5,angle=0,
                    margin=margin(t=20,l=0,r=0,b =0),face="bold"),
                ## axis.title.y = element_blank(), 
                axis.title.y = element_text(family = base_family, size = size.t, angle = 90, vjust = 0.5,
                    margin=margin(t=0,l=0,r=10,b =0),face="bold"), 
                ## axis.ticks.length = unit(0.3, "lines"), 
                ## axis.ticks.margin = unit(0.5, "lines"), 
                legend.background = element_rect(colour = NA),
                legend.margin = unit(0.4, "cm"), legend.key = element_rect(colour = NA),
                legend.key.size = unit(8.0, "lines"), 
                legend.key.width = NULL,
                legend.key.height = unit(2, "lines"), 
                legend.text = element_text(family = base_family, size = size.l * 0.8, angle=0,face="bold"),      
                legend.text.align = NULL,
                legend.title = element_blank(),
                ## legend.title = element_text(family = base_family, size = size.l * 0.8,
                ##     face = "bold", hjust = 0, angle=0), 
                legend.title.align = NULL,
                legend.position = "top", 
                legend.direction = "horizontal", 
                legend.justification = "center",
                legend.box.just = "top",
                
                panel.background = element_rect(fill = "white", colour = NA), 
                panel.border = element_rect(fill = NA, colour = "white"), 
                panel.grid.major = element_line(colour = "white", size = 0.2), 
                panel.grid.minor = element_line(colour = "white", size = 0.5), 
                panel.margin = unit(0.25, "lines"),
                
                strip.background = element_rect(fill = NA, colour = NA),
                strip.text.x = element_text(family = base_family, size = base_size * 0.8, angle=0, margin(100)), 
                strip.text.y = element_text(family = base_family, size = base_size * 0.8, angle = -90), 
                plot.background = element_rect(colour = NA),
                plot.title = element_text(family = base_family, size = size.t * 1.2, angle=0), 
                plot.margin = unit(c(1, 1, 2, 2), "lines")
            ) 
        }


    TwowayPlot <- function(two.way,i,main,orderx,ordery,VarN1,VarN2){
        
        min.v <- min(two.way[[i]])
        max.v <- max(two.way[[i]])
        
        data.plot <- as.data.frame(cbind(expand.grid(orderx,ordery),c(two.way[[i]])))
        colnames(data.plot) <- c("Var1","Var2","Coefficients")
        
        ## orderx <- as.character(unique(two.way[[i]]$name.1))
        ## ordery <- as.character(unique(two.way[[i]]$name.2))
        
        p <- ggplot(as.data.frame(data.plot), 
                    aes_string(x='Var1',
                        y='Var2',
                        fill='Coefficients')) +
                            ggtitle(main) + 
                                geom_tile(color="black") +
                                    xlab(VarN1) + ylab(VarN2)+
                                        scale_fill_gradient2(low = "#3794bf",high = "#df8640",mid="white",
                                                             limit=c(min.v, max.v))  +
                                                                 theme(legend.text=element_text(size=20),
                                                                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                                                                       axis.title=element_text(size=20,face="bold"),
                                                                       legend.position="top")  + 
                                                                           theme_minimalN(size.x=20,
                                                                                          size.y=20,
                                                                                          size.l=20,size.t=20) +                                           
                                                                                              scale_x_discrete(limits=orderx) + 
                                                                                                  scale_y_discrete(limits=ordery)
        return(p)
    }

    object <- x
    rm(x)

    if(missing(fac.name)){
        warning("Specify two columns of interest.")
    }

    if(class(object)[2]=="stab"){
        fit <- object$fit        
    }else{
        fit <- object
    }
        
    n.fac <- length(fit$fac.level)
    data <- fit$data
    formula <- fit$formula
    
    ## require(ggplot2)
    ## source("TwowayPlot1125.R")
    
    ## TwoPlotR2
    ## if(fac.int[2]<=fac.int[1]){
    ##     warning("fac.int[1] should be smaller than fac.int[2]")
    ## }
    X <- model.frame(formula, data=data)[,-1]
    fac.int <- which(colnames(X) %in% fac.name)
    if(length(fac.int)!=2){
        warning("Use two column names in the data.")
    }
    
    combTwo <- combn(n.fac,2)
    fac.int.ind <- which(apply(combTwo == fac.int,2,function(x) all(x)==TRUE))
    ind.u <- n.fac + fac.int.ind
    
    ## #################
    ## Setup
    ## #################
    
    orderx <- levels(X[,fac.int[1]])
    ordery <- levels(X[,fac.int[2]])

    if(missing(xlab)){
        xlab <- colnames(X)[fac.int[1]]
    }
    if(missing(ylab)){
        ylab <- colnames(X)[fac.int[2]]
    }
    
    ## ####### 
    ## Plot
    ## ####### 
    plotTwo <- TwowayPlot(fit$coefs,i=ind.u,main=main,orderx=orderx,ordery=ordery,
                          VarN1=xlab,VarN2=ylab)

    print(plotTwo)
    cat("Note: Baseline is the mean of the AMIEs as in the figures of Egami and Imai(2016+). \n")
}
