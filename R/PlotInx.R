#' Internal function:
#'
#'
#' @export

FilterInx_step1 <- function(INX,cellTypeA,cellTypeB,proteinTypeA,proteinTypeB) {
  if (cellTypeA == cellTypeB) {
    INX <- INX[[which(sapply(strsplit(names(INX),"~"),function(X) all(X == cellTypeA)))]]
  } else {
    INX <- INX[[which(sapply(strsplit(names(INX),"~"),function(X) cellTypeA %in% X & cellTypeB %in% X))]]
  }

  nodesA <- INX$nodes[INX$nodes$cellType == cellTypeA &
                        grepl(proteinTypeA,INX$nodes$proteinType),]
  if (nrow(nodesA) < 1) {
    stop(paste0("No ",proteinTypeA," in ",cellTypeA,"."))
  }
  nodesA$side <- "A"
  rownames(nodesA) <- paste(rownames(nodesA),"A",sep="_")
  nodesB <- INX$nodes[INX$nodes$cellType == cellTypeB &
                        grepl(proteinTypeB,INX$nodes$proteinType),]
  if (nrow(nodesA) < 1) {
    stop(paste0("No ",proteinTypeB," in ",cellTypeB,"."))
  }
  nodesB$side <- "B"
  rownames(nodesB) <- paste(rownames(nodesB),"B",sep="_")
  INX$nodes <- rbind(nodesA,nodesB)

  edgesA <- INX$edges[(INX$edges$nodeA %in% INX$nodes$node[INX$nodes$side == "A"] &
                         INX$edges$nodeB %in% INX$nodes$node[INX$nodes$side == "B"]),]
  if (nrow(edgesA) > 0) {
    edgesA$nodeA <- paste(edgesA$nodeA,"A",sep="_")
    edgesA$nodeB <- paste(edgesA$nodeB,"B",sep="_")
  }
  edgesB <- INX$edges[(INX$edges$nodeB %in% INX$nodes$node[INX$nodes$side == "A"] &
                         INX$edges$nodeA %in% INX$nodes$node[INX$nodes$side == "B"]),]
  if (nrow(edgesB) > 0) {
    tempB <- paste(edgesB$nodeA,"B",sep="_")
    tempA <- paste(edgesB$nodeB,"A",sep="_")
    edgesB$nodeA <- tempA
    edgesB$nodeB <- tempB
  }
  INX$edges <- rbind(edgesA,edgesB)

  INX$nodes <- INX$nodes[rownames(INX$nodes) %in% unique(c(INX$edges$nodeA,INX$edges$nodeB)),]
  INX$nodes <- INX$nodes[,-which(colnames(INX$nodes) == "node")]

  attr(INX,"cellType") <- list(A=cellTypeA,B=cellTypeB)
  attr(INX,"proteinType") <- list(A=proteinTypeA,B=proteinTypeB)
  return(INX)
}


#' Internal function:
#'
#'
#' @export

FilterInx_topN <- function(INX,topN) {
  INX$edges <- INX$edges[head(order(abs(INX$edges$meanDEscaled),decreasing=T),topN),]
  INX$nodes <- INX$nodes[unique(c(INX$edges$nodeA,INX$edges$nodeB)),]
  return(INX)
}


#' Internal function:
#'
#'
#' @export

FilterInx_GeneStatistic <- function(INX,statThresh) {
  temp_nodes <- rownames(INX$nodes)[INX$nodes[,attr(INX,"GeneStatistic")] <= statThresh]
  INX$edges <- INX$edges[INX$edges$nodeA %in% temp_nodes | INX$edges$nodeB %in% temp_nodes,]
  INX$nodes <- INX$nodes[unique(c(INX$edges$nodeA,INX$edges$nodeB)),]
  return(INX)
}


#' Internal function:
#'
#'
#' @export

FilterInx_GeneMagnitude <- function(INX,magnThresh) {
  temp_nodes <- rownames(INX$nodes)[abs(INX$nodes[,attr(INX,"GeneMagnitude")]) > magnThresh]
  INX$edges <- INX$edges[INX$edges$nodeA %in% temp_nodes | INX$edges$nodeB %in% temp_nodes,]
  INX$nodes <- INX$nodes[unique(c(INX$edges$nodeA,INX$edges$nodeB)),]
  return(INX)
}


#' Internal function:
#'
#'
#' @export

FilterInx_genenames <- function(INX,genenames) {
  temp_nodes <- rownames(INX$nodes)[INX$nodes$gene %in% genenames]
  INX$edges <- INX$edges[INX$edges$nodeA %in% temp_nodes | INX$edges$nodeB %in% temp_nodes,]
  INX$nodes <- INX$nodes[unique(c(INX$edges$nodeA,INX$edges$nodeB)),]
  return(INX)
}


#' Internal function:
#'
#'
#' @export

DoPlotInx <- function(INX,ySpacing) {
  yoCol <- colorRampPalette(rgb(red=c(247,220,87),green=c(78,220,87),blue=c(214,220,249),
                                names=c("old","none","young"),maxColorValue=255),
                            bias=1,interpolate="linear")
  colourScheme <- yoCol(100)
  # colourScheme <- colorspace::diverge_hcl(100)

  if (missing(ySpacing)) {
    ySpacing <- "relative"
  }
  if (!ySpacing %in% c("absolute","relative")) {
    warning("ySpacing must be one of: 'absolute' or 'relative'. Set to 'relative'.")
    ySpacing <- "relative"
  }
  if (nrow(INX$nodes) < 1 | nrow(INX$edges) < 1) {
    stop("No genes passed filters.")
  }

  INX$nodes$x[INX$nodes$side == "A"] <- 1
  INX$nodes$x[INX$nodes$side == "B"] <- 3

  if (ySpacing == "relative") {
    temp <- INX$nodes[,attr(INX,"GeneMagnitude")]
    temp[temp == Inf] <- max(temp[temp < Inf]) * 1.1
    temp[temp == -Inf] <- min(temp[temp > -Inf]) * 1.1
    p <- INX$nodes$side == "A"
    temp[p] <- seq(min(temp),max(temp),length.out=sum(p))[rank(temp[p],ties.method="first")]
    p <- INX$nodes$side == "B"
    temp[p] <- seq(min(temp),max(temp),length.out=sum(p))[rank(temp[p],ties.method="first")]
    INX$nodes$y <- temp
  } else {
    INX$nodes$y <- INX$nodes[,attr(INX,"GeneMagnitude")]
    INX$nodes$y[INX$nodes$y == Inf] <- max(INX$nodes$y[INX$nodes$y < Inf]) * 1.1
    INX$nodes$y[INX$nodes$y == -Inf] <- min(INX$nodes$y[INX$nodes$y > -Inf]) * 1.1
  }

  temp_b <- INX$nodes[,attr(INX,"GeneMagnitude")] <= 0
  if (any(temp_b)) {
    if (any(is.infinite(INX$nodes[temp_b,attr(INX,"GeneMagnitude")]))) {
      temp_bc <- INX$nodes[temp_b,attr(INX,"GeneMagnitude")]
      temp_bc[is.infinite(temp_bc)] <- min(temp_bc[!is.infinite(temp_bc)]) * 1.1
      temp_bc <- cut(c(0,temp_bc),50,labels=F)[-1]
    } else {
      temp_bc <- cut(c(0,INX$nodes[temp_b,attr(INX,"GeneMagnitude")]),50,labels=F)[-1]
    }
    INX$nodes$col[temp_b] <- temp_bc
  }
  temp_a <- INX$nodes[,attr(INX,"GeneMagnitude")] > 0
  if (any(temp_a)) {
    if (any(is.infinite(INX$nodes[temp_a,attr(INX,"GeneMagnitude")]))) {
      temp_ac <- INX$nodes[temp_a,attr(INX,"GeneMagnitude")]
      temp_ac[is.infinite(temp_ac)] <- max(temp_ac[!is.infinite(temp_ac)]) * 1.1
      temp_ac <- cut(c(0,temp_ac),50,labels=F)[-1]
    } else {
      temp_ac <- cut(c(0,INX$nodes[temp_a,attr(INX,"GeneMagnitude")]),50,labels=F)[-1]
    }
    INX$nodes$col[temp_a] <- temp_ac + 50
  }

  INX$nodes$signif <- cut(INX$nodes[,attr(INX,"GeneStatistic")],
                          breaks=c(1,.05,.01,.001,.0001,0),
                          right=T,include.lowest=T)

  INX$edges$col <- cut(c(1,-1,INX$edges$meanDEscaled),100,labels=F)[-1:-2]
  INX$edges$lwd <- seq(2,6)[cut(c(0,1,abs(INX$edges$meanDEscaled)),5,labels=F)[-1:-2]]


  par(mar=c(3,3,1,1),mgp=2:0)
  plot(x=NULL,y=NULL,xlim=c(0,7),ylim=range(INX$nodes$y),
       xaxs="i",xaxt="n",yaxt="n",bty="n",
       xlab=NA,ylab=NA)
  temp_junk <- sapply(rownames(INX$edges),function(X)
    lines(x=INX$nodes[unlist(INX$edges[X,c("nodeA","nodeB")]),"x"],
          y=INX$nodes[unlist(INX$edges[X,c("nodeA","nodeB")]),"y"],
          col=scales::alpha(colourScheme,
                            seq(.3,.8,length.out=100))[INX$edges[X,"col"]],
          lwd=INX$edges[X,"lwd"])
  )
  points(x=INX$nodes$x,y=INX$nodes$y,
         pch=19,cex=2,
         col=colourScheme[INX$nodes$col])
  points(x=INX$nodes$x,y=INX$nodes$y,
         pch=1,cex=2,lwd=2,
         col=c("gray0","gray25","gray50","gray75","gray100")[INX$nodes$signif])
  text(x=INX$nodes$x[INX$nodes$side == "A"],
       y=INX$nodes$y[INX$nodes$side == "A"],
       labels=INX$nodes$gene[INX$nodes$side == "A"],
       pos=2,col="black")
  text(x=INX$nodes$x[INX$nodes$side == "B"],
       y=INX$nodes$y[INX$nodes$side == "B"],
       labels=INX$nodes$gene[INX$nodes$side == "B"],
       pos=4,col="black")
  mtext(unlist(attr(INX,"cellType")),side=1,line=.5,font=2,
        at=c(unique(INX$nodes$x[INX$nodes$side == "A"]),
             unique(INX$nodes$x[INX$nodes$side == "B"])))
  mtext(unlist(attr(INX,"proteinType")),side=1,line=1.5,font=2,
        at=c(unique(INX$nodes$x[INX$nodes$side == "A"]),
             unique(INX$nodes$x[INX$nodes$side == "B"])))
  if (ySpacing == "absolute") {
    axis(2,pos=0)
    mtext(attr(INX,"GeneMagnitude"),font=2,side=2,line=2)
  }

  legend(x=4.5,y=par("usr")[4] - (par("usr")[4] - par("usr")[3]) * .1,
         bty="n",pch=1,pt.lwd=2,pt.cex=2,
         legend=c("> 0.05","0.01 to 0.05","0.001 to 0.01","0.0001 to 0.001","< 0.0001"),
         col=c("gray100","gray75","gray50","gray25","gray0"))
  text(labels=attr(INX,"GeneStatistic"),x=4.5,
       y=par("usr")[4] - (par("usr")[4] - par("usr")[3]) * .1,
       font=2,adj=c(-.1,.5))

  temp_fc <- INX$nodes$col[!is.infinite(INX$nodes[,attr(INX,"GeneMagnitude")])]
  if (any(INX$nodes[,attr(INX,"GeneMagnitude")] <= 0)) {
    rect(xleft=4.6,xright=5,
         ybottom=seq(from=par("usr")[3] + (par("usr")[4] - par("usr")[3]) * .1,
                     to=par("usr")[3] + (par("usr")[4] - par("usr")[3]) * .3,
                     length.out=50 - min(temp_fc) + 1)[seq(1,50 - min(temp_fc))],
         ytop=seq(from=par("usr")[3] + (par("usr")[4] - par("usr")[3]) * .1,
                  to=par("usr")[3] + (par("usr")[4] - par("usr")[3]) * .3,
                  length.out=50 - min(temp_fc) + 1)[seq(1,50 - min(temp_fc)) + 1],
         col=colourScheme[seq(min(temp_fc),50)],border=NA)
    text(x=5,y=par("usr")[3] + (par("usr")[4] - par("usr")[3]) * .1,
         labels=round(min(INX$nodes[,attr(INX,"GeneMagnitude")][
           !is.infinite(INX$nodes[,attr(INX,"GeneMagnitude")])]),2),
         adj=c(-.1,0))
  }
  if (any(INX$nodes[,attr(INX,"GeneMagnitude")] > 0)) {
    rect(xleft=4.6,xright=5,
         ybottom=seq(from=par("usr")[3] + (par("usr")[4] - par("usr")[3]) * .3,
                     to=par("usr")[3] + (par("usr")[4] - par("usr")[3]) * .5,
                     length.out=max(temp_fc) - 51 + 1)[seq(1,max(temp_fc) - 51)],
         ytop=seq(from=par("usr")[3] + (par("usr")[4] - par("usr")[3]) * .3,
                  to=par("usr")[3] + (par("usr")[4] - par("usr")[3]) * .5,
                  length.out=max(temp_fc) - 51 + 1)[seq(1,max(temp_fc) - 51) + 1],
         col=colourScheme[seq(51,max(temp_fc))],border=NA)
    text(x=5,y=par("usr")[3] + (par("usr")[4] - par("usr")[3]) * .5,
         labels=round(max(INX$nodes[,attr(INX,"GeneMagnitude")][
           !is.infinite(INX$nodes[,attr(INX,"GeneMagnitude")])]),2),
         adj=c(-.1,1))
  }
  text(x=5,y=par("usr")[3] + (par("usr")[4] - par("usr")[3]) * .3,
       labels=0,adj=c(-0.5,0.5))
  if (any(INX$nodes[,attr(INX,"GeneMagnitude")] == Inf)) {
    rect(xleft=4.6,xright=5,
         ybottom=par("usr")[3] + (par("usr")[4] - par("usr")[3]) * .51,
         ytop=par("usr")[3] + (par("usr")[4] - par("usr")[3]) * .54,
         col=colourScheme[100],border=NA)
    text(x=5,y=par("usr")[3] + (par("usr")[4] - par("usr")[3]) * .54,
         labels=Inf,adj=c(-.2,1))
  }
  if (any(INX$nodes[,attr(INX,"GeneMagnitude")] == -Inf)) {
    rect(xleft=4.6,xright=5,
         ybottom=par("usr")[3] + (par("usr")[4] - par("usr")[3]) * .09,
         ytop=par("usr")[3] + (par("usr")[4] - par("usr")[3]) * .06,
         col=colourScheme[1],border=NA)
    text(x=5,y=par("usr")[3] + (par("usr")[4] - par("usr")[3]) * .06,
         labels=-Inf,adj=c(-.2,0))
  }
  text(x=4.6,y=par("usr")[3] + (par("usr")[4] - par("usr")[3]) *
         switch(as.character(any(INX$nodes[,attr(INX,"GeneMagnitude")] == Inf)),
                "TRUE"=.54,
                "FALSE"=.5),
       labels=attr(INX,"GeneMagnitude"),font=2,adj=c(0,-1))
}


#' Plot cell-cell interactions as bipartite graph
#'
#'
#' @export

PlotCCInx <- function(INX,cellTypeA,cellTypeB,proteinTypeA,proteinTypeB,
                      GeneMagnitudeThreshold,GeneStatisticThreshold,
                      TopEdges,GeneNames,YSpacing="relative") {
  if (missing(INX) |
      missing(cellTypeA) |
      missing(cellTypeB) |
      missing(proteinTypeA) |
      missing(proteinTypeB)) {
    stop("The following arguments are required: INX, cellTypeA, cellTypeB, proteinTypeA, proteinTypeB.")
  }

  INX <- FilterInx_step1(INX,
                         cellTypeA=cellTypeA,
                         cellTypeB=cellTypeB,
                         proteinTypeA=proteinTypeA,
                         proteinTypeB=proteinTypeB)

  if (!missing(GeneMagnitudeThreshold)) {
    INX <- FilterInx_GeneMagnitude(INX,GeneMagnitudeThreshold)
  } else if (!missing(GeneStatisticThreshold)) {
    INX <- FilterInx_GeneStatistic(INX,GeneStatisticThreshold)
  } else if (!missing(TopEdges)) {
    INX <- FilterInx_topN(INX,TopEdges)
  } else if (!missing(GeneNames)) {
    INX <- FilterInx_genenames(INX,GeneNames)
  }

  DoPlotInx(INX,YSpacing)
}

