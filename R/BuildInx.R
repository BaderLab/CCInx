
#' Score differential expression gene statistics from cell type clusters
#'
#' This function takes differential expression gene statistics from scRNAseq
#' data representing a cell type as a dataframe, and assigns scaled scores to
#' the statistic of choice. This is used to rank nodes and edges by differential
#' expression when viewing the bipartite ligand-receptor plots.
#'
#' @param gdb A data frame representing gene statistics from a cell type, where
#'   each row is a gene with official gene symbols as row names. Variables
#'   should be appropriately named statistics or annotations to be included in
#'   the resulting node metadata.
#' @param DEmagn A character vector of length 1 representing the variable name
#'   in the GeneStatList data frames carrying information on the magnitude and
#'   direction of the change of expression for the node (gene) in each cell
#'   type. This is generally a signed logFC or gene expression ratio.
#' @param DEstat A character vector of length 1 representing the variable name
#'   in the GeneStatList data frames carrying information on the statistical
#'   significance of expression change. This is generally a corrected p-value.
#'
#' @export

ScoreDE <- function(gdb,DEmagn,DEstat) {
  if (any(gdb[[DEstat]] <= 0)) {
    gdb[[DEstat]][gdb[[DEstat]] == 0] <- 10^-ceiling(-log10(min(gdb[[DEstat]][gdb[[DEstat]] > 0])))
  }
  gdb$score <- sign(gdb[[DEmagn]]) * -log10(gdb[[DEstat]])
  gdb$scoreScaled <- sapply(gdb$score,function(X)
    X / switch(as.character(X >= 0),
               "TRUE"=max(gdb$score),
               "FALSE"=min(gdb$score) * -1))
  gdb$scoreScaled[is.na(gdb$scoreScaled)] <- 0
  if (length(unique(gdb$scoreScaled)) == 1) {
    gdb$scaledFactor <- 6
  } else {
    gdb$scaledFactor <- 6
    if (any(gdb$scoreScaled < 0)) {
      gdb$scaledFactor[gdb$scoreScaled < 0] <-
        cut(gdb$scoreScaled[gdb$scoreScaled < 0],5,labels=F)
    }
    if (any(gdb$scoreScaled > 0)) {
      gdb$scaledFactor[gdb$scoreScaled > 0] <-
        gdb$scaledFactor[gdb$scoreScaled > 0] +
        cut(gdb$scoreScaled[gdb$scoreScaled > 0],5,labels=F)
    }
  }
  return(gdb)
}


ScoreExpr <- function(gdb,expr) {
  stop("Only DE scoring is supported right now")
}


#' Score gene nodes from all cell type clusters
#'
#' This function takes cluster-wise gene statistics from scRNAseq data as a list
#' of dataframes, and assigns scaled scores to the statistic of choice. Two
#' types of networks can be created with CCInx, and this is where the decision
#' is made. Networks can be created based on detection of genes, or differential
#' expression between conditions. If the latter, \code{GeneStatistic} should be
#' specified.
#'
#' @param GeneStatList A named list of dataframes. Each list element should
#'   represent a cell type / cluster to be included in the interaction network,
#'   and should be named accordingly. List elements should contain data frames
#'   where each row is a gene with official gene symbols as row names. Variables
#'   should be appropriately named statistics or annotations to be included in
#'   the resulting node metadata. Variable names should be consistent between
#'   list elements.
#' @param GeneMagnitude A character vector of length 1 representing the variable
#'   name in the GeneStatList data frames carrying information on the magnitude
#'   (and direction of the change) of expression for the node (gene) in each
#'   cell type. This is either a measure of expression (generally mean
#'   expression or detection rate) or a measure of change (signed log expression
#'   ratio a.k.a. logFC).
#' @param GeneStatistic Optional. A character vector of length 1 representing
#'   the variable name in the GeneStatList data frames carrying information on
#'   the statistical significance of expression change. This is generally a
#'   corrected p-value.
#'
#' @export

ScoreClustGenes <- function(GeneStatList,GeneMagnitude,GeneStatistic) {
  if (any(sapply(GeneStatList,function(X) any(is.na(X[,GeneMagnitude])))) |
      any(sapply(GeneStatList,function(X) any(is.na(X[,GeneStatistic]))))) {
    stop(paste("This function can't handle missing values (NA) in",GeneMagnitude,GeneStatistic))
  }
  if (missing(GeneStatistic)) {
    message("Scoring gene expression for all cell types...")
    return(pbapply::pbsapply(GeneStatList,
                             ScoreExpr,
                             expr=GeneMagnitude,
                             simplify=F))
  } else {
    message("Scoring differential expression for all cell types...")
    return(pbapply::pbsapply(GeneStatList,
                             ScoreDE,
                             DEmagn=GeneMagnitude,
                             DEstat=GeneStatistic,
                             simplify=F))
  }
}



#' Build cell-cell interaction edge list from scored cell-type statistics
#'
#' This function uses the output of \code{\link{ScoreClustGenes}} to build a
#' predicted cell-cell interaction network. It returns a list of bipartite
#' graphs (as edge lists) for each combination of cell types.
#'
#' @param SCL The output of \code{\link{ScoreClustGenes}}.
#' @param NodeScoreThreshold Default=0. The minimum absolute score required to
#'   include a node in the cell-cell interaction network. The default does not
#'   filter any nodes.
#' @param EdgeScoreThreshold Default=0. The minimum absolute score required to
#'   include an edge in the cell-cell interaction network. Absolute edge score
#'   is the summed magnitude of node scores. The default does not filter any
#'   edges.
#' @param species Default='hsapiens'. The species of the source data. One of
#'   'hsapiens' or 'mmusculus'. Note that the ligand-receptor database was built
#'   for human, and the mouse version is generated by homology mapping (only
#'   using uniquely mapped homologues).
#'
#' @export

BuildInx <- function(SCL,NodeScoreThreshold=0,EdgeScoreThreshold=0,species="hsapiens") {
  switch(species,
         hsapiens=load(system.file("LigRecDB_RData/BaderCCIeditedbyBI.RData",
                                   package="CCInx")),
         mmusculus=load(system.file("LigRecDB_RData/BaderCCIeditedbyBI_mouse.RData",
                                    package="CCInx")),
         stop("species must be one of 'hsapiens' or 'mmusculus'."))

  tempCN <- c()
  names(SCL) <- gsub("~","_",names(SCL))
  for (a in names(SCL)) {
    for (b in names(SCL)) {
      temp <- paste(sort(c(a,b)),collapse="~")
      if (!temp %in% tempCN) {
        tempCN <- append(tempCN,temp)
      }
    }
  }
  rm(a,b)
  tempComp <- strsplit(tempCN,"~")
  names(tempComp) <- tempCN

  message("Building predicted cell-cell interaction networks...")
  inxL <- pbapply::pbsapply(tempComp[1:3],function(Z) {
    a <- Z[1]; b <- Z[2]
    temp_a <- rownames(SCL[[a]])[abs(SCL[[a]]$score) >= NodeScoreThreshold]
    temp_b <- rownames(SCL[[b]])[abs(SCL[[b]]$score) >= NodeScoreThreshold]
    if (length(temp_a) < 1 | length(temp_b) < 1) { return(NULL) }
    inx <- c()
    keysAB <- inxDB$key[inxDB$nodeA %in% temp_a & inxDB$nodeB %in% temp_b]
    sumAbsAB <- sapply(strsplit(keysAB,"_"),function(X)
      sum(abs(c(SCL[[a]][X[1],"score"],SCL[[b]][X[2],"score"]))))
    sumScaledAB <- sapply(strsplit(keysAB,"_"),function(X)
      sum(c(SCL[[a]][X[1],"scoreScaled"],SCL[[b]][X[2],"scoreScaled"])))
    if (any(sumAbsAB >= EdgeScoreThreshold)) {
      keysAB <- keysAB[sumAbsAB >= EdgeScoreThreshold]
      namesAB <- sapply(strsplit(keysAB,"_"),function(X)
        paste(paste(X[1],a,sep="_"),paste(X[2],b,sep="_"),sep="~"))
      temp_inx <- inxDB[keysAB,c("key","nodeA","nodeB")]
      rownames(temp_inx) <- temp_inx$key <- namesAB
      colnames(temp_inx)[2:3] <- c("geneA","geneB")
      temp_inx$nodeA <- paste(temp_inx$geneA,a,sep="_")
      temp_inx$nodeB <- paste(temp_inx$geneB,b,sep="_")
      temp_inx <- temp_inx[,c(1,4,5,2,3)]
      temp_inx$cellTypeA <- a
      temp_inx$cellTypeB <- b
      temp_inx$proteinTypeA <- geneInfo[temp_inx$geneA,"protein_type"]
      temp_inx$proteinTypeB <- geneInfo[temp_inx$geneB,"protein_type"]
      temp_inx$sumAbsScore <- sumAbsAB[sumAbsAB >= EdgeScoreThreshold]
      temp_inx$sumScaledScore <- sumScaledAB[sumAbsAB >= EdgeScoreThreshold]
      inx <- rbind(inx,temp_inx)
    }
    keysBA <- inxDB$key[inxDB$nodeA %in% temp_b & inxDB$nodeB %in% temp_a]
    sumAbsBA <- sapply(strsplit(keysBA,"_"),function(X)
      sum(abs(c(SCL[[b]][X[1],"score"],SCL[[a]][X[2],"score"]))))
    sumScaledBA <- sapply(strsplit(keysBA,"_"),function(X)
      sum(c(SCL[[b]][X[1],"scoreScaled"],SCL[[a]][X[2],"scoreScaled"])))
    if (any(sumAbsBA >= EdgeScoreThreshold)) {
      keysBA <- keysBA[sumAbsBA >= EdgeScoreThreshold]
      namesBA <- sapply(strsplit(keysBA,"_"),function(X)
        paste(paste(X[2],a,sep="_"),paste(X[1],b,sep="_"),sep="~"))
      temp_inx <- inxDB[keysBA,c("key","nodeB","nodeA")]
      rownames(temp_inx) <- temp_inx$key <- namesBA
      colnames(temp_inx)[2:3] <- c("geneA","geneB")
      temp_inx$nodeA <- paste(temp_inx$geneA,a,sep="_")
      temp_inx$nodeB <- paste(temp_inx$geneB,b,sep="_")
      temp_inx <- temp_inx[,c(1,4,5,2,3)]
      temp_inx$cellTypeA <- a
      temp_inx$cellTypeB <- b
      temp_inx$proteinTypeA <- geneInfo[temp_inx$geneA,"protein_type"]
      temp_inx$proteinTypeB <- geneInfo[temp_inx$geneB,"protein_type"]
      temp_inx$sumAbsScore <- sumAbsBA[sumAbsBA >= EdgeScoreThreshold]
      temp_inx$sumScaledScore <- sumScaledBA[sumAbsBA >= EdgeScoreThreshold]
      inx <- rbind(inx,temp_inx)
    }
    if (is.null(inx)) { return(NULL) }
    if (all(inx$sumAbsScore == 0)) {
      inx$sumAbsFactor <- 1
    } else {
      inx$sumAbsFactor <- cut(inx$sumAbsScore,10,labels=F)
    }
    inx$sumScaledFactor <- 6
    if (any(inx$sumScaledScore < 0)) {
      inx$sumScaledFactor[inx$sumScaledScore < 0] <-
        cut(inx$sumScaledScore[inx$sumScaledScore < 0],5,labels=F)
    }
    if (any(inx$sumScaledScore > 0)) {
      inx$sumScaledFactor[inx$sumScaledScore > 0] <-
        inx$sumScaledFactor[inx$sumScaledScore > 0] +
        cut(inx$sumScaledScore[inx$sumScaledScore > 0],5,labels=F)
    }
    inx$direction <- apply(apply(inx[,c("proteinTypeA","proteinTypeB")],1,function(X)
      c(grepl("Ligand",X),grepl("Receptor",X))),2,function(Y) {
        if (all(1:4 %in% which(Y))) {
          return("bothLR")
        } else if (all(c(1,4) %in% which(Y))) { # ligand from EC, with receptor in other cell type
          return("LtoR")
        } else if (all(c(2,3) %in% which(Y))) { # ligand from other cell type, receptor in EC
          return("RtoL")
        } else { return("tryECM") } # both ligands or both receptors
      })
    if (any(inx$direction == "tryECM")) {
      inx$direction[inx$direction == "tryECM"] <- apply(apply(inx[
        inx$direction == "tryECM",c("proteinTypeA","proteinTypeB")],1,function(X)
          c(grepl("ECM",X),grepl("Receptor",X))),2,function(Y) {
            if (all(1:4 %in% which(Y))) {
              return("bothER")
            } else if (all(c(1,4) %in% which(Y))) { # ligand from EC, with receptor in other cell type
              return("EtoR")
            } else if (all(c(2,3) %in% which(Y))) { # ligand from other cell type, receptor in EC
              return("RtoE")
            } else { return("tryEL") } # both ligands or both receptors
          })
    }
    if (any(inx$direction == "tryEL")) {
      inx$direction[inx$direction == "tryEL"] <- apply(apply(inx[
        inx$direction == "tryEL",c("proteinTypeA","proteinTypeB")],1,function(X)
          c(grepl("ECM",X),grepl("Ligand",X))),2,function(Y) {
            if (all(1:4 %in% which(Y))) {
              return("bothEL")
            } else if (all(c(1,4) %in% which(Y))) { # ligand from EC, with receptor in other cell type
              return("EtoL")
            } else if (all(c(2,3) %in% which(Y))) { # ligand from other cell type, receptor in EC
              return("LtoE")
            } else { return("None") } # both ligands or both receptors
          })
    }
    return(inx)
  },simplify=F)
  return(inxL)
}


#' Write node metadata for predicted cell-cell interaction networks
#'
#' This function uses the output of \code{\link{BuildInx}} and
#' \code{\link{ScoreClustGenes}} to write node metadata tables for all predicted
#' cell-cell interaction bipartite graphs.
#'
#' @param INX The output of \code{\link{BuildInx}}.
#' @param SCL The output of \code{\link{ScoreClustGenes}}.
#'
#' @export

BuildInxNode <- function(INX,SCL) {
  load(system.file("LigRecDB_RData/BaderCCIeditedbyBI.RData",package="CCInx"))
  message("Writing node metadata for cell-cell interaction networks...")
  return(pbapply::pbsapply(INX,function(X) {
    temp <- data.frame(node=c(unique(X$nodeA),unique(X$nodeB)),
                               gene=c(unique(X$geneA),unique(X$geneB)),stringsAsFactors=F)
    temp$cellType <- sapply(strsplit(temp$node,"_"),function(X) X[2])
    temp$proteinType <- geneInfo[temp$gene,"protein_type"]
    temp <- cbind(temp,
                  rbind(SCL[[unique(X$cellTypeA)]][unique(X$geneA),],
                        SCL[[unique(X$cellTypeB)]][unique(X$geneB),]))
    temp <- temp[!duplicated(temp$node),]
    rownames(temp) <- temp$node
    return(temp)
  },simplify=F))
}


#' Internal Fx - Prep for plotting
#'

prepBipartiteGraphSpread <- function(compName,INX,INXNODE,direction,topN,topPer) {
  if (!direction %in% c("LtoR","RtoL","RtoE","EtoR","LtoE","EtoL")) {
    stop(paste("direction must be one of:",
               paste(c("LtoR","RtoL","RtoE","EtoR","LtoE","EtoL"),collapse=";")))
  }
  co <- strsplit(compName,"~")[[1]]
  names(co) <- strsplit(direction,"to")[[1]]

  temp_both <- unique(grep("both",INX$direction,value=T))
  temp_both <- temp_both[grepl(names(co)[1],temp_both) & grepl(names(co)[2],temp_both)]
  edl <- INX$direction %in% c(direction,temp_both)
  if (!any(edl)) { return(NA) }
  INX <- INX[edl,]
  INXNODE <- INXNODE[unique(c(INX$nodeA,INX$nodeB)),]

  if (is.null(topN) & is.null(topPer)) {
    # warning(paste("Unless either 'topN' or 'topPer' is set, all edges",
    #               "will be shown in figure (which may be more than you want)."))
  } else if (is.null(topPer)) {
    INX <- INX[order(INX$sumAbsScore,decreasing=T),]
    INX <- head(INX,topN)
    INXNODE <- INXNODE[unique(c(INX$nodeA,INX$nodeB)),]
  } else if (is.null(topN)) {
    INX <- INX[order(INX$sumScaledScore,decreasing=T),]
    INX <- rbind(head(INX,topPer),tail(INX,topPer))
    INXNODE <- INXNODE[unique(c(INX$nodeA,INX$nodeB)),]
  } else {
    stop(paste("Only one of 'topN' (top # of edges by absolute DE score) or",
               "'topPer' (top DE edges up/down) can be set."))
  }

  temp <- list()
  z1 <- paste(co[1],names(co)[1],sep="_")
  z2 <- paste(co[2],names(co)[2],sep="_")
  temp[[z1]] <- INXNODE[INX$nodeA,
                        c("gene","score","scoreScaled","scaledFactor")]
  temp[[z1]] <- temp[[z1]][!duplicated(temp[[z1]]),]
  temp[[z1]]$ordered <- rank(
    temp[[z1]]$score,ties.method="random")/(nrow(temp[[z1]]) + 1)*100

  temp[[z2]] <- INXNODE[INX$nodeB,
                        c("gene","score","scoreScaled","scaledFactor")]
  temp[[z2]] <- temp[[z2]][!duplicated(temp[[z2]]),]
  temp[[z2]]$ordered <- rank(
    temp[[z2]]$score,ties.method="random")/(nrow(temp[[z2]]) + 1)*100

  temp$e <- INX[,c("nodeA","nodeB","sumAbsFactor","sumScaledFactor")]
  colnames(temp$e)[1:2] <- c(z1,z2)
  return(temp)
}


#' Internal Fx - Plotting
#'

plotBipartiteGraphSpread <- function(temp) {
  co <- sapply(strsplit(colnames(temp$e)[1:2],"_"),function(X) X[1])
  names(co) <- sapply(strsplit(colnames(temp$e)[1:2],"_"),function(X) X[2])
  co <- co[order(names(co))]
  z1 <- paste(co[1],names(co)[1],sep="_")
  z2 <- paste(co[2],names(co)[2],sep="_")
  names(co) <- sapply(names(co),function(X) switch(X,L="Ligands",R="Receptors",E="ECM"))
  temp_main <- paste(names(co)[1],"from",co[1],"to",names(co)[2],"from",co[2])

  plot(x=NULL,y=NULL,xlim=c(.1,8.9),ylim=c(0,100),yaxs="i",xaxs="i",
       xaxt="n",yaxt="n",xlab=NA,ylab=NA,main=temp_main,bty="n")
  points(c(rep(2,nrow(temp[[z1]])),rep(7,nrow(temp[[z2]]))),
         c(temp[[z1]]$ordered,temp[[z2]]$ordered),pch=21,cex=2,
         col=RColorBrewer::brewer.pal(11,"PRGn")[c(temp[[z1]]$scaledFactor,temp[[z2]]$scaledFactor)],
         bg=scales::alpha(RColorBrewer::brewer.pal(11,"PRGn")[c(temp[[z1]]$scaledFactor,temp[[z2]]$scaledFactor)],.5))

  temp_junk <- apply(temp$e,1,function(X)
    lines(x=c(2,7),y=c(temp[[z1]][X[[z1]],"ordered"],
                       temp[[z2]][X[[z2]],"ordered"]),
          lwd=seq(1,4,length.out=10)[as.integer(X["sumAbsFactor"])],
          col=scales::alpha(RColorBrewer::brewer.pal(11,"PRGn")[as.integer(X["sumScaledFactor"])],
                    seq(.6,1,length.out=10)[as.integer(X["sumAbsFactor"])])))

  text(x=2,y=temp[[z1]]$ordered,labels=temp[[z1]]$gene,pos=2,col="black")
  text(x=7,y=temp[[z2]]$ordered,labels=temp[[z2]]$gene,pos=4,col="black")

  mtext(text=co,side=1,at=c(2,7),line=.5,font=2)
  mtext(text=c("Negative","Differential Gene Expression","Positive"),
        side=2,at=c(10,50,90),line=1.2,font=2,
        col=c(RColorBrewer::brewer.pal(11,"PRGn")[1],"black",RColorBrewer::brewer.pal(11,"PRGn")[11]))
  barplot(rep(-.1,100),space=0,border=NA,horiz=T,xaxt="n",add=T,
          col=RColorBrewer::brewer.pal(11,"PRGn")[cut(1:100,11,labels=F)])
  barplot(rep(.1,100),space=0,border=NA,horiz=T,xaxt="n",add=T,
          col=RColorBrewer::brewer.pal(11,"PRGn")[cut(1:100,11,labels=F)])
}


#' Plot bipartite graphs of predicted cell-cell interaction networks.
#'
#' This function generates figures displaying bipartite graphs for predicted
#' cell-cell interaction networks between pairs of cell types.
#'
#' @param INX The output of \code{\link{BuildInx}}.
#' @param INXNODE The output of \code{\link{BuildInxNode}}.
#' @param comparisons Optional, default is to plot all pairs of cell types. A
#'   character vector representing the names of \code{inx} to plot.
#' @param directions Optional, default is to plot all combinations of ligands,
#'   receptors, and ECM components. A character vector representing sets of
#'   interactions to consider. Must be a subset of
#'   \code{c('LtoR','RtoL','RtoE','EtoR','LtoE','EtoL')}. 'L' refers to ligands,
#'   'R' to receptors, and 'E' to ECM components. The order is the same as that
#'   in \code{comparisons} (the names of \code{inx}), such that 'RtoL' for
#'   \code{inx[["A~B"]]} will plot the receptors from cell type A and their
#'   ligands from cell type B.
#' @param topN Optional. Number of edges to include in the plot, ranked by
#'   absolute edge score. Only one of \code{topN} and \code{topPer} can be
#'   specified.
#' @param topPer Optional. Number of most differentially expressed edges to
#'   include from each side of the comparison. Only one of \code{topN} and
#'   \code{topPer} can be specified.
#' @param imageFileType Default="pdf". The file format for saved figures. One of
#'   \code{"pdf"} (generated with \code{\link[grDevices]{cairo_pdf}}),
#'   \code{"eps"} (generated with \code{\link[grDevices]{cairo_ps}}),
#'   \code{"tiff"} (generated with \code{\link[grDevices]{tiff}}), or
#'   \code{"png"} (generated with \code{\link[grDevices]{png}}).
#' @param imageFilePath Default is working directory. Character vector
#'   indicating where to save image files.
#'
#' @export

PlotCCInx <- function(INX,INXNODE,
                      comparisons,
                      directions=c('LtoR','RtoL','RtoE','EtoR','LtoE','EtoL'),
                      topN,topPer,
                      imageFileType="pdf",
                      imageFilePath="./") {
  if (missing(comparisons)) { comparisons <- names(INX) }
  if (missing(topN)) { topN <- NULL }
  if (missing(topPer)) { topPer <- NULL }
  if (!imageFileType %in% c("pdf","eps","png","tiff")) {
    warning('imageFileType must be one of c("pdf","eps","png","tiff"). Setting to "pdf".')
    imageFileType <- "pdf"
  }

  garbage <- pbapply::pbsapply(comparisons,function(l) {
    if (is.null(INX[[l]])) { return(NULL) }
    tempL <- sapply(directions,function(tf)
      prepBipartiteGraphSpread(compName=l,
                               INX=INX[[l]],
                               INXNODE=INXNODE[[l]],
                               direction=tf,
                               topN=topN,
                               topPer=topPer),
      simplify=F)
    tempL <- tempL[sapply(tempL,typeof) == "list"]
    if (length(tempL) == 0) { return(NULL) }
    # temp_figH <- 0.20 * max(sapply(tempL,function(X) sapply(X[1:2],nrow)))
    # if (temp_figH < 5) { temp_figH <- 5 }
    # ^ If you want all figures from a comparison to be the same size.
    for (tf in names(tempL)) {
      temp_figH <- 0.20 * max(sapply(tempL[[tf]][1:2],nrow))
      if (temp_figH < 5) { temp_figH <- 5 }
      if (is.null(topN) & is.null(topPer)) {
        temp_FP <- paste0(imageFilePath,sub("~","_",l),"_",tf,"_all.",imageFileType)
      } else if (is.null(topPer)) {
        temp_FP <- paste0(imageFilePath,sub("~","_",l),"_",tf,"_topN",topN,".",imageFileType)
      } else if (is.null(topN)) {
        temp_FP <- paste0(imageFilePath,sub("~","_",l),"_",tf,"_topPer",topPer,".",imageFileType)
      }
      switch(imageFileType,
             "pdf"=grDevices::cairo_pdf(temp_FP,height=temp_figH,width=6,fallback_resolution=300),
             "eps"=grDevices::cairo_ps(temp_FP,height=temp_figH,width=6,fallback_resolution=300),
             "tiff"=grDevices::tiff(temp_FP,height=temp_figH,width=6,units="in",res=300),
             "png"=grDevices::png(temp_FP,height=temp_figH,width=6,units="in",res=300))
      par(mar=c(2,3,2,1),mgp=c(1,1,0))
      plotBipartiteGraphSpread(tempL[[tf]])
      dev.off()
    }
  })
}
