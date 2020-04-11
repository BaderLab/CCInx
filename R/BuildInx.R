#' Build cluster-wise list of gene expression statistics from scRNAseq data
#'
#' This function takes a \code{Seurat} or \code{SingleCellExperiment} object and
#' builds a list of dataframes containing gene expression statistics for all
#' genes of each cluster. This can be used as the input to
#' \code{\link{BuildCCInx}} for generating cell-cell interaction predictions
#' between cell-type clusters.
#'
#' @param inD The input dataset. An object of class \code{\link[Seurat]{seurat}}
#'   or \code{\link[SingleCellExperiment]{SingleCellExperiment}}. Other data
#'   classes are not currently supported.
#'   \href{https://github.com/BaderLab/scClustViz/issues}{Please submit requests
#'   for other data objects here!}
#' @param cl a factor where each value is the cluster assignment for a cell
#'   (column) in the input gene expression matrix.
#' @param assayType Default = "" (for Seurat v1/2). A length-one character
#'   vector representing the assay slot in which the expression data is stored
#'   in the input object. This is not required for Seurat v1 or v2 objects. See
#'   \code{\link[scClustViz]{getExpr}} for details.
#' @param assaySlot An optional length-one character vector representing
#'   the slot of the Seurat v3 \code{\link[Seurat]{Assay}} object to use. Not
#'   used for other single-cell data objects. The default is to use the
#'   normalized data in the "data" slot, but you can also use the
#'   \code{\link[Seurat]{SCTransform}}-corrected counts by setting
#'   \code{assayType = "SCT"} and \code{assaySlot = "counts"}. This is
#'   recommended, as it will speed up differential expression
#'   calculations. See \code{\link{getExpr}} for details.
#' @param exponent Default = 2. A length-one numeric vector representing the
#'   base of the log-normalized gene expression data to be processed. Generally
#'   gene expression data is transformed into log2 space when normalizing (set
#'   this to 2), though \code{Seurat} uses the natural log (set this to exp(1)).
#' @param pseudocount Default = 1. A length-one numeric vector representing the
#'   pseudocount added to all log-normalized values in your input data. Most
#'   methods use a pseudocount of 1 to eliminate log(0) errors.
#'
#' @seealso \code{\link[scClustViz]{CalcCGS}}
#'
#' @references Mean gene expression calculations
#'   Innes BT and Bader GD. scClustViz – Single-cell RNAseq cluster
#'   assessment and visualization [version 2; peer review: 2 approved].
#'   F1000Research 2019, 7:1522 (\url{https://doi.org/10.12688/f1000research.16198.2})
#'
#' @export

BuildGeneStatList <- function(inD,
                              cl,
                              assayType="",
                              assaySlot="",
                              exponent=2,
                              pseudocount=1) {
  if (!require(scClustViz)) {
    stop(paste("scClustViz is required for this function. Install from github:",
               "    devtools::install_github('Baderlab/scClustViz')",sep="\n"))
  }
  if (!is(inD)[1] %in% methods::findMethodSignatures(scClustViz::getExpr)) {
    stop(paste(
      paste0("Input data object must be one of: ",
             paste(methods::findMethodSignatures(scClustViz::getExpr),collapse=", "),
             "."),
      paste("Other input objects are not supported at this time,",
            "but please let me know what object class"),
      paste("you'd like supported at",
            "https://github.com/BaderLab/scClustViz/issues, thanks!"),
      sep="\n  "))
  }
  if (is.null(colnames(getExpr(inD,assayType,assaySlot))) |
      is.null(rownames(getExpr(inD,assayType,assaySlot)))) {
    stop("Gene expression matrix returned by 'getExpr(inD,assayType,assaySlot)' is missing col/rownames.")
  }
  if (length(cl) != ncol(scClustViz::getExpr(inD,assayType,assaySlot))) {
    stop(paste("cl must be a factor where each value is the cluster assignment",
               "for a cell (column) in the input gene expression matrix.",
               sep="\n  "))
  }
  if (is.character(cl) | is.numeric(cl)) {
    cl <- as.factor(cl)
  }
  if (!all(names(cl) == colnames(scClustViz::getExpr(inD,assayType,assaySlot))) |
      is.null(names(cl))) {
    names(cl) <- colnames(scClustViz::getExpr(inD,assayType,assaySlot))
  }
  if (any(grepl("_",levels(cl)))) {
    stop("Cluster names cannot contain '_' due to internal naming conventions.")
  }
  if (any(grepl("~",levels(cl)))) {
    stop("Cluster names cannot contain '~' due to internal naming conventions.")
  }

  temp <- scClustViz:::fx_calcCGS(nge=scClustViz::getExpr(inD,
                                                          assayType,
                                                          assaySlot),
                                  cl=cl,
                                  exponent=2,
                                  pseudocount=1)
  return(
    sapply(temp,function(X) {
      X <- X[X$DR > 0,]
      names(X)[names(X) == "DR"] <- "DetectRate"
      names(X)[names(X) == "MDGE"] <- "MeanDetectGeneExpr"
      names(X)[names(X) == "MGE"] <- "MeanNormGeneExpr"
      return(X)
    },simplify=F)
  )
}





#' Check input and score and scale DE stats.
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

CalcDiffExprScaled <- function(gdb,DEmagn,DEstat) {
  if (any( is.na(gdb[[DEmagn]]) )) {
    stop(paste("This function doesn't tolerate missing",
               DEmagn,"values."))
  }
  if (any( is.na(gdb[[DEstat]]) )) {
    stop(paste("This function doesn't tolerate missing",
               DEstat,"values."))
  }
  return(sapply(gdb[[DEmagn]],function(X) {
    if (X == Inf) {
      1.1
    } else if (X == -Inf) {
      -1.1
    } else {
      X / switch(as.character(X >= 0),
                 "TRUE"=max(gdb[[DEmagn]][!is.infinite(gdb[[DEmagn]])]),
                 "FALSE"=min(gdb[[DEmagn]][!is.infinite(gdb[[DEmagn]])]) * -1)
    }
  }))
}


#' Check input and score and scale gene expression.
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
#' @param expr A character vector of length 1 representing the variable name
#'   in the GeneStatList data frames carrying information on the magnitude and
#'   direction of the change of expression for the node (gene) in each cell
#'   type. This is generally a signed logFC or gene expression ratio.

CalcExprScaled <- function(gdb,expr) {
  if (any( is.na(gdb[[expr]]) )) {
    stop(paste("This function doesn't tolerate missing",
               expr,"values."))
  }
  return( ( gdb[[expr]] - min(gdb[[expr]]) ) /
            max( gdb[[expr]] - min(gdb[[expr]]) ) )
}

#' Build cell-cell interaction predictions between cell types
#'
#' This function takes a list of gene statistics per cluster to predict
#' cell-cell interactions between each cell-type (cluster). If the
#' \code{GeneStatistic} argument is provided, this function will assume the gene
#' statistics represent differential expression between experimental conditions,
#' and will weight the predicted interactions accordingly. Otherwise,
#' predictions will be weighted by expression magnitude per cell type. The
#' output of this function can be explored interactively with
#' \code{\link{ViewCCInx}}, or static figures can be generated with
#' \code{\link{PlotCCInx}}.
#'
#' @param GeneStatList A named list of dataframes. Each list element should
#'   represent a cell type / cluster to be included in the interaction network,
#'   and should be named accordingly. List elements should contain data frames
#'   where each row is a gene with official gene symbols as row names. Variables
#'   should be appropriately named statistics or annotations to be included in
#'   the resulting node metadata. Variable names should be consistent between
#'   list elements. The function \code{\link{BuildGeneStatList}} can be used to
#'   generate this list from \code{Seurat} or \code{SingleCellExperiment}
#'   objects when generating predictions not involving differential gene
#'   expression.
#' @param GeneMagnitude Default = "MeanNormGeneExpr". A character vector of length 1
#'   representing the variable name in the GeneStatList data frames carrying
#'   information on the magnitude (and direction of the change) of expression
#'   for the node (gene) in each cell type. This is either a measure of
#'   expression (generally mean expression or detection rate) or a measure of
#'   change (signed log expression ratio a.k.a. logFC). Default assumes
#'   \code{GeneStatList} is output from \code{\link{BuildGeneStatList}}, and
#'   uses mean normalized gene expression to weight nodes and edges.
#' @param GeneStatistic Optional. A character vector of length 1 representing
#'   the variable name in the GeneStatList data frames carrying information on
#'   the statistical significance of expression change. This is generally a
#'   corrected p-value.
#' @param Species Default='hsapiens'. The species of the source data. One of
#'   'hsapiens' or 'mmusculus'. Note that the ligand-receptor database was built
#'   for human, and the mouse version is generated by homology mapping (only
#'   using uniquely mapped homologues).
#'
#' @export

BuildCCInx <- function(GeneStatList,
                       GeneMagnitude="MeanNormGeneExpr",
                       GeneStatistic,
                       Species="hsapiens") {
  if (length(names(GeneStatList)) != length(GeneStatList)) {
    stop("GeneStatList must be a named list where names are cell types.")
  }
  if (any(duplicated(names(GeneStatList)))) {
    stop("GeneStatList names must be unique.")
  }
  if (any(grepl("_",names(GeneStatList)))) {
    stop("GeneStatList names cannot contain '_' due to internal naming conventions.")
  }
  if (any(grepl("~",names(GeneStatList)))) {
    stop("GeneStatList names cannot contain '~' due to internal naming conventions.")
  }
  message("Scaling node weights per cell type...")
  if (missing(GeneStatistic)) {
    temp_scaled <- pbapply::pbsapply(X=GeneStatList,
                                     FUN=CalcExprScaled,
                                     expr=GeneMagnitude,
                                     simplify=F)
  } else {
    temp_scaled <- pbapply::pbsapply(X=GeneStatList,
                                     FUN=CalcDiffExprScaled,
                                     DEmagn=GeneMagnitude,
                                     DEstat=GeneStatistic,
                                     simplify=F)
  }
  inx <- list()
  message("Building node metadata...")
  temp_cellNames <- names(GeneStatList)
  inx$nodes <- do.call(rbind,GeneStatList)
  temp_gene <- unlist(lapply(GeneStatList,rownames),use.names=F)
  temp_cellType <- unlist(mapply(function(N,X) rep(N,X),
                                 N=temp_cellNames,
                                 X=sapply(GeneStatList,nrow),
                                 SIMPLIFY=F),
                          use.names=F)

  switch(Species,
         hsapiens=load(system.file("LigRecDB_RData/BaderCCIeditedbyBI_human.RData",
                                   package="CCInx")),
         mmusculus=load(system.file("LigRecDB_RData/BaderCCIeditedbyBI_mouse.RData",
                                    package="CCInx")),
         MillerKaplan=load(system.file("LigRecDB_RData/MillerKaplan_mouse.RData",
                                       package="CCInx")),
         FANTOM5=load(system.file("LigRecDB_RData/FANTOM5_human.RData",
                                  package="CCInx")),
         stop("Species must be one of 'hsapiens' or 'mmusculus'."))
  if (sum(rownames(geneInfo) %in% temp_gene) < 20) {
    warning(paste("Less than 20 genes from GeneStatList were detected in the CCInx database.",
                  "  Please ensure that you've set the Species argument correctly.",
                  "  Rownames of each entry in GeneStatList must be official gene symbols.",
                  sep="\n"))
  }
  temp_proteinType <- geneInfo[temp_gene,"protein_type"]
  inx$nodes <- cbind(data.frame(node=paste(temp_gene,temp_cellType,sep="_"),
                                gene=temp_gene,
                                cellType=temp_cellType,
                                proteinType=temp_proteinType,
                                nodeWeight=unlist(temp_scaled),
                                stringsAsFactors=F),
                     inx$nodes)
  rownames(inx$nodes) <- inx$nodes$node
  inx$nodes <- inx$nodes[!is.na(inx$nodes$proteinType),]

  tempCN <- c()
  for (a in temp_cellNames) {
    for (b in temp_cellNames) {
      temp <- paste(sort(c(a,b)),collapse="~")
      if (!temp %in% tempCN) {
        tempCN <- append(tempCN,temp)
      }
    }
  }
  rm(a,b)
  tempComp <- strsplit(tempCN,"~")
  names(tempComp) <- tempCN

  message("Building edge list...")
  inx$edges <- pbapply::pbsapply(tempComp,function(Z) {
    a <- Z[1]; b <- Z[2]
    if (sum(inx$nodes$cellType == a) < 1 |
        sum(inx$nodes$cellType == b) < 1) {
      return(NULL)
    }

    keysAB <- inxDB$key[inxDB$nodeA %in% inx$nodes$gene[inx$nodes$cellType == a] &
                          inxDB$nodeB %in% inx$nodes$gene[inx$nodes$cellType == b]]
    if (length(keysAB) < 1) {
      edgesAB <- data.frame(row.names=c("key","nodeA","nodeB"))
    } else {
      edgesAB <- data.frame(
        sapply(strsplit(keysAB,"_"),function(X)
          paste(paste(X[1],a,sep="_"),paste(X[2],b,sep="_"),sep="~")),
        t(sapply(strsplit(keysAB,"_"),function(X)
          c(paste(X[1],a,sep="_"),paste(X[2],b,sep="_")))),
        stringsAsFactors=F)
      colnames(edgesAB) <- c("key","nodeA","nodeB")
      rownames(edgesAB) <- edgesAB$key
    }

    keysBA <- inxDB$key[inxDB$nodeA %in% inx$nodes$gene[inx$nodes$cellType == b] &
                          inxDB$nodeB %in% inx$nodes$gene[inx$nodes$cellType == a]]
    if (length(keysBA) < 1) {
      edgesBA <- data.frame(row.names=c("key","nodeA","nodeB"))
    } else {
      edgesBA <- data.frame(
      sapply(strsplit(keysBA,"_"),function(X)
        paste(paste(X[2],a,sep="_"),paste(X[1],b,sep="_"),sep="~")),
      t(sapply(strsplit(keysBA,"_"),function(X)
        c(paste(X[2],a,sep="_"),paste(X[1],b,sep="_")))),
      stringsAsFactors=F)
    colnames(edgesBA) <- c("key","nodeA","nodeB")
    rownames(edgesBA) <- edgesBA$key
    }

    temp <- rbind(edgesAB,edgesBA)
    if (nrow(temp) < 1) {
      return(NULL)
    } else {
      return(temp)
    }
  },simplify=F)
  inx$edges <- do.call(rbind,inx$edges)
  rownames(inx$edges) <- inx$edges$key
  inx$edges <- inx$edges[,2:3]
  inx$edges$edgeWeight <- rowMeans(cbind(inx$nodes[inx$edges$nodeA,"nodeWeight"],
                                         inx$nodes[inx$edges$nodeB,"nodeWeight"]))

  # inx$nodes <- inx$nodes[inx$nodes$node %in% inx$edges$nodeA | inx$nodes$node %in% inx$edges$nodeB,]
  # Orphan nodes should be left in for clarity - the genes were in the database, just without interacting partners.
  # They get filtered from the viewer in FilterInx_step1 anyway.

  attr(inx,"GeneMagnitude") <- GeneMagnitude
  if (!missing(GeneStatistic)) {
    attr(inx,"GeneStatistic") <- GeneStatistic
  }

  return(inx)
}

