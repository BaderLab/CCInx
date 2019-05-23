# Here put fx for cleaning up the CCI data from web, homology mapping fx.


#### Homology map to make mouse db ####
load(system.file("LigRecDB_RData/BaderCCIeditedbyBI.RData",package="CCInx"))
geneInfo_pf <- geneInfo
inxDB_pf <- inxDB
rm(geneInfo,inxDB)

require(biomaRt)
martM <- useMart("ensembl","mmusculus_gene_ensembl")
martH <- useMart("ensembl","hsapiens_gene_ensembl")
h2m <- getLDS(attributes="hgnc_symbol",mart=martH,
              filters="hgnc_symbol",values=geneInfo_pf$hgnc_symbol,
              attributesL="mgi_symbol",martL=martM)

## using only uniquely mapping homologs ##
temp <- unique(h2m$HGNC.symbol[duplicated(h2m$HGNC.symbol)])
h2m <- h2m[!h2m$HGNC.symbol %in% temp,]
temp <- unique(h2m$MGI.symbol[duplicated(h2m$MGI.symbol)])
h2m <- h2m[!h2m$MGI.symbol %in% temp,]
#### ^ THIS METHOD SUCKS ####

geneInfo <- geneInfo_pf[geneInfo_pf$hgnc_symbol %in% h2m$HGNC.symbol,]
geneInfo$mgi_symbol <- sapply(geneInfo$hgnc_symbol,function(X)
  h2m$MGI.symbol[h2m$HGNC.symbol == X])
rownames(geneInfo) <- geneInfo$mgi_symbol

inxDB <- inxDB_pf[inxDB_pf$nodeA %in% geneInfo$hgnc_symbol & inxDB_pf$nodeB %in% geneInfo$hgnc_symbol,]
inxDB$nodeA <- sapply(inxDB$nodeA,function(X) geneInfo$mgi_symbol[geneInfo$hgnc_symbol == X])
inxDB$nodeB <- sapply(inxDB$nodeB,function(X) geneInfo$mgi_symbol[geneInfo$hgnc_symbol == X])
rownames(inxDB) <- inxDB$key <- paste(inxDB$nodeA,inxDB$nodeB,sep="_")

save(inxDB,geneInfo,file="inst/LigRecDB_RData/BaderCCIeditedbyBI_mouse.RData")



#### Format Miller/Kaplan db for use ####
mkdb <- read.table("../CCCnetResources/MillerKaplan_LRdatabase.txt",
                   sep="\t",header=T,as.is=T)

# apparently despite the all caps symbols, these are mouse genes
# converting them to MGI style sentence case gene names
mkdb$Ligand <- sapply(strsplit(tolower(mkdb$Ligand),""),function(X)
  paste(c(toupper(X[1]),X[-1]),collapse=""))
mkdb$Receptor <- sapply(strsplit(tolower(mkdb$Receptor),""),function(X)
  paste(c(toupper(X[1]),X[-1]),collapse=""))

temp_ligand <- data.frame(symbol=unique(mkdb$Ligand),
                          protein_type="Ligand",
                          stringsAsFactors=F)
temp_receptor <- data.frame(symbol=unique(mkdb$Receptor),
                            protein_type="Receptor",
                            stringsAsFactors=F)
temp_ligand$protein_type[temp_ligand$symbol %in% temp_receptor$symbol] <- "Receptor/Ligand"
temp_receptor$protein_type[temp_receptor$symbol %in% temp_ligand$symbol] <- "Receptor/Ligand"

geneInfo <- rbind(temp_ligand,temp_receptor)
geneInfo <- geneInfo[!duplicated(geneInfo$symbol),]
rownames(geneInfo) <- geneInfo$symbol
geneInfo <- geneInfo[geneInfo$symbol != "None",]

inxDB <- data.frame(key=paste(mkdb$Ligand,mkdb$Receptor,sep="_"),
                    nodeA=mkdb$Ligand,
                    nodeB=mkdb$Receptor,
                    stringsAsFactors=F)
rownames(inxDB) <- inxDB$key
inxDB <- inxDB[!grepl("None",inxDB$key),]

save(geneInfo,inxDB,file="inst/LigRecDB_RData/MillerKaplan_mouse.RData")

