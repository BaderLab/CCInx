# Here put fx for cleaning up the CCI data from web, homology mapping fx.


load(system.file("LigRecDB_RData/BaderCCIeditedbyBI.RData",package="CCInx"))
geneInfo_pf <- geneInfo
inxDB_pf <- inxDB
rm(geneInfo,inxDB)

require(biomaRt)
martM <- useMart("ensembl","mmusculus_gene_ensembl")
martH <- useMart("ensembl","hsapiens_gene_ensembl")
h2m <- getLDS(attributes="hgnc_symbol",mart=martH,
              filters="hgnc_symbol",values=geneInfo$hgnc_symbol,
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
