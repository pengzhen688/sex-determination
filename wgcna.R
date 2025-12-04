# ===============================
# WGCNA on gene expression
# ===============================

# --------- step 1. Environmental and Packaging Inspection ---------
options(stringsAsFactors = FALSE)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
needed <- c("WGCNA", "DESeq2", "edgeR", "pheatmap")
for (pkg in needed) {
  if (!suppressWarnings(require(pkg, character.only = TRUE))) {
    if (pkg %in% c("DESeq2", "edgeR")) BiocManager::install(pkg)
    else install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
library(WGCNA); library(DESeq2); library(edgeR); library(pheatmap);library(readxl);library(openxlsx)
allowWGCNAThreads()
enableWGCNAThreads(nThreads = 12)

# --------- step 2. Read and process files ---------
fn <- "count.txt"
expr_raw <- read.table(fn, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
cat("Read expression matrix:", dim(expr_raw)[1], "genes x", dim(expr_raw)[2], "samples\n")
print(colnames(expr_raw))

# --------- step 3. Constructing sample information（Sex） ---------
sampleInfo <- data.frame(
  sample = colnames(expr_raw),
  Sex = factor(c(rep("F", 6), rep("M", 6)), levels = c("F", "M"))
)
rownames(sampleInfo) <- sampleInfo$sample

# --------- step 4. Low expression gene filtering ---------
dge <- DGEList(counts = expr_raw)
# Retain genes when CPM >= 1 in at least 3 samples (threshold adjustable)
keep <- rowSums(cpm(dge) >= 1) >= 3
expr_filt <- expr_raw[keep, ]
cat("After low-count filter:", nrow(expr_filt), "genes retained\n")
write.table(expr_filt, file = "expr_filtered_counts.txt", sep = "\t", quote = FALSE)

# --------- step 5. Normalization / Transformation ---------
dds <- DESeqDataSetFromMatrix(countData = expr_filt, colData = sampleInfo, design = ~ Sex)
dds <- estimateSizeFactors(dds) 
# Choose a transformation: vst is suitable for most cases; if ngenes is very small, rlog can be used.
nGenes <- nrow(expr_filt)
cat("Number of genes after filter:", nGenes, "\n")
if (nGenes < 2000) {
  cat("Using rlog (recommended for small gene sets)\n")
  rld <- rlog(dds, blind = TRUE)
  expr_vst <- assay(rld)
} else {
  cat("Using vst\n")
  vst_res <- varianceStabilizingTransformation(dds, blind = TRUE)
  expr_vst <- assay(vst_res)
}
# expr_vst: rows = genes, cols = samples
write.table(round(expr_vst,4), file = "expr_vst_matrix.txt", sep = "\t", quote = FALSE)

# --------- step 6. Data format conversion ---------
datExpr <- t(expr_vst)
datExpr <- as.data.frame(datExpr)
datExpr[] <- lapply(datExpr, as.numeric)
rownames(datExpr) <- colnames(expr_vst)

# --------- step 7. Data quality check ---------
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) cat("Removing bad genes:", sum(!gsg$goodGenes), "\n")
  if (sum(!gsg$goodSamples) > 0) cat("Removing bad samples:", sum(!gsg$goodSamples), "\n")
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Sample clustering diagram (for outlier detection)
sampleTree <- hclust(dist(datExpr), method = "average")
pdf("sample_clustering.pdf", width = 10, height = 7)
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.2)
dev.off()

# --------- step 8. Choosing soft-threshold power ---------
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# picture scale-free fit
pdf("scale_free_fit.pdf", width = 14, height = 7)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2",
     type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=cex1, col="red")
abline(h=0.7, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()
#Based on the above drawings, we choose power = 5
softPower <- 5

# --------- step 9. Build the network and identify blockwise modules.（blockwiseModules） ---------
cor <- WGCNA::cor
net <- blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = "signed",
  minModuleSize = 30,
  maxBlockSize = 20000,
  deepSplit = 3,
  reassignThreshold = 0,
  mergeCutHeight = 0.15,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  verbose = 3
)

table(net$colors)
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
table(moduleColors)
# picture gene dendrogram and module colors
pdf("dendro_modules.pdf", width = 15, height = 8)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, main = "Gene dendrogram and module colors")
dev.off()

# --------- step 9. Merge similar modules ---------
# Calculate the module feature vector
MEList <- moduleEigengenes(datExpr, colors = moduleColors)
MEs <- MEList$eigengenes

# Calculate the distance between modules and cluster them.
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

# Visualization of the correlation tree between modules
pdf("Clustering_of_module_eigengenes.pdf", width = 7, height = 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
MEDissThres <- 0.3
abline(h = MEDissThres, col = "red")  
dev.off()

# Automatically merge modules with similar expression patterns
merge <- mergeCloseModules(datExpr, moduleColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors
table(mergedColors)
mergedMEs <- merge$newMEs
color<-unique(mergedColors)
for (i  in 1:length(color)) {
  y=t(assign(paste(color[i],"expr",sep = "."),datExpr[mergedColors==color[i]]))
  write.csv(y,paste('hebing',color[i],"csv",sep = "."),quote = F)
}

# Gene heatmaps include comparisons before and after merging.
pdf("Merged_Module_Dendrogram.pdf", width = 10, height = 6)
plotDendroAndColors(
  net$dendrograms[[1]],
  cbind(moduleColors, mergedColors),
  c("Original Modules", "Merged Modules"),
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05
)
dev.off()

# Reassign new values ​​to the module
moduleColors <- mergedColors
MEs <- mergedMEs
colorOrder = c("grey", standardColors(50)) 
moduleLabels = match(moduleColors, colorOrder)-1 
names(moduleColors) <- names(datExpr)
table(moduleColors)

# Module reordering and renaming
moduleSizes <- table(moduleColors)
moduleSizes <- moduleSizes[names(moduleSizes) != "grey"]
moduleSizes <- sort(moduleSizes, decreasing = TRUE)  # 计算每个模块的大小
orderedColors <- names(moduleSizes)  # 按模块大小排序的颜色顺序
color_to_newID <- setNames(seq_along(orderedColors), orderedColors) # 按模块大小生成新的连续编号
moduleLabels <- as.numeric(factor(moduleColors, levels = orderedColors)) # 为每个基因赋予新的模块编号
moduleColors <- labels2colors(moduleLabels) # 根据新编号生成新的颜色
MEs <- MEs[paste0("ME", orderedColors)]   # 重新排列模块特征向量（MEs）列顺序
colnames(MEs) <- paste0("ME", seq_along(orderedColors))   # 重新命名 MEs 列为新的编号
module_info <- data.frame(newID = seq_along(orderedColors),color = orderedColors,MEname = paste0("ME", seq_along(orderedColors)))
print(head(module_info))

# --------- step 10. Relating modules to external traits ---------
stopifnot(all(colnames(expr_vst) %in% rownames(sampleInfo)))  
traitData <- data.frame(
  Male   = ifelse(sampleInfo[colnames(expr_vst), "Sex"] == "M", 1, 0),
  Female = ifelse(sampleInfo[colnames(expr_vst), "Sex"] == "F", 1, 0)
)
rownames(traitData) <- colnames(expr_vst) 

moduleTraitCor <- cor(MEs, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))
write.csv(moduleTraitCor, file = "module_trait_correlation.csv")
write.csv(moduleTraitPvalue, file = "module_trait_pvalues.csv")

pdf("module_trait_heatmap.pdf", width = 8, height = 10)
textMatrix <- paste(signif(moduleTraitCor, 2), "(", signif(moduleTraitPvalue, 2), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(traitData),    # 性状
  yLabels = colnames(MEs),          # 模块名
  ySymbols = colnames(MEs),         # 模块名
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  main = "Module-trait relationships (correlation (p-value))"
)
dev.off()

# --------- step 11. Calculate the correlation between genes and modules (MM) and the correlation between genes and phenotypes (GS) ---------
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples = nrow(datExpr)))

geneTraitSignificance <- as.data.frame(cor(datExpr, traitData$Male, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples = nrow(datExpr)))

write.csv(geneModuleMembership, "geneModuleMembership.csv")
write.csv(geneTraitSignificance, "geneTraitSignificance.csv")

# --------- step 12. Identify hub genes in significant modules and export the information for picturing ---------
sigModules <- rownames(moduleTraitCor)[apply(abs(moduleTraitCor) >= 0.45 & moduleTraitPvalue <= 0.1, 1, any)]
cat("Significant modules:\n")
print(sigModules)

for (mod in sigModules) {
  modGenes <- moduleColors == substring(mod, 3) 
  moduleGenes <- names(datExpr)[modGenes]
  kMEvalues <- geneModuleMembership[moduleGenes, mod,drop=FALSE]
  GSvalues <- geneTraitSignificance[moduleGenes, 1,drop=FALSE]
  pdf(paste0("MM_vs_GS_", mod, ".pdf"), width = 6, height = 6)
  verboseScatterplot(
    abs(kMEvalues[,1]), abs(GSvalues[,1]),
    xlab = paste("Module Membership in", mod, "module"),
    ylab = "Gene significance for phenotype",
    main = paste("Module membership vs. gene significance\n"),
    abline = TRUE,
    pch = 21,
    bg = moduleColors, 
    col = "black"
  )
  dev.off()
  
  hubData <- data.frame(
    Gene = rownames(kMEvalues),
    ModuleMembership = kMEvalues[, 1],
    GeneSignificance = GSvalues[rownames(kMEvalues), 1]
  )
  hubData <- hubData[hubData$ModuleMembership > 0.8 & hubData$GeneSignificance > 0.5, ]
  if (nrow(hubData) == 0) {
    cat(paste0("Module ", mod, ": No genes meet MM>0.8 & GS>0.5 criteria\n"))
    next
  }
  hubData <- hubData[order(abs(hubData$ModuleMembership), decreasing = TRUE), ]
  
  write.table(hubData$Gene, file = paste0("HubGenes_", mod, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.csv(hubData, file = paste0("HubGene_Details_", mod, ".csv"), row.names = FALSE)
  
  cat(paste0("Module ", mod, ": ", nrow(hubData), " hub genes saved (MM>0.8 & GS>0.5)\n"))
  
  inModule <- colnames(datExpr) %in% hubData$Gene
  modTOM <- TOM[inModule, inModule]
  dimnames(modTOM) <- list(hubData$Gene, hubData$Gene)
  IMConn <- softConnectivity(datExpr[, hubData$Gene])
  names(IMConn) <- hubData$Gene
  topN <- min(100, length(hubData$Gene))
  topGenes <- names(sort(IMConn, decreasing = TRUE))[1:topN]

  filterTOM <- modTOM[topGenes, topGenes]
  vis <- exportNetworkToVisANT(filterTOM, file = paste("visANTinput-", mod, ".txt", sep = ""),
                               weighted = TRUE, threshold = 0)
  
  cyt <- exportNetworkToCytoscape(filterTOM,
                                  edgeFile = paste("CytoscapeInput-edges-", mod, ".txt", sep = ""),
                                  nodeFile = paste("CytoscapeInput-nodes-", mod, ".txt", sep = ""),
                                  weighted = TRUE, threshold = 0.2,
                                  nodeNames = topGenes, nodeAttr = rep(mod, length(topGenes)))
}


# --------- step 13. Overlap analysis between modules and known gene sets ---------
DEGs_M1vsF1 <- read_excel("M1vsF1_deg.xlsx")
UP_DEGs_M1vsF1 <- as.character(DEGs_M1vsF1$gene.id[DEGs_M1vsF1$pvalue < 0.01 & DEGs_M1vsF1$log2FoldChange > 1])
DOWN_DEGs_M1vsF1 <- as.character(DEGs_M1vsF1$gene.id[DEGs_M1vsF1$pvalue < 0.01 & DEGs_M1vsF1$log2FoldChange < -1])
DEGs_M2vsF2 <- read_excel("M2vsF2_deg.xlsx")
UP_DEGs_M2vsF2 <- as.character(DEGs_M2vsF2$gene.id[DEGs_M2vsF2$pvalue < 0.01 & DEGs_M2vsF2$log2FoldChange > 1])
DOWN_DEGs_M2vsF2 <- as.character(DEGs_M2vsF2$gene.id[DEGs_M2vsF2$pvalue < 0.01 & DEGs_M2vsF2$log2FoldChange < -1])
SDR_genes <- read_excel("SDR.xlsx", col_names = FALSE)[[1]]

results <- list()
for (mod in sigModules) {
  cat("\n=== Module", mod, "===\n")
  modGenes <- moduleColors == substring(mod, 3)
  moduleGenes <- names(datExpr)[modGenes]

  overlap_UP_M1 <- length(intersect(moduleGenes, UP_DEGs_M1vsF1))
  overlap_DOWN_M1 <- length(intersect(moduleGenes, DOWN_DEGs_M1vsF1))
  overlap_UP_M2 <- length(intersect(moduleGenes, UP_DEGs_M2vsF2))
  overlap_DOWN_M2 <- length(intersect(moduleGenes, DOWN_DEGs_M2vsF2))
  overlap_SDR <- length(intersect(moduleGenes, SDR_genes))

  results[[mod]] <- data.frame(
    Module = mod,
    GeneCount = length(moduleGenes),
    UP_M1vsF1 = overlap_UP_M1,
    DOWN_M1vsF1 = overlap_DOWN_M1,
    UP_M2vsF2 = overlap_UP_M2,
    DOWN_M2vsF2 = overlap_DOWN_M2,
    SDR_overlap = overlap_SDR,
    stringsAsFactors = FALSE
  )
}

final_results <- do.call(rbind, results)
write.xlsx(final_results, file = "All_Modules_Overlap_Results.xlsx", rowNames = FALSE)

# --------- step 14. Functional enrichment analysis of significant modules ---------
# Load required packages
library(clusterProfiler)
library(dplyr)
library(stringr)
library(jsonlite)
library(stringr)
library(devtools)
library(KEGGREST)
BiocManager::install("GO.db")
library(GO.db)

# KEGG annotation preparation
pathway2name <- tibble(Pathway = character(), Name = character())
ko2pathway <- tibble(KO = character(), Pathway = character())
kegg <- fromJSON("ko00001.json")

for (a in seq_along(kegg[["children"]][["children"]])) {
  for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
    for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
      pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
      pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
      pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
      pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
      kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
      kos <- str_match(kos_info, "K[0-9]*")[,1]
      ko2pathway <- rbind(ko2pathway, tibble(KO = kos, Pathway = rep(pathway_id, length(kos))))
    }
  }
}
colnames(ko2pathway) <- c("KO","Pathway")
write.table(pathway2name,"KEGG.library",sep="\t",row.names = F)

options(stringsAsFactors = F)
egg <- read.delim("female.emapper.annotations",header = T,sep="\t")
egg[egg==""]<-NA
gene2ko <- egg %>%
  dplyr::select(GID = query, KO = KEGG_ko) %>%
  na.omit()
pathway2name <- read.delim("KEGG.library")
colnames(pathway2name)<-c("Pathway","Name")

ko2gene <- tibble(Ko=character(),GID=character())
for (Query in gene2ko$GID){
  ko_list <- strsplit(gene2ko$KO[which(gene2ko[,1]==Query)],split = ',')
  for (ko in ko_list){
    if (length(which(ko2gene[,1]==ko))==0){
      tmp <- data.frame(Ko=ko,GID=Query)
      ko2gene <- rbind(ko2gene,tmp)
    }
    else{
      old_Query <- ko2gene$GID[which(ko2gene[,1]==ko)]
      ko2gene$GID[which(ko2gene[,1]==ko)] <- paste(old_Query,Query,sep = ',')
    }
  }
}
pathway2gene <- tibble(Pathway = character(), GID = character())

for (ko in ko2pathway$KO){
  pathway_list <- ko2pathway$Pathway[which(ko2pathway[,1]==ko)]
  for (pathway in pathway_list){
    if (paste('ko:',ko,sep='') %in% ko2gene$Ko){
      ko <- paste('ko:',ko,sep='')
      if (length(which(pathway2gene[,1]==pathway))==0 ){
        ko2gene$GID[which(ko2gene[,1]==ko)]
        tmp <- data.frame(pathway=pathway,GID=ko2gene$GID[which(ko2gene[,1]==ko)])
        pathway2gene <- rbind(pathway2gene,tmp)
      }
      else{
        old_Query <- pathway2gene$GID[which(pathway2gene[,1]==pathway)]
        Query <- ko2gene$GID[which(ko2gene[,1]==ko)]
        pathway2gene$GID[which(pathway2gene[,1]==pathway)] <- paste(old_Query,Query,sep=',')
      }
    }
  }
}
new_pathway2gene <- data.frame()	

for (i in 1:nrow(pathway2gene)) {
  pathway <- pathway2gene$pathway[i]
  genes <- strsplit(as.character(pathway2gene$GID[i]), ",")[[1]]
  for (gene in genes) {
    new_row <- data.frame(pathway = pathway, gene = gene)
    new_pathway2gene <- rbind(new_pathway2gene, new_row)
  }
}

# GO annotation preparation
egg_go <- read.delim("female.emapper.annotations", header = TRUE, sep = "\t")
egg_go[egg_go == ""] <- NA
gene_ids <- egg_go$query
eggnog_lines_with_go <- egg_go$GOs != ""
eggnog_annotations_go <- str_split(egg_go[eggnog_lines_with_go, ]$GOs, ",")
gene2go <- data.frame(
  gene = rep(gene_ids[eggnog_lines_with_go], times = sapply(eggnog_annotations_go, length)),
  term = unlist(eggnog_annotations_go)
)
go2name <- read.delim("GO.library", header = FALSE, stringsAsFactors = FALSE)
names(go2name) <- c("ID", "Description", "Ontology")

# GO enrichment analysis
gene_select <- read.delim(file = "HubGenes_ME1.txt", stringsAsFactors = FALSE, header = FALSE)$V1
go_rich <- enricher(gene = gene_select,
                    TERM2GENE = gene2go[c('term','gene')], 
                    TERM2NAME = go2name[c('ID', 'Description')], 
                    pvalueCutoff = 1, 
                    pAdjustMethod = 'BH', 
                    qvalueCutoff = 1)
tmp <- merge(go_rich, go2name[c('ID', 'Ontology')], by = 'ID')
tmp <- tmp[order(tmp$pvalue), ]
write.xlsx(tmp, file = 'ME1.GO_enrichment.xlsx', sep = '\t', rowNames = FALSE)
# KEGG enrichment analysis
kegg_rich <- enricher (gene = gene_select,
                       TERM2GENE = new_pathway2gene[c ('pathway','gene')], 
                       TERM2NAME = pathway2name[c ('Pathway','Name')], 
                       pvalueCutoff = 1, 
                       pAdjustMethod = 'BH', 
                       qvalueCutoff = 1, 
                       maxGSSize = 500)
write.xlsx(kegg_rich, file = 'ME1.KEGG_enrichment.xlsx', sep = '\t', rowNames = FALSE)

save.image(file = "WGCNA.RData")
