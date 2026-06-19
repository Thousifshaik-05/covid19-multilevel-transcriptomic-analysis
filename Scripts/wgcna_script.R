setwd("C://B.TECH//covid19-gene-expression-analysis//WGCNA")

library("WGCNA")
library(readxl)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

norm_counts <- read.csv("normalized_counts.csv")

# Check structure
dim(norm_counts)
norm_counts[1:3, 1:5]

trait_data <- read.csv("trait_data.csv")

dim(trait_data)
head(trait_data)

#1.EXPRESSION MATRIX
# Set gene names as rownames and remove the ID column
rownames(norm_counts) <- norm_counts$X
norm_counts <- norm_counts[, -1]  # remove first column

dim(norm_counts)  # should be 36148 x 34

#2.FILTER LOW VARIANCE GENES(KEEP TOP 5000)
library(WGCNA)

# Calculate variance for each gene
gene_var <- apply(norm_counts, 1, var)

# Keep top 5000 most variable genes
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:5000]
expr_filtered <- norm_counts[top_genes, ]

dim(expr_filtered)  # should be 5000 x 34

#3.TRANSPOSE
datExpr <- as.data.frame(t(expr_filtered))

dim(datExpr)  # should be 34 x 5000

# Check for bad genes/samples
gsg <- goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK  # should say TRUE


#4.PICK SOFT THRESHOLD POWER

powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit",
     type = "n",
     main = "Scale independence")
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, col = "red")
abline(h = 0.90, col = "red")


#5.BUILD NETWORK
png("Soft_Threshold_Plot.png", width = 800, height = 600)

plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit",
     type = "n",
     main = "Scale independence")
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, col = "red")
abline(h = 0.90, col = "red")

dev.off()
moduleColors <- net$colors
table(moduleColors)


#6.MODULE TRAIT CORRELATION
# Get module eigengenes
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs)

# Prepare trait
rownames(trait_data) <- trait_data$Sample
trait_data <- trait_data[rownames(datExpr), ]  # match sample order
trait_vector <- as.data.frame(trait_data$Condition)
rownames(trait_vector) <- rownames(datExpr)
colnames(trait_vector) <- "Condition"

# Correlation
moduleTraitCor  <- cor(MEs, trait_vector, use = "p")
moduleTraitPval <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# View results
print(moduleTraitCor)
print(moduleTraitPval)



#7.HEATMAP
png("Module_Trait_Heatmap.png", width = 800, height = 1000)

par(mar = c(6, 10, 3, 3))
labeledHeatmap(
  Matrix        = moduleTraitCor,
  xLabels       = "Condition",
  yLabels       = names(MEs),
  ySymbols      = names(MEs),
  colorLabels   = FALSE,
  colors        = blueWhiteRed(50),
  textMatrix    = textMatrix,
  setStdMargins = FALSE,
  cex.text      = 0.8,
  zlim          = c(-1, 1),
  main          = "Module-Trait Relationships"
)

dev.off()


8.#EXTRACT HUB GENES-TORQUOISE
module_of_interest <- "turquoise"

# All genes in this module
module_genes <- names(moduleColors)[moduleColors == module_of_interest]
cat("Genes in turquoise module:", length(module_genes), "\n")

# Gene Significance
GS <- as.numeric(cor(datExpr, trait_vector$Condition, use = "p"))
names(GS) <- colnames(datExpr)

# Module Membership
MM <- as.numeric(cor(datExpr,
                     MEs[, "MEturquoise"],
                     use = "p"))
names(MM) <- colnames(datExpr)

# Hub genes: MM > 0.8 and GS > 0.2
hub_genes_wgcna <- names(MM)[
  names(MM) %in% module_genes &
    abs(MM) > 0.8 &
    abs(GS) > 0.2
]

cat("Hub genes found:", length(hub_genes_wgcna), "\n")
print(hub_genes_wgcna)




#INTERSECTING WITH DEG AND FORMING PROTEIN LIST
#1.CONVERT ENSEMBLE IDS TO GENE SYMBOLS
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(AnnotationDbi)

# Convert WGCNA hub genes
wgcna_symbols <- mapIds(org.Hs.eg.db,
                        keys      = hub_genes_wgcna,
                        column    = "SYMBOL",
                        keytype   = "ENSEMBL",
                        multiVals = "first")

wgcna_symbols <- na.omit(wgcna_symbols)
wgcna_symbols <- unique(wgcna_symbols)
cat("WGCNA hub genes converted:", length(wgcna_symbols), "\n")


#2.CONVERT DEG ENSEMBLE IDS TOO

# Load DEG results
all_deg <- read.csv("C:/B.TECH/covid19-gene-expression-analysis/Data/All_DEG_results.csv")

# Filter significant DEGs
sig_deg <- all_deg[all_deg$padj < 0.05 & abs(all_deg$log2FoldChange) >= 1, ]
sig_deg <- na.omit(sig_deg)

cat("Significant DEGs:", nrow(sig_deg), "\n")

# Your sig_deg$gene contains Ensembl IDs
deg_symbols <- mapIds(org.Hs.eg.db,
                      keys      = sig_deg$gene,
                      column    = "SYMBOL",
                      keytype   = "ENSEMBL",
                      multiVals = "first")

deg_symbols <- na.omit(deg_symbols)
deg_symbols <- unique(deg_symbols)
cat("DEG genes converted:", length(deg_symbols), "\n")

#3.INTERSECT
final_protein_list <- intersect(deg_symbols, wgcna_symbols)
cat("Final protein list for STRING:", length(final_protein_list), "\n")
print(final_protein_list)


#4.PROTEIN FINAL LIST
write.table(final_protein_list,
            "C:/B.TECH/covid19-gene-expression-analysis/PIP_ANALYSIS/protein_list.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote     = FALSE)

cat("✓ Saved! Ready for STRING\n")




#GO ANALYSIS
# Install if needed
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
install.packages("ggplot2")

# Load libraries
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)

# Your hub genes
hub_genes_final <- c("CDC45", "SPAG5", "TTK", "CCNA2", "BIRC5",
                     "BUB1B", "PBK", "CDK1", "KIF2C", "MELK")

setwd("C://B.TECH//covid19-gene-expression-analysis//PIP_ANALYSIS")

#GO PLOT
go_results <- enrichGO(
  gene          = hub_genes_final,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# Check results
head(go_results)

go_plot <- dotplot(go_results, showCategory = 15, title = "GO Enrichment - Hub Genes")

ggsave("GO_Enrichment_HubGenes.png",
       plot   = go_plot,
       width  = 10,
       height = 8,
       dpi    = 300)

cat("✓ GO plot saved\n")


#KEGG ANALYSIS
# Convert to Entrez IDs
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys    = hub_genes_final,
                     column  = "ENTREZID",
                     keytype = "SYMBOL")
entrez_ids <- na.omit(entrez_ids)

# KEGG Analysis
kegg_results <- enrichKEGG(
  gene         = entrez_ids,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

# Save plot
kegg_plot <- dotplot(kegg_results, showCategory = 15, title = "KEGG Enrichment - Hub Genes")

ggsave("KEGG_Enrichment_HubGenes.png",
       plot   = kegg_plot,
       width  = 10,
       height = 8,
       dpi    = 300)

cat("✓ KEGG plot saved\n")